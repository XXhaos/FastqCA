import os
import shutil
from datetime import datetime

from flask import Flask, g, jsonify, request, send_file
from flask_cors import CORS
from werkzeug.utils import secure_filename

from auth import admin_required, create_token, login_required
from config import Config
from extensions import db
from models import FastqFile, PerformanceMetric, Task, User
from tasks import compress_task

try:
    import psutil
except Exception:  # pragma: no cover
    psutil = None


def ensure_dirs(app):
    for key in ["STORAGE_ROOT", "CHUNK_DIR", "UPLOAD_DIR", "OUTPUT_DIR", "REPORT_DIR"]:
        os.makedirs(app.config[key], exist_ok=True)


def task_payload(t: Task):
    return {
        "id": t.id,
        "file_id": t.file_id,
        "quality_mode": t.quality_mode,
        "threads": t.threads,
        "status": t.status,
        "progress": t.progress,
        "result_path": t.result_path,
        "report_path": t.report_path,
        "error_message": t.error_message,
        "compression_time_sec": t.compression_time_sec,
        "throughput_mb_s": t.throughput_mb_s,
        "created_at": t.created_at.isoformat() if t.created_at else None,
        "started_at": t.started_at.isoformat() if t.started_at else None,
        "finished_at": t.finished_at.isoformat() if t.finished_at else None,
    }


def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)
    db.init_app(app)
    CORS(app)
    ensure_dirs(app)

    with app.app_context():
        db.create_all()

    @app.post("/api/auth/register")
    def register():
        data = request.get_json() or {}
        if not data.get("username") or not data.get("password"):
            return {"error": "username and password required"}, 400
        if User.query.filter_by(username=data["username"]).first():
            return {"error": "username exists"}, 400
        user = User(username=data["username"], role=data.get("role", "researcher"))
        user.set_password(data["password"])
        db.session.add(user)
        db.session.commit()
        return {"message": "ok"}

    @app.post("/api/auth/login")
    def login():
        data = request.get_json() or {}
        user = User.query.filter_by(username=data.get("username", "")).first()
        if not user or not user.verify_password(data.get("password", "")):
            return {"error": "invalid credentials"}, 401
        return {"token": create_token(user.id), "role": user.role, "username": user.username}

    @app.get("/api/auth/me")
    @login_required
    def me():
        return {"id": g.current_user.id, "username": g.current_user.username, "role": g.current_user.role}

    @app.get("/api/users")
    @admin_required
    def list_users():
        users = User.query.order_by(User.created_at.desc()).all()
        return jsonify([
            {
                "id": u.id,
                "username": u.username,
                "role": u.role,
                "created_at": u.created_at.isoformat(),
            }
            for u in users
        ])

    @app.put("/api/users/<int:user_id>/role")
    @admin_required
    def change_user_role(user_id: int):
        data = request.get_json() or {}
        role = data.get("role")
        if role not in ["admin", "researcher"]:
            return {"error": "invalid role"}, 400
        user = User.query.get_or_404(user_id)
        user.role = role
        db.session.commit()
        return {"message": "updated"}

    @app.post("/api/files/chunk")
    @login_required
    def upload_chunk():
        upload_id = request.form["upload_id"]
        idx = request.form["chunk_index"]
        chunk = request.files["chunk"]

        folder = os.path.join(app.config["CHUNK_DIR"], upload_id)
        os.makedirs(folder, exist_ok=True)
        chunk.save(os.path.join(folder, f"{idx}.part"))
        return {"message": "chunk received", "chunk_index": int(idx)}

    @app.post("/api/files/merge")
    @login_required
    def merge_chunks():
        data = request.get_json() or {}
        upload_id = data["upload_id"]
        filename = secure_filename(data["filename"])
        total_chunks = int(data["total_chunks"])

        folder = os.path.join(app.config["CHUNK_DIR"], upload_id)
        out_path = os.path.join(app.config["UPLOAD_DIR"], f"{datetime.utcnow().timestamp()}_{filename}")

        with open(out_path, "wb") as out:
            for i in range(total_chunks):
                part_path = os.path.join(folder, f"{i}.part")
                if not os.path.exists(part_path):
                    return {"error": f"missing chunk {i}"}, 400
                with open(part_path, "rb") as part:
                    shutil.copyfileobj(part, out)

        shutil.rmtree(folder, ignore_errors=True)

        size = os.path.getsize(out_path)
        file_obj = FastqFile(filename=filename, filepath=out_path, file_size=size, user_id=g.current_user.id)
        db.session.add(file_obj)
        db.session.commit()
        return {"file_id": file_obj.id, "size": size, "filename": filename}

    @app.get("/api/files")
    @login_required
    def list_files():
        rows = FastqFile.query.filter_by(user_id=g.current_user.id).order_by(FastqFile.created_at.desc()).all()
        return jsonify([
            {
                "id": row.id,
                "filename": row.filename,
                "size": row.file_size,
                "compressed_size": row.compressed_size,
                "compression_ratio": row.compression_ratio,
                "status": row.status,
                "created_at": row.created_at.isoformat(),
            }
            for row in rows
        ])

    @app.post("/api/tasks")
    @login_required
    def create_task():
        data = request.get_json() or {}
        file_obj = FastqFile.query.filter_by(id=data["file_id"], user_id=g.current_user.id).first_or_404()
        task = Task(
            mode="compress",
            quality_mode=data.get("quality_mode", "lossless"),
            threads=int(data.get("threads", 4)),
            user_id=g.current_user.id,
            file_id=file_obj.id,
            status="queued",
        )
        db.session.add(task)
        db.session.commit()

        celery_task = compress_task.delay(task.id)
        task.task_uuid = celery_task.id
        db.session.commit()
        return {"task_id": task.id, "task_uuid": task.task_uuid}

    @app.get("/api/tasks")
    @login_required
    def list_tasks():
        rows = Task.query.filter_by(user_id=g.current_user.id).order_by(Task.created_at.desc()).all()
        return jsonify([task_payload(t) for t in rows])

    @app.get("/api/tasks/<int:task_id>")
    @login_required
    def get_task(task_id: int):
        task = Task.query.filter_by(id=task_id, user_id=g.current_user.id).first_or_404()
        return task_payload(task)

    @app.get("/api/tasks/<int:task_id>/download")
    @login_required
    def download_result(task_id: int):
        task = Task.query.filter_by(id=task_id, user_id=g.current_user.id).first_or_404()
        if not task.result_path or not os.path.exists(task.result_path):
            return {"error": "compressed file not ready"}, 404
        return send_file(task.result_path, as_attachment=True)

    @app.get("/api/reports/<int:task_id>")
    @login_required
    def download_report(task_id: int):
        task = Task.query.filter_by(id=task_id, user_id=g.current_user.id).first_or_404()
        if not task.report_path or not os.path.exists(task.report_path):
            return {"error": "report not ready"}, 404
        return send_file(task.report_path, as_attachment=True)

    @app.get("/api/history")
    @login_required
    def history():
        rows = Task.query.filter_by(user_id=g.current_user.id).filter(Task.status == "finished").order_by(Task.finished_at.desc()).all()
        return jsonify([
            {
                "task_id": t.id,
                "filename": t.source_file.filename,
                "quality_mode": t.quality_mode,
                "original_size": t.original_size,
                "fastqca_size": t.fastqca_size,
                "compression_ratio": (1 - (t.fastqca_size / t.original_size)) if t.original_size and t.fastqca_size else None,
                "compression_time_sec": t.compression_time_sec,
                "throughput_mb_s": t.throughput_mb_s,
                "finished_at": t.finished_at.isoformat() if t.finished_at else None,
            }
            for t in rows
        ])

    @app.get("/api/analytics/<int:task_id>")
    @login_required
    def analytics(task_id: int):
        task = Task.query.filter_by(id=task_id, user_id=g.current_user.id).first_or_404()
        if task.status != "finished":
            return {"error": "task not finished"}, 400
        return {
            "bar": {
                "labels": ["Original", "FastqCA"],
                "values": [task.original_size or 0, task.fastqca_size or 0],
            },
            "quality": {
                "x": [10, 20, 30, 40],
                "y": [1200, 6500, 9800, 4300],
            },
        }

    @app.get("/api/performance/system")
    @login_required
    def performance_system():
        cpu = psutil.cpu_percent(interval=0.1) if psutil else 0.0
        mem = psutil.virtual_memory().percent if psutil else 0.0
        active = Task.query.filter(Task.status == "running").count()
        queued = Task.query.filter(Task.status == "queued").count()

        metric = PerformanceMetric(cpu_percent=cpu, memory_percent=mem, active_tasks=active, queued_tasks=queued)
        db.session.add(metric)
        db.session.commit()

        return {
            "cpu_percent": cpu,
            "memory_percent": mem,
            "active_tasks": active,
            "queued_tasks": queued,
        }

    @app.get("/api/performance/overview")
    @login_required
    def performance_overview():
        finished = Task.query.filter(Task.status == "finished").all()
        failed = Task.query.filter(Task.status == "failed").count()
        all_count = Task.query.count()

        avg_time = (sum(t.compression_time_sec or 0 for t in finished) / len(finished)) if finished else 0
        avg_ratio = (
            sum((1 - (t.fastqca_size / t.original_size)) for t in finished if t.original_size and t.fastqca_size)
            / max(1, sum(1 for t in finished if t.original_size and t.fastqca_size))
        )
        success_rate = ((len(finished) / all_count) * 100) if all_count else 0

        recent_metrics = PerformanceMetric.query.order_by(PerformanceMetric.recorded_at.desc()).limit(20).all()
        return {
            "task_total": all_count,
            "task_finished": len(finished),
            "task_failed": failed,
            "avg_compression_time_sec": avg_time,
            "avg_compression_ratio": avg_ratio,
            "success_rate": success_rate,
            "system_recent": [
                {
                    "cpu_percent": m.cpu_percent,
                    "memory_percent": m.memory_percent,
                    "recorded_at": m.recorded_at.isoformat(),
                }
                for m in reversed(recent_metrics)
            ],
        }

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(host="0.0.0.0", port=5000, debug=True)
