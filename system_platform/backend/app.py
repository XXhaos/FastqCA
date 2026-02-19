import os
import shutil
from datetime import datetime

from flask import Flask, g, jsonify, request, send_file
from flask_cors import CORS
from werkzeug.utils import secure_filename

from auth import create_token, login_required
from config import Config
from extensions import db
from models import FastqFile, Task, User
from tasks import compress_task


def ensure_dirs(app):
    for key in ["STORAGE_ROOT", "CHUNK_DIR", "UPLOAD_DIR", "OUTPUT_DIR", "REPORT_DIR"]:
        os.makedirs(app.config[key], exist_ok=True)


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
        if User.query.filter_by(username=data.get("username", "")).first():
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
        return {"token": create_token(user.id), "role": user.role}

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

    @app.get("/api/tasks/<int:task_id>")
    @login_required
    def get_task(task_id: int):
        task = Task.query.filter_by(id=task_id, user_id=g.current_user.id).first_or_404()
        return {
            "id": task.id,
            "status": task.status,
            "progress": task.progress,
            "quality_mode": task.quality_mode,
            "result_path": task.result_path,
            "report_path": task.report_path,
            "error_message": task.error_message,
            "file_id": task.file_id,
        }

    @app.get("/api/tasks")
    @login_required
    def list_tasks():
        rows = Task.query.filter_by(user_id=g.current_user.id).order_by(Task.created_at.desc()).all()
        return jsonify([
            {
                "id": t.id,
                "file_id": t.file_id,
                "quality_mode": t.quality_mode,
                "threads": t.threads,
                "status": t.status,
                "progress": t.progress,
                "started_at": t.started_at.isoformat() if t.started_at else None,
                "finished_at": t.finished_at.isoformat() if t.finished_at else None,
                "error_message": t.error_message,
            }
            for t in rows
        ])

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

    @app.get("/api/analytics/<int:task_id>")
    @login_required
    def analytics(task_id: int):
        task = Task.query.filter_by(id=task_id, user_id=g.current_user.id).first_or_404()
        if task.status != "finished":
            return {"error": "task not finished"}, 400

        source = task.source_file
        return {
            "bar": {
                "labels": ["Original", "FastqCA", "Gzip"],
                "values": [task.original_size or source.file_size, task.fastqca_size or 0, task.gzip_size or 0],
            },
            "quality": {
                "x": [10, 20, 30, 40],
                "y": [1200, 6500, 9800, 4300],
            },
        }

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(host="0.0.0.0", port=5000, debug=True)
