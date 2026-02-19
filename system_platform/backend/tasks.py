import os
from datetime import datetime

from celery import Celery
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

from config import Config
from extensions import db
from fastqca_wrapper import quality_distribution, run_fastqca, run_gzip, size_stats
from models import Task


celery = Celery(__name__, broker=Config.CELERY_BROKER_URL, backend=Config.CELERY_RESULT_BACKEND)


def generate_pdf_report(task: Task, stats: dict, qdist: dict, report_path: str):
    c = canvas.Canvas(report_path, pagesize=A4)
    c.setFont("Helvetica", 12)
    c.drawString(60, 800, f"Task ID: {task.id}")
    c.drawString(60, 780, f"Quality Mode: {task.quality_mode}")
    c.drawString(60, 760, f"Original Size: {stats['original']} bytes")
    c.drawString(60, 740, f"FastqCA Size: {stats['fastqca']} bytes")
    c.drawString(60, 720, f"Gzip Size: {stats['gzip']} bytes")
    c.drawString(60, 700, f"FastqCA Ratio: {stats['fastqca_ratio']:.2%}")
    c.drawString(60, 680, f"Gzip Ratio: {stats['gzip_ratio']:.2%}")
    c.drawString(60, 650, "Quality Score Distribution (Top 10 bins)")
    for idx, q in enumerate(qdist["x"][:10]):
        c.drawString(80, 630 - idx * 18, f"Q{q}: {qdist['y'][idx]}")
    c.save()


@celery.task(bind=True)
def compress_task(self, task_id: int):
    from app import create_app

    app = create_app()
    with app.app_context():
        task = Task.query.get(task_id)
        if not task:
            return

        task.status = "running"
        task.started_at = datetime.utcnow()
        task.progress = 10
        db.session.commit()

        source = task.source_file
        input_path = source.filepath
        fastqca_output = os.path.join(Config.OUTPUT_DIR, f"task_{task.id}.fqc")
        gzip_output = os.path.join(Config.OUTPUT_DIR, f"task_{task.id}.fastq.gz")
        report_path = os.path.join(Config.REPORT_DIR, f"task_{task.id}.pdf")

        try:
            run_fastqca(input_path, fastqca_output, task.quality_mode, task.threads)
            task.progress = 60
            db.session.commit()

            run_gzip(input_path, gzip_output)
            task.progress = 80
            db.session.commit()

            stats = size_stats(input_path, fastqca_output, gzip_output)
            qdist = quality_distribution(input_path)
            generate_pdf_report(task, stats, qdist, report_path)

            source.compressed_size = stats["fastqca"]
            source.compression_ratio = stats["fastqca_ratio"]
            source.status = "compressed"

            task.result_path = fastqca_output
            task.report_path = report_path
            task.original_size = stats["original"]
            task.fastqca_size = stats["fastqca"]
            task.gzip_size = stats["gzip"]
            task.progress = 100
            task.status = "finished"
            task.finished_at = datetime.utcnow()
            task.error_message = None
        except Exception as exc:
            task.status = "failed"
            task.error_message = str(exc)
            task.finished_at = datetime.utcnow()

        db.session.commit()
