import os


class Config:
    SECRET_KEY = os.getenv("SECRET_KEY", "fastqca-dev-secret")
    SQLALCHEMY_DATABASE_URI = os.getenv("DATABASE_URL", "sqlite:///fastqca_platform.db")
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    STORAGE_ROOT = os.getenv("STORAGE_ROOT", os.path.abspath("storage"))
    CHUNK_DIR = os.path.join(STORAGE_ROOT, "chunks")
    UPLOAD_DIR = os.path.join(STORAGE_ROOT, "uploads")
    OUTPUT_DIR = os.path.join(STORAGE_ROOT, "outputs")
    REPORT_DIR = os.path.join(STORAGE_ROOT, "reports")

    CELERY_BROKER_URL = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0")
    CELERY_RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")
