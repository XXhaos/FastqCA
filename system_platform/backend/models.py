from datetime import datetime

from werkzeug.security import check_password_hash, generate_password_hash

from extensions import db


class User(db.Model):
    __tablename__ = "users"

    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), unique=True, nullable=False)
    password_hash = db.Column(db.String(256), nullable=False)
    role = db.Column(db.String(16), default="researcher", nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    files = db.relationship("FastqFile", backref="owner", lazy=True)
    tasks = db.relationship("Task", backref="creator", lazy=True)

    def set_password(self, password: str) -> None:
        self.password_hash = generate_password_hash(password)

    def verify_password(self, password: str) -> bool:
        return check_password_hash(self.password_hash, password)


class FastqFile(db.Model):
    __tablename__ = "files"

    id = db.Column(db.Integer, primary_key=True)
    filename = db.Column(db.String(255), nullable=False)
    filepath = db.Column(db.String(512), nullable=False)
    file_size = db.Column(db.BigInteger, nullable=False)
    compressed_size = db.Column(db.BigInteger)
    compression_ratio = db.Column(db.Float)
    status = db.Column(db.String(32), default="uploaded")
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    user_id = db.Column(db.Integer, db.ForeignKey("users.id"), nullable=False)


class Task(db.Model):
    __tablename__ = "tasks"

    id = db.Column(db.Integer, primary_key=True)
    task_uuid = db.Column(db.String(64), unique=True)
    mode = db.Column(db.String(32), nullable=False)
    quality_mode = db.Column(db.String(16), default="lossless")
    threads = db.Column(db.Integer, default=4)
    status = db.Column(db.String(32), default="queued")
    progress = db.Column(db.Integer, default=0)
    result_path = db.Column(db.String(512))
    report_path = db.Column(db.String(512))
    error_message = db.Column(db.String(512))
    original_size = db.Column(db.BigInteger)
    fastqca_size = db.Column(db.BigInteger)
    gzip_size = db.Column(db.BigInteger)
    started_at = db.Column(db.DateTime)
    finished_at = db.Column(db.DateTime)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    user_id = db.Column(db.Integer, db.ForeignKey("users.id"), nullable=False)
    file_id = db.Column(db.Integer, db.ForeignKey("files.id"), nullable=False)

    source_file = db.relationship("FastqFile", backref="tasks", lazy=True)
