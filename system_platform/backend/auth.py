from functools import wraps

from flask import current_app, g, request
from itsdangerous import BadSignature, URLSafeTimedSerializer

from models import User


def _serializer() -> URLSafeTimedSerializer:
    return URLSafeTimedSerializer(current_app.config["SECRET_KEY"])


def create_token(user_id: int) -> str:
    return _serializer().dumps({"uid": user_id})


def parse_token(token: str, max_age: int = 7 * 24 * 3600):
    try:
        payload = _serializer().loads(token, max_age=max_age)
        return payload.get("uid")
    except BadSignature:
        return None


def login_required(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        token = request.headers.get("Authorization", "").replace("Bearer ", "")
        if not token:
            return {"error": "missing token"}, 401
        user_id = parse_token(token)
        if not user_id:
            return {"error": "invalid token"}, 401
        user = User.query.get(user_id)
        if not user:
            return {"error": "user not found"}, 401
        g.current_user = user
        return fn(*args, **kwargs)

    return wrapper
