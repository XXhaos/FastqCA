#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FRONTEND_DIR="$ROOT_DIR/frontend"
OUT_DIR="${1:-$ROOT_DIR/screenshots}"
PORT="${PORT:-4173}"

mkdir -p "$OUT_DIR"

if ! command -v python3 >/dev/null 2>&1; then
  echo "[ERROR] python3 未安装"
  exit 1
fi

if ! python3 -c "import playwright" >/dev/null 2>&1; then
  echo "[INFO] 正在安装 Python Playwright..."
  python3 -m pip install --user playwright
fi

if ! python3 -c "import playwright" >/dev/null 2>&1; then
  echo "[ERROR] Playwright 安装失败"
  exit 1
fi

echo "[INFO] 安装浏览器驱动（仅首次较慢）..."
python3 -m playwright install chromium >/dev/null

cleanup() {
  if [[ -n "${HTTP_PID:-}" ]] && kill -0 "$HTTP_PID" >/dev/null 2>&1; then
    kill "$HTTP_PID" >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT

echo "[INFO] 启动前端静态服务: http://127.0.0.1:${PORT}"
python3 -m http.server "$PORT" --directory "$FRONTEND_DIR" >/tmp/fastqca_screenshot_http.log 2>&1 &
HTTP_PID=$!
sleep 1

if ! curl -sf "http://127.0.0.1:${PORT}" >/dev/null; then
  echo "[ERROR] 前端服务启动失败，请检查 /tmp/fastqca_screenshot_http.log"
  exit 1
fi

OUT_DIR="$OUT_DIR" PORT="$PORT" python3 - <<'PY'
import os
from playwright.sync_api import sync_playwright

port = os.environ["PORT"]
out_dir = os.environ["OUT_DIR"]

shots = [
    ("压缩任务", "tasks.png"),
    ("压缩历史", "history.png"),
    ("性能监测", "performance.png"),
    ("用户管理", "users.png"),
]

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={"width": 1680, "height": 980})
    page.goto(f"http://127.0.0.1:{port}", wait_until="networkidle", timeout=90000)
    page.wait_for_timeout(1200)

    for text, filename in shots:
        page.get_by_text(text, exact=True).click()
        page.wait_for_timeout(800)
        page.screenshot(path=os.path.join(out_dir, filename), full_page=True)

    browser.close()

print("截图完成:")
for _, filename in shots:
    print(os.path.join(out_dir, filename))
PY

echo "[OK] 所有页面截图已生成到: $OUT_DIR"
