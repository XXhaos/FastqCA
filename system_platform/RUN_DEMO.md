# FastqCA Web 系统演示（Ubuntu）

## 1. 安装依赖
```bash
cd /workspace/FastqCA
python3 -m venv .venv
source .venv/bin/activate
pip install -r system_platform/backend/requirements.txt
```

## 2. 启动 Redis
```bash
sudo apt-get update
sudo apt-get install -y redis-server
sudo systemctl enable --now redis-server
redis-cli ping
```

## 3. 启动后端 + Worker
终端A：
```bash
source /workspace/FastqCA/.venv/bin/activate
cd /workspace/FastqCA/system_platform/backend
python app.py
```

终端B：
```bash
source /workspace/FastqCA/.venv/bin/activate
cd /workspace/FastqCA/system_platform/backend
celery -A tasks.celery worker --loglevel=info
```

## 4. 启动前端
终端C：
```bash
cd /workspace/FastqCA/system_platform/frontend
python3 -m http.server 4173
```
浏览器打开：`http://127.0.0.1:4173`

## 5. 模块演示建议（左侧菜单切换四个页面）
1. **压缩任务页面**：上传 FASTQ，选无损/有损模式，启动任务并观察进度与压缩效果分析图。
2. **压缩历史页面**：展示压缩率、压缩耗时、吞吐率历史记录。
3. **性能监测页面**：展示 CPU/内存曲线和任务成功率等指标（示例数据内存占用约 10%~20%）。
4. **用户管理页面**（admin）：查看用户列表并修改角色。

> 说明：即便后端未连接成功，前端也会显示完整演示数据（已按真实压缩结果填充），便于论文答辩展示界面与流程。

## 6. 注意事项
- 若数据库是旧版本结构，删除 `system_platform/backend/instance/fastqca_platform.db` 后重启。
- FastqCA 压缩调用 `main_new.py`，请确保仓库根目录 `lpaq8` 可执行。
