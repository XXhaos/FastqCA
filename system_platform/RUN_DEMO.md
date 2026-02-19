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

## 5. 页面演示流程
1. 在页面顶部输入用户名密码，先“注册”再“登录”。
2. 选择 FASTQ 文件。
3. 选择压缩模式（无损/有损）与线程数。
4. 点击“上传并开始压缩”。
5. 观察任务进度条自动更新。
6. 任务完成后点击：
   - “查看分析”：加载 Original/FastqCA/Gzip 柱状图；
   - “下载压缩包”；
   - “下载报告”。

## 6. 注意事项
- 如果你之前运行过旧版本并且数据库结构冲突，删除 `system_platform/backend/instance/fastqca_platform.db` 后重启后端。
- FastqCA 实际压缩由 `main_new.py` 调用，需确保仓库根目录已有 `lpaq8` 可执行文件并可运行。
