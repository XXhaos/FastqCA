import os
import subprocess
import time
import psutil
import csv
from pathlib import Path
from datetime import datetime

# 解压实验线程数（需与压缩侧保持一致）
THREAD_COUNT = "4"

def get_file_size(file_path):
    """返回文件大小（字节）"""
    return os.path.getsize(file_path)

def monitor_process(process):
    """监控进程CPU与内存（记录平均CPU与最大RSS内存MB）"""
    cpu_percentages = []
    memory_usages = []
    try:
        p = psutil.Process(process.pid)
        # 先调用一次，避免第一次调用总是0
        p.cpu_percent(interval=None)
        while process.poll() is None:
            cpu_percentages.append(p.cpu_percent(interval=0.1))
            memory_usages.append(p.memory_info().rss / 1024 / 1024)
        max_memory = max(memory_usages) if memory_usages else 0.0
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        max_memory = 0.0

    return {
        "avg_cpu": sum(cpu_percentages) / len(cpu_percentages) if cpu_percentages else 0.0,
        "max_memory": max_memory
    }

def decompress_and_collect_metrics(input_file, output_dir, decompressed_dir):
    """
    解压缩单个文件，并收集指标：
    - 解压前(压缩文件)大小
    - 解压后(还原fastq)大小
    - 解压用时(s)
    - 解压速度(MB/s)：以“解压后输出字节数/时间”计算
    - 平均CPU占用(%)
    - 最大RSS内存(MB)
    """
    # file_name 不带扩展名（压缩脚本将输出文件命名为 output_dir/<原文件名>）
    file_name = os.path.basename(input_file)
    # 解压输出到 decompressed_dir
    os.makedirs(decompressed_dir, exist_ok=True)

    # 对于 main.py 的接口：--input_path 指向压缩文件；--output_path 指向解压输出目录
    cmd = [
        "python", "main_new.py",
        "--compressor", "Lossy",
        "--input_path", input_file,
        "--output_path", decompressed_dir,
        "--mode", "d",
        "--threads", THREAD_COUNT,
    ]

    start_time = time.time()
    process = subprocess.Popen(cmd)
    metrics = monitor_process(process)
    process.wait()
    end_time = time.time()
    decompress_time = end_time - start_time

    # 等待IO落盘
    time.sleep(1)

    # 统计大小
    compressed_size = get_file_size(input_file) if os.path.exists(input_file) else 0

    # 根据压缩脚本的约定：压缩时以原“<name>.fastq”去掉后缀得到输出文件名 <name>
    # 因此我们在解压后尝试在 decompressed_dir 中寻找 <name>.fastq
    # 1) 若 input_file 名为 <name>（无扩展），则目标推断为 <name>.fastq
    # 2) 若 input_file 自身有扩展名，则去掉扩展得到 <stem>.fastq
    in_stem = Path(file_name).stem
    candidate_fastq = Path(decompressed_dir) / f"{in_stem}.fastq"

    # 兜底：若找不到，尝试直接统计 decompressed_dir 下最大文件，作为输出fastq
    if not candidate_fastq.exists():
        fastq_candidates = sorted(
            [p for p in Path(decompressed_dir).glob("*.fastq") if p.is_file()],
            key=lambda p: p.stat().st_size,
            reverse=True
        )
        if fastq_candidates:
            candidate_fastq = fastq_candidates[0]

    if candidate_fastq.exists():
        output_size = get_file_size(str(candidate_fastq))
    else:
        print(f"Warning: Decompressed FASTQ not found for {input_file}")
        output_size = 0

    # 解压速度：以“解压后输出字节数/时间”计（更符合解压透出速率）
    decompress_speed = (output_size / 1024 / 1024) / decompress_time if decompress_time > 0 else 0.0

    return {
        "file_name": file_name,
        "compressed_size_mb": compressed_size / 1024 / 1024,
        "decompressed_size_mb": output_size / 1024 / 1024,
        "decompress_time_s": decompress_time,
        "decompress_speed_mbs": decompress_speed,
        "avg_cpu_percent": metrics["avg_cpu"],
        "max_memory_mb": metrics["max_memory"]
    }

def main():
    # 与压缩脚本一致：压缩输出目录
    input_dir = "/media/compress/新加卷/output/test_FastCA-Lossy"     # 解压的输入目录（= 压缩脚本的输出目录）
    decompressed_dir = "/media/compress/新加卷/output/test_FastCA-Lossy_decompressed"  # 解压后的FASTQ目录
    Path(input_dir).mkdir(parents=True, exist_ok=True)
    os.makedirs(decompressed_dir, exist_ok=True)

    # 结果CSV
    csv_path = os.path.join(
        input_dir,
        f"decompression_metrics_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    )
    Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
    csv_header = [
        "file_name",
        "compressed_size_mb",
        "decompressed_size_mb",
        "decompress_time_s",
        "decompress_speed_mbs",
        "avg_cpu_percent",
        "max_memory_mb"
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(csv_header)

    # 选择需要解压的“压缩产物”
    # 压缩脚本的产物命名为 output_dir/<原文件名去掉.fastq>（无扩展名），同目录还会有metrics CSV。
    # 这里过滤：排除目录与 .csv，仅处理常规文件。
    candidates = [
        p for p in Path(input_dir).iterdir()
        if p.is_file() and p.suffix.lower() != ".csv"
    ]
    # 若你的 main.py 产物有专用扩展名，可改成类似：Path(input_dir).glob("*.lz") 等。

    if not candidates:
        print(f"No compressed files found in {input_dir}.")
        return

    for comp_path in candidates:
        print(f"Decompressing {comp_path} ...")
        metrics = decompress_and_collect_metrics(str(comp_path), input_dir, decompressed_dir)
        with open(csv_path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([metrics[k] for k in csv_header])
        print(f"Completed: {comp_path}")

    print(f"\nAll done. Metrics saved to: {csv_path}")
    print(f"Decompressed FASTQ saved under: {decompressed_dir}")

if __name__ == "__main__":
    main()
