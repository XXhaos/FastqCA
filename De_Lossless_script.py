import os
import subprocess
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
import csv

import psutil

# 与 LossLess_script_new 保持一致的监控逻辑


def get_file_size(file_path: Path) -> int:
    return file_path.stat().st_size if file_path.exists() else 0


def monitor_process(process: subprocess.Popen) -> dict:
    cpu_percentages = []
    memory_usages = []
    max_total_memory = 0
    try:
        p = psutil.Process(process.pid)
        step_count = 0
        while process.poll() is None:
            try:
                all_procs = [p] + p.children(recursive=True)
                total_rss = 0
                total_cpu = 0
                mem_breakdown = defaultdict(float)
                for proc in all_procs:
                    try:
                        mem_info = proc.memory_info()
                        name = proc.name()
                        rss_mb = mem_info.rss / 1024 / 1024
                        total_rss += rss_mb
                        total_cpu += proc.cpu_percent()
                        if "python" in name.lower():
                            group = "Python"
                        elif "lpaq8" in name.lower():
                            group = "lpaq8"
                        else:
                            group = "Other"
                        mem_breakdown[group] += rss_mb
                    except Exception:
                        continue

                memory_usages.append(total_rss)
                cpu_percentages.append(total_cpu)
                max_total_memory = max(max_total_memory, total_rss)

                if step_count % 20 == 0 and len(all_procs) > 1:
                    breakdown_str = " | ".join([f"{k}: {v:.1f}MB" for k, v in mem_breakdown.items()])
                    print(f"[监控] {breakdown_str} (总: {total_rss:.1f}MB)")

                step_count += 1
            except Exception:
                break
            time.sleep(0.05)
    except Exception:
        pass

    return {
        "avg_cpu": sum(cpu_percentages) / len(cpu_percentages) if cpu_percentages else 0.0,
        "max_memory": max_total_memory,
    }


def decompress_and_collect_metrics(input_file: Path, output_dir: Path, decompressed_dir: Path) -> dict:
    file_name = input_file.name
    decompressed_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python",
        "main_new.py",
        "--compressor",
        "Lossless",
        "--input_path",
        str(input_file),
        "--output_path",
        str(decompressed_dir),
        "--mode",
        "d",
    ]

    start_time = time.time()
    process = subprocess.Popen(cmd)
    metrics = monitor_process(process)
    process.wait()
    end_time = time.time()
    decompress_time = end_time - start_time

    time.sleep(1)

    compressed_size = get_file_size(input_file)
    in_stem = input_file.stem
    candidate_fastq = decompressed_dir / f"{in_stem}.fastq"
    if not candidate_fastq.exists():
        fastq_candidates = sorted(
            [p for p in decompressed_dir.glob("*.fastq") if p.is_file()],
            key=lambda p: p.stat().st_size,
            reverse=True,
        )
        if fastq_candidates:
            candidate_fastq = fastq_candidates[0]

    output_size = get_file_size(candidate_fastq)
    decompress_speed = (output_size / 1024 / 1024) / decompress_time if decompress_time > 0 else 0.0

    return {
        "file_name": file_name,
        "compressed_size_mb": compressed_size / 1024 / 1024,
        "decompressed_size_mb": output_size / 1024 / 1024,
        "decompress_time_s": decompress_time,
        "decompress_speed_mbs": decompress_speed,
        "avg_cpu_percent": metrics["avg_cpu"],
        "max_memory_mb": metrics["max_memory"],
    }


def main() -> None:
    input_dir = "/media/compress/新加卷/output/test_FastCA_thread"
    decompressed_dir = "/media/compress/新加卷/output/test_FastCA_thread_decompressed"
    Path(decompressed_dir).mkdir(parents=True, exist_ok=True)

    csv_path = os.path.join(
        input_dir, f"decompression_metrics_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    )
    csv_header = [
        "file_name",
        "compressed_size_mb",
        "decompressed_size_mb",
        "decompress_time_s",
        "decompress_speed_mbs",
        "avg_cpu_percent",
        "max_memory_mb",
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(csv_header)

    candidates = [p for p in Path(input_dir).iterdir() if p.is_file() and p.suffix.lower() != ".csv"]
    if not candidates:
        print(f"No compressed files found in {input_dir}.")
        return

    for comp_path in candidates:
        print(f"Decompressing {comp_path} ...")
        metrics = decompress_and_collect_metrics(comp_path, Path(input_dir), Path(decompressed_dir))
        with open(csv_path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([metrics[k] for k in csv_header])
        print(f"Completed: {comp_path}")

    print(f"\nAll done. Metrics saved to: {csv_path}")
    print(f"Decompressed FASTQ saved under: {decompressed_dir}")


if __name__ == "__main__":
    main()
