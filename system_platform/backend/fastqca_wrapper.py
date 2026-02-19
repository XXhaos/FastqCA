import gzip
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO


REPO_ROOT = Path(__file__).resolve().parents[2]


def run_fastqca(input_path: str, output_path: str, quality_mode: str, threads: int):
    compressor = "LossLess" if quality_mode == "lossless" else "Lossy"
    cmd = [
        sys.executable,
        str(REPO_ROOT / "main_new.py"),
        "--compressor",
        compressor,
        "--input_path",
        input_path,
        "--output_path",
        output_path,
        "--mode",
        "compress",
        "--threads",
        str(threads),
    ]
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def run_gzip(input_path: str, output_path: str):
    with open(input_path, "rb") as src, gzip.open(output_path, "wb") as dst:
        dst.writelines(src)


def quality_distribution(fastq_path: str, sample_limit: int = 50000):
    counter = Counter()
    read_count = 0
    for record in SeqIO.parse(fastq_path, "fastq"):
        counter.update(record.letter_annotations["phred_quality"])
        read_count += 1
        if read_count >= sample_limit:
            break
    x = sorted(counter.keys())
    y = [counter[v] for v in x]
    return {"x": x, "y": y}


def size_stats(original: str, fastqca_out: str, gzip_out: str):
    original_size = os.path.getsize(original)
    fastqca_size = os.path.getsize(fastqca_out)
    gzip_size = os.path.getsize(gzip_out)
    return {
        "original": original_size,
        "fastqca": fastqca_size,
        "gzip": gzip_size,
        "fastqca_ratio": 1 - fastqca_size / original_size,
        "gzip_ratio": 1 - gzip_size / original_size,
    }
