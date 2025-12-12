"""实验性多线程有损压缩脚本。

基于 LossLess_script_new 风格的命令行参数，直接调用 `Lossy_thread` 中
的多线程压缩/解压实现，并在压缩阶段记录 read 数量、在解压阶段校验
read 数量，确保不会出现缺失或重复。
"""
import argparse
import os
import sys
from typing import Optional

from Bio import SeqIO

# 确保可以从当前脚本目录加载依赖模块
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from Lossy_thread import (  # noqa: E402
    compress_multithread,
    decompress,
    delete_temp_files,
    get_output_path,
)


def resolve_lpaq8_path() -> str:
    """返回 lpaq8 可执行文件的绝对路径。

    会依次检查脚本所在目录下的 `lpaq8` 与 `lpaq8.exe`，以便在不同平台
    上均可直接运行。
    """
    candidates = [os.path.join(SCRIPT_DIR, "lpaq8"), os.path.join(SCRIPT_DIR, "lpaq8.exe")]
    for path in candidates:
        if os.path.exists(path):
            return path
    # 默认返回第一个候选路径，调用方若不存在会得到清晰的报错
    return candidates[0]


def count_reads(fastq_path: str) -> int:
    """统计 FASTQ 文件中的 read 数量。"""
    return sum(1 for _ in SeqIO.parse(fastq_path, "fastq"))


def write_manifest(manifest_path: str, read_count: int) -> None:
    """将 read 数量写入清单文件。"""
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True) if os.path.dirname(manifest_path) else None
    with open(manifest_path, "w", encoding="utf-8") as handle:
        handle.write(str(read_count))


def read_manifest(manifest_path: str) -> Optional[int]:
    """读取清单文件中的 read 数量。"""
    if not os.path.exists(manifest_path):
        return None
    with open(manifest_path, "r", encoding="utf-8") as handle:
        content = handle.read().strip()
        try:
            return int(content)
        except ValueError:
            return None


def main() -> None:
    parser = argparse.ArgumentParser(description="fastq lossy compress optimized (multithread experimental)")
    parser.add_argument("--input_path", type=str, required=True, help="输入 FASTQ 或压缩文件路径")
    parser.add_argument("--output_path", type=str, required=True, help="输出路径：压缩文件或还原 FASTQ 目录/文件")
    parser.add_argument("--mode", type=str, required=True, help="compress(c) 或 decompress(d)")
    parser.add_argument("--save", type=str, default="False", help="是否保留中间文件 True/False")
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help="工作线程数")
    parser.add_argument("--block_size", type=int, default=256 * 1024 * 1024, help="块大小（字节）")
    parser.add_argument(
        "--manifest",
        type=str,
        default=None,
        help="read 计数清单文件路径，默认为压缩文件名追加 .readcount",
    )
    args = parser.parse_args()

    save_flag = args.save.lower() == "true"
    lpaq8_path = resolve_lpaq8_path()

    if args.mode in ["compress", "c"]:
        archive_path = get_output_path(args.input_path, args.output_path)
        manifest_path = args.manifest or f"{archive_path}.readcount"

        read_count = count_reads(args.input_path)
        compress_multithread(
            args.input_path,
            args.output_path,
            lpaq8_path,
            save_flag,
            args.block_size,
            args.threads,
        )
        write_manifest(manifest_path, read_count)
        if not save_flag:
            delete_temp_files(args.output_path)

    elif args.mode in ["decompress", "d"]:
        archive_path = args.input_path
        manifest_path = args.manifest or f"{archive_path}.readcount"

        decompress(archive_path, args.output_path, lpaq8_path, save_flag, None)
        restored_path = get_output_path(archive_path, args.output_path)
        restored_reads = count_reads(restored_path)
        expected_reads = read_manifest(manifest_path)
        if expected_reads is not None and restored_reads != expected_reads:
            raise RuntimeError(
                f"解压后的 read 数量 ({restored_reads}) 与清单记录的数量 ({expected_reads}) 不一致"
            )
        if not save_flag:
            delete_temp_files(args.output_path)
    else:
        raise SystemExit("Unknown mode, 请指定 compress/c 或 decompress/d")


if __name__ == "__main__":
    main()
