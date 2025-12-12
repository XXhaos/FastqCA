import argparse
import os
from typing import Optional

from Bio import SeqIO

from Lossy_thread import (
    compress_multithread,
    decompress,
    delete_temp_files,
    get_output_path,
)


def count_reads(fastq_path: str) -> int:
    return sum(1 for _ in SeqIO.parse(fastq_path, "fastq"))


def write_manifest(manifest_path: str, read_count: int) -> None:
    with open(manifest_path, "w") as handle:
        handle.write(str(read_count))


def read_manifest(manifest_path: str) -> Optional[int]:
    if not os.path.exists(manifest_path):
        return None
    with open(manifest_path, "r") as handle:
        content = handle.read().strip()
        try:
            return int(content)
        except ValueError:
            return None


def main() -> None:
    lpaq8_path = f"{os.getcwd()}/lpaq8"
    parser = argparse.ArgumentParser(description="实验性多线程有损压缩脚本")
    parser.add_argument("--input_path", type=str, required=True, help="输入 FASTQ 路径或压缩文件路径")
    parser.add_argument("--output_path", type=str, required=True, help="输出路径")
    parser.add_argument("--mode", type=str, required=True, help="compress(c) 或 decompress(d)")
    parser.add_argument("--save", type=str, default="False", help="是否保留中间文件 (True/False)")
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help="工作进程数")
    parser.add_argument("--block_size", type=int, default=256 * 1024 * 1024, help="块大小 (字节)")
    parser.add_argument(
        "--manifest",
        type=str,
        default=None,
        help="存储/读取 read 数量的清单文件路径，默认与压缩文件同名并追加 .readcount",
    )
    args = parser.parse_args()

    save_flag = args.save.lower() == "true"

    if args.mode in ["compress", "c"]:
        archive_path = get_output_path(args.input_path, args.output_path)
        manifest_path = args.manifest or f"{archive_path}.readcount"
        read_count = count_reads(args.input_path)
        compress_multithread(args.input_path, args.output_path, lpaq8_path, save_flag, args.block_size, args.threads)
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
