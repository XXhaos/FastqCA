import argparse
import os

# 假设 Lossy 保持不变
from Lossy import main as lossy

# 【修改】这里改为导入 LossLess_thread
# 注意：确保 LossLess_thread.py 和 main.py 在同一个目录下
from LossLess_thread import compress_multithread, decompress

if __name__ == '__main__':
    lpaq8_path = f"{os.getcwd()}/lpaq8"

    # 创建参数解析器
    parser = argparse.ArgumentParser(description='fastq compress')

    # 添加参数
    parser.add_argument('--compressor', type=str, default='Lossy', help='Lossy or LossLess?')
    parser.add_argument('--input_path', type=str, required=True, help='input_path')
    parser.add_argument('--output_path', type=str, required=True, help='output_path')
    parser.add_argument('--mode', type=str, required=True, help='mode: c (compress) or d (decompress)')
    parser.add_argument('--save', type=str, default='False', help='save intermediate files (default: False)')

    # 多线程相关参数
    parser.add_argument('--threads', type=int, default=4, help='number of threads/workers')
    parser.add_argument('--block_size', type=int, default=256 * 1024 * 1024, help='block size in bytes')

    # 解析参数
    args = parser.parse_args()

    # 处理 save 参数
    save_flag = args.save.lower() == 'true'

    Lossy_commands = ['Lossy', 'lossy']
    # 扩充命令列表，包含 'Lossless_thread'
    LossLess_commands = ['LossLess', 'lossless', 'lossLess', 'Lossless', 'Lossless_thread']

    if args.compressor in Lossy_commands:
        lossy(args.mode, args.input_path, args.output_path, lpaq8_path, args.save, None)

    elif args.compressor in LossLess_commands:
        # 调用 LossLess_thread 中的函数
        if args.mode in ['c', 'compress']:
            compress_multithread(
                args.input_path,
                args.output_path,
                lpaq8_path,
                save_flag,
                args.block_size,
                args.threads
            )
        elif args.mode in ['d', 'decompress']:
            decompress(
                args.input_path,
                args.output_path,
                lpaq8_path,
                save_flag,
                None
            )
        else:
            print(f"错误：未知的模式 {args.mode}")
            exit(1)

    else:
        print("错误：没有指定正确的压缩器")
        exit(1)