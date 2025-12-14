import ctypes
import io
import os
import subprocess
import threading
import time
from contextlib import contextmanager

_LIB_NAME = "liblpaq8.so"
_MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
_LIB_PATH = os.path.join(_MODULE_DIR, _LIB_NAME)
_lpaq8_lib = None


def _build_shared_library():
    if os.path.exists(_LIB_PATH):
        return
    cmd = [
        "g++",
        "-std=c++11",
        "-O3",
        "-DNDEBUG",
        "-fPIC",
        "-shared",
        "-o",
        _LIB_NAME,
        "lpaq8.cpp",
    ]
    subprocess.run(cmd, check=True, cwd=_MODULE_DIR)


def _load_library():
    global _lpaq8_lib
    if _lpaq8_lib is None:
        _build_shared_library()
        _lpaq8_lib = ctypes.CDLL(_LIB_PATH)
        _lpaq8_lib.lpaq8_main.restype = ctypes.c_int
    return _lpaq8_lib


@contextmanager
def _memfd(name, data=None):
    if hasattr(os, "memfd_create"):
        fd = os.memfd_create(name)
    else:
        fd = os.open(os.path.join(_MODULE_DIR, f".{name}"), os.O_RDWR | os.O_CREAT | os.O_TRUNC)
    try:
        if data:
            os.write(fd, data)
            os.lseek(fd, 0, os.SEEK_SET)
        yield fd, f"/proc/self/fd/{fd}"
    finally:
        try:
            os.close(fd)
        except OSError:
            pass
        # cleanup fallback file if created on disk
        fallback_path = os.path.join(_MODULE_DIR, f".{name}")
        if os.path.exists(fallback_path):
            try:
                os.remove(fallback_path)
            except OSError:
                pass


def _run_lpaq8(argv):
    lib = _load_library()
    argc = len(argv)
    c_argv = (ctypes.c_char_p * argc)(*argv)
    ret = lib.lpaq8_main(argc, c_argv)
    if ret != 0:
        raise RuntimeError(f"lpaq8 native 调用失败，返回码 {ret}")


def compress_bytes(data: bytes, compression_level: str = '9') -> bytes:
    argv = [
        b"lpaq8",
        compression_level.encode(),
    ]
    with _memfd("lpaq8_in", data) as (in_fd, in_path), _memfd("lpaq8_out") as (out_fd, out_path):
        argv.extend([in_path.encode(), out_path.encode()])
        _run_lpaq8(argv)
        os.lseek(out_fd, 0, os.SEEK_SET)
        return os.read(out_fd, os.fstat(out_fd).st_size)


def decompress_bytes(data: bytes) -> bytes:
    argv = [
        b"lpaq8",
        b"d",
    ]
    with _memfd("lpaq8_in", data) as (in_fd, in_path), _memfd("lpaq8_out") as (out_fd, out_path):
        argv.extend([in_path.encode(), out_path.encode()])
        _run_lpaq8(argv)
        os.lseek(out_fd, 0, os.SEEK_SET)
        return os.read(out_fd, os.fstat(out_fd).st_size)


def compress_file(input_file, output_file, lpaq8_path=None, compression_level='9'):
    with open(input_file, "rb") as f:
        data = f.read()
    compressed = compress_bytes(data, compression_level)
    with open(output_file, "wb") as f:
        f.write(compressed)
    return None


def decompress_file(input_file, output_file, lpaq8_path=None):
    with open(input_file, "rb") as f:
        data = f.read()
    decompressed = decompress_bytes(data)
    with open(output_file, "wb") as f:
        f.write(decompressed)
    return None


def compress_all_files_in_directory(input_directory, output_directory, lpaq8_path, compression_level='9'):
    """
    压缩目录中的所有文件。

    参数:
    - input_directory: 输入目录路径。
    - output_directory: 输出目录路径。
    - lpaq8_path: lpaq8压缩器的完整路径。
    - compression_level: 压缩级别（默认为9，范围0-9）。
    """
    # 记录开始时间
    start_time = time.time()

    # 确保输出目录存在
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for root, dirs, files in os.walk(input_directory):
        for file in files:
            input_file_path = os.path.join(root, file)
            compressed_filename = f"{os.path.splitext(os.path.basename(input_file_path))[0]}.lpaq8"
            output_path = os.path.join(output_directory, compressed_filename)

            compress_file(input_file_path, output_path, lpaq8_path, compression_level)

    # 记录结束时间
    end_time = time.time()

    print(f"所有文件已压缩完成。总共耗时: {(end_time - start_time) / 60} 分钟。")


def decompress_all_files_in_directory(input_directory, output_directory, lpaq8_path):
    """
    解压目录中的所有文件。

    参数:
    - input_directory: 输入目录路径。
    - output_directory: 输出目录路径。
    - lpaq8_path: lpaq8压缩器的完整路径。
    - compression_level: 压缩级别（默认为9，范围0-9）。
    """
    # 记录开始时间
    start_time = time.time()

    # 确保输出目录存在
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    output_path_mapping = {
        "id_regex.lpaq8": "txt",
        "id_tokens.lpaq8": "txt",
        "base_g_prime.lpaq8": "tiff",
        "quality.lpaq8": "tiff"
    }

    for root, dirs, files in os.walk(input_directory):
        for file in files:
            error = True
            input_file_path = os.path.join(root, file)
            base_filename = os.path.splitext(os.path.basename(input_file_path))[0]

            for suffix in output_path_mapping.keys():
                if file.endswith(suffix):
                    output_filename = f"{base_filename}.{output_path_mapping[suffix]}"
                    output_path = os.path.join(output_directory, output_filename)
                    decompress_file(input_file_path, output_path, lpaq8_path)
                    error = False

            if error:
                print(f"未知文件类型: {file}")

    # 记录结束时间
    end_time = time.time()

    print(f"所有文件已解压完成。总共耗时: {(end_time - start_time) / 60} 分钟。")


def get_file_size(file_path):
    file_size = os.path.getsize(file_path)

    if file_size < 1024:
        return f"{file_size} bytes"
    elif file_size < 1024 * 1024:
        return f"{file_size / 1024:.2f} KB"
    elif file_size < 1024 * 1024 * 1024:
        return f"{file_size / (1024 * 1024):.2f} MB"
    else:
        return f"{file_size / (1024 * 1024 * 1024):.2f} GB"


def get_directory_size(directory_path):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)

    if total_size < 1024:
        return f"{total_size} bytes"
    elif total_size < 1024 * 1024:
        return f"{total_size / 1024:.2f} KB"
    elif total_size < 1024 * 1024 * 1024:
        return f"{total_size / (1024 * 1024):.2f} MB"
    else:
        return f"{total_size / (1024 * 1024 * 1024):.2f} GB"

def monitor_output_file(output_file):
    while True:
        file_size = os.path.getsize(output_file)
        print(f"Output file size: {file_size} bytes")
        time.sleep(1)  # 每隔一秒检查一次文件大小


if __name__ == '__main__':
    # 示例用法
    input_directory1 = r"D:\pythonProject\fastqtobmp\input\change_to_gray" # 定义需要压缩的文件路径
    destination_directory1 = r'D:\pythonProject\fastqtobmp\input\change_to_gray_lpaq8'  # 定义输出目录
    lpaq8_exe_path = f"{os.getcwd()}\lpaq8.exe"  # 确保这是正确的lpaq8路径


    input_destination = r"D:\pythonProject\fastqtobmp\input"
    output_destination = r"D:\pythonProject\fastqtobmp\output\1"

    output_file = os.path.join("output", "SRR554369")
    monitor_thread = threading.Thread(target=monitor_output_file, args=(output_file, ))
    monitor_thread.start()

    compress_file(os.path.join(os.getcwd(), "input", "SRR554369.fastq"), os.path.join(os.getcwd(), "output", "SRR554369"), lpaq8_exe_path)

    monitor_thread.join()

    # input_directory2 = r"D:\pythonProject\fastqtobmp\input\compressed" # 定义需要压缩的文件路径
    # destination_directory2 = r'D:\pythonProject\fastqtobmp\input\compressed_lpaq8'  # 定义输出目录

    # 压缩目录中的所有文件
    # compress_all_files_in_directory(input_directory1, destination_directory1, lpaq8_exe_path)
    # compress_all_files_in_directory(input_directory2, destination_directory2, lpaq8_exe_path)

    # 计算输出目录的大小，并转换为MB
    # size1 = get_directory_size(destination_directory1)
    # size2 = get_directory_size(destination_directory2)

    # 输出两个目录大小的比较结果
    # difference = size1 - size2
    # print(f"{destination_directory1} 比 {destination_directory2} 大了 {difference:.2f} MB。")
