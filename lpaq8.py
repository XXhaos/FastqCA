"""Python wrappers around the external LPAQ8 executable used by FastqCA."""

import os
import subprocess
import threading
import time


def compress_lpaq8(lpaq8_path, input_path, output_path, compression_level='9'):
    """Start LPAQ8 compression for one input stream.

    The default compression level is 9, matching the high-density setting used
    by FastqCA in the benchmark workflow.
    """
    command = [lpaq8_path, compression_level, input_path, output_path]
    try:
        process = subprocess.Popen(command)
        return process
        # print(f"文件 {input_path} 压缩成功，保存为 {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"压缩过程中出错: {e}")
    except Exception as e:
        print(f"发生未知错误: {str(e)}")

def compress_lpaq8_test(lpaq8_path, input_stream, output_path, compression_level='9'):
    """Run LPAQ8 compression synchronously for a single input stream."""
    command = [lpaq8_path, compression_level, input_stream, output_path]
    try:
        subprocess.run(command, check=True)
        print(f"文件 {input_stream} 压缩成功，保存为 {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"压缩过程中出错: {e}")
    except Exception as e:
        print(f"发生未知错误: {str(e)}")


def decompress_lpaq8(lpaq8_path, input_path, output_path):
    """Start LPAQ8 decompression for one compressed stream."""
    command = [lpaq8_path, 'd', input_path, output_path]
    try:
        process = subprocess.Popen(command)
        # print(f"文件 {input_path} 解压成功，保存为 {output_path}")
        return process
    except subprocess.CalledProcessError as e:
        print(f"解压过程中出错: {e}")
    except Exception as e:
        print(f"发生未知错误: {str(e)}")


def compress_file(input_file, output_file, lpaq8_path, compression_level='9'):
    """Compress one file with the configured LPAQ8 executable."""
    # Delegate the actual entropy coding to the external LPAQ8 executable.
    return compress_lpaq8(lpaq8_path, input_file, output_file, compression_level)


def decompress_file(input_file, output_file, lpaq8_path):
    """Decompress one LPAQ8-compressed stream."""
    return decompress_lpaq8(lpaq8_path, input_file, output_file)


def compress_all_files_in_directory(input_directory, output_directory, lpaq8_path, compression_level='9'):
    """Compress all files in a directory with LPAQ8."""
    # Record elapsed time for console reporting only.
    start_time = time.time()

    # Ensure the destination directory exists before launching LPAQ8 jobs.
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for root, dirs, files in os.walk(input_directory):
        for file in files:
            input_file_path = os.path.join(root, file)
            compressed_filename = f"{os.path.splitext(os.path.basename(input_file_path))[0]}.lpaq8"
            output_path = os.path.join(output_directory, compressed_filename)

            compress_file(input_file_path, output_path, lpaq8_path, compression_level)

    # Report total wall-clock time for this batch helper.
    end_time = time.time()

    print(f"所有文件已压缩完成。总共耗时: {(end_time - start_time) / 60} 分钟。")


def decompress_all_files_in_directory(input_directory, output_directory, lpaq8_path):
    """Decompress all recognized LPAQ8 stream files in a directory."""
    # Record elapsed time for console reporting only.
    start_time = time.time()

    # Ensure the destination directory exists before writing restored streams.
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

    # Report total wall-clock time for this batch helper.
    end_time = time.time()

    print(f"所有文件已解压完成。总共耗时: {(end_time - start_time) / 60} 分钟。")


def get_file_size(file_path):
    """Return a human-readable 1024-based file-size string."""
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
    """Return a recursive human-readable directory-size string."""
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
    """Print the size of an output file periodically for manual monitoring."""
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
