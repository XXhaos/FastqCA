import argparse
import io
import os
import shutil
import sys
import time
from tqdm import tqdm
from PIL import Image, UnidentifiedImageError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import defaultdict
from itertools import product
from PIL import Image
from Bio import SeqIO
from tqdm import tqdm
import re
from lpaq8 import compress_file, decompress_file

# 定义碱基到灰度值的映射
base_to_gray = {'A': 32, 'T': 64, 'G': 192, 'C': 224, 'N': 0}
# 定义碱基从灰度值的映射
gray_to_base = {32: 'A', 64: 'T', 192: 'G', 224: 'C', 0: 'N'}


def find_delimiters(identifier):
    delimiters = re.findall(r'[.:_\s=/-]', identifier)
    return delimiters


def split_identifier(identifier, delimiters):
    tokens = re.split(r'[.:_\s=/-]', identifier)
    return tokens


def generate_regex(delimiters):
    regex = ""
    count = 1
    for delimiter in delimiters:
        regex += f"{{{{T{count}}}}}"  # 使用{{T1}}格式
        regex += delimiter
        count += 1
    regex += f"{{{{T{count}}}}}"
    return regex


# 初始化规则字典，记录规则的使用次数
def init_rules_dict():
    # 定义规则表
    rule_table = []

    # 添加所有可能的规则到规则库
    values = [0, 32, 64, 192, 224]
    combinations = list(product(values, repeat=4))

    for combination in combinations:
        rule_table.append(combination)

    # 初始化规则字典，记录规则的使用次数
    rules_dict = defaultdict(int)

    return rules_dict


# 获取分块需要读取的reads数
def get_reads_num_per_block(fastq_path, block_size):
    # 使用 SeqIO.parse() 预读取第一个 read 来确定 read 长度
    with open(fastq_path, 'r') as file:
        first_record = next(SeqIO.parse(file, "fastq"))
        read_length = len(first_record.seq)

    # 计算每个read所需的字节数（DNA序列和质量分数各占一半）
    bytes_per_read = read_length * 2

    # 根据block_size和bytes_per_read计算每个图像块应包含的reads数目
    reads_per_block = block_size // bytes_per_read
    total_reads = os.path.getsize(fastq_path) // bytes_per_read

    return reads_per_block, total_reads

def generate_g_prime(G, rules_dict, p_bar, gr_bar):
    # 初始化G_prime为与G相同大小的零矩阵
    G_prime = np.zeros_like(G)

    # 按照规则表进行预测并生成G_prime矩阵
    for i in range(0, G.shape[0]):
        for j in range(0, G.shape[1]):
            center = G[i, j]
            if i == 0:
                up = 0
            else:
                up = G[i - 1, j]
            if j == 0:
                left = 0
            else:
                left = G[i, j - 1]
            if i != 0 and j != 0:
                left_up = G[i - 1, j - 1]
            else:
                left_up = 0
            matched_rule = (up, left_up, left, center)
            matched_rule1 = (up, left_up, left, 32)
            matched_rule2 = (up, left_up, left, 224)
            matched_rule3 = (up, left_up, left, 192)
            matched_rule4 = (up, left_up, left, 64)
            matched_rule5 = (up, left_up, left, 0)

            matched_rules = [matched_rule1, matched_rule2, matched_rule3, matched_rule4, matched_rule5]
            max_freq = -1
            top_rule = None

            for rule in matched_rules:
                freq = rules_dict[rule]
                if freq > max_freq:
                    max_freq = freq
                    top_rule = rule

            # 如果规则中心与当前元素相同，则将G_prime[i, j]赋值为1，否则赋值为当前元素
            G_prime[i, j] = 1 if top_rule[3] == center else center
            # 更新规则字典中该规则的频率
            rules_dict[matched_rule] += 1
        p_bar.update(0.5)

        if gr_bar:
            gr_bar.update(5)

    return G_prime


# 对records做压缩
def process_records(records, rules_dict, p_bar, gr_bar):
    # 对原始数据进行处理
    id_block = []
    base_image_block = []
    quality_block = []

    start_time = time.time()

    # 生成g_prime，quality_block， id_block
    for i, record in enumerate(records):
        id_str = record.description
        delimiters = find_delimiters(id_str)
        tokens = split_identifier(id_str, delimiters)
        regex = generate_regex(delimiters)
        id_block.append((tokens, regex))

        base_gray_values = [base_to_gray.get(base, 0) for base in record.seq]
        base_image_block.append(base_gray_values)

        quality_gray_values = [q * 2 for q in record.letter_annotations["phred_quality"]]
        quality_block.append(quality_gray_values)

        p_bar.update(0.1)
        if gr_bar:
            gr_bar.update(1)
        # # 计算已经过去的平均时间
        # elapsed_time = time.time() - start_time
        # avg_time_per_iteration = elapsed_time / (i + 1)
        # # 计算剩余时间
        # remaining_iterations = len(records) - (i + 1)
        # estimated_remaining_time = remaining_iterations * avg_time_per_iteration
        # # 更新剩余时间信息
        # p_bar.set_postfix(ETA=f"预计前端压缩剩余时间：{estimated_remaining_time:.1f}s")

    g_prime = generate_g_prime(np.array(base_image_block, dtype=np.uint8), rules_dict, p_bar, gr_bar)
    quality_block = np.array(quality_block, dtype=np.uint8)

    p_bar.update(p_bar.total * 0.5 - p_bar.n)
    if gr_bar:
        gr_bar.update(np.ceil(10 * (p_bar.total * 0.5 - p_bar.n)))
    g_prime_img = Image.fromarray(g_prime.astype(np.uint8))
    quality_img = Image.fromarray(quality_block.astype(np.uint8))

    return g_prime_img, quality_img, id_block


def save_G_block(g_prime_img, output_path, block_count):
    # 构建输出文件路径
    g_prime_file_path = os.path.join(os.path.dirname(output_path), "front_compressed", f'chunk_{block_count}_base.tiff')

    # 将G_prime矩阵转换为图像并保存为TIFF文件
    g_prime_img.save(g_prime_file_path)


def save_quality_block(quality_img, output_path, block_count):
    # 构建输出文件路径
    quality_file_path = os.path.join(os.path.dirname(output_path), "front_compressed",
                                     f'chunk_{block_count}_quality.tiff')

    # 将G矩阵转换为图像并保存为TIFF文件
    quality_img.save(quality_file_path)


def save_id_block(id_block, output_path, block_count):
    tokens_file_path = os.path.join(os.path.dirname(output_path), "front_compressed",
                                    f"chunk_{block_count}_id_tokens.txt")
    regex_file_path = os.path.join(os.path.dirname(output_path), "front_compressed",
                                   f"chunk_{block_count}_id_regex.txt")

    with open(tokens_file_path, 'w') as tokens_file, open(regex_file_path, 'w') as regex_file:
        for tokens, regex in id_block:
            tokens_file.write(' '.join(tokens) + '\n')
            regex_file.write(regex + '\n')


def front_compress(records, rules_dict, block_count, output_path, save, p_bar, gr_bar):
    tqdm.write(f"info：开始前端压缩，当前压缩的是第{block_count}块")
    g_block, quality_block, id_block = process_records(records, rules_dict, p_bar, gr_bar)
    tqdm.write(f"info：前端压缩完成，当前压缩的是第{block_count}块")

    if save:
        tqdm.write("info：正在保存前端压缩文件")

        # 确保输出路径存在
        output_dir = os.path.join(os.path.dirname(output_path), "front_compressed")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        save_G_block(g_block, output_path, block_count)
        save_quality_block(quality_block, output_path, block_count)
        save_id_block(id_block, output_path, block_count)

        tqdm.write("info：前端压缩文件保存完成")

    return g_block, quality_block, id_block

    # tqdm.write(f"处理的reads数量: {read_count}，每个图像块的reads数量: {reads_per_block}，一共有: {total_blocks}个图像块。")
    # tqdm.write(f"前端压缩完成，压缩前文件大小：{get_file_size(fastq_path)}")


# def first_decompress(fastq_path, first_compressed_path):
#     id_token_files = sorted(
#         glob.glob(f"{first_compressed_path}/chunk_*_id_tokens.txt"),
#         key=lambda x: int(re.search(r'chunk_(\d+)_id_tokens.txt', x).group(1))
#     )
#     id_regex_files = sorted(
#         glob.glob(f"{first_compressed_path}/chunk_*_id_regex.txt"),
#         key=lambda x: int(re.search(r'chunk_(\d+)_id_regex.txt', x).group(1))
#     )
#     base_files = sorted(
#         glob.glob(f"{first_compressed_path}/chunk_*_base_g_prime.tiff"),
#         key=lambda x: int(re.search(r'chunk_(\d+)_base_g_prime.tiff', x).group(1))
#     )
#     quality_files = sorted(
#         glob.glob(f"{first_compressed_path}/chunk_*_quality.tiff"),
#         key=lambda x: int(re.search(r'chunk_(\d+)_quality.tiff', x).group(1))
#     )
#     output_fastq_path = r"D:\pythonProject\fastqtobmp\input\reconstructed.fastq"
#     reconstruct_fastq(id_token_files, id_regex_files, base_files, quality_files, output_fastq_path)
#     tqdm.write("FASTQ文件已还原。")

def monitor(process, temp_input_path, temp_output_path, p_bar, gr_bar, progress):
    while True:

        # 可能此时temp_input_path和temp_output_path还没有生成
        if os.path.exists(temp_input_path) and os.path.exists(temp_output_path):
            input_file_size = os.path.getsize(temp_input_path)
            output_file_size = os.path.getsize(temp_output_path)

            if p_bar:
                p_bar.update(output_file_size / input_file_size * p_bar.total * progress / 100)

            if gr_bar:
                gr_bar.update(np.ceil(output_file_size / input_file_size * p_bar.total * progress / 100 * 10))

            if process.poll() is not None:
                break


def compress_with_monitor(temp_input_path, temp_output_path, lpaq8_path, p_bar, gr_bar, progress):
    process = compress_file(temp_input_path, temp_output_path, lpaq8_path)

    monitor(process, temp_input_path, temp_output_path, p_bar, gr_bar, progress)


def decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path, p_bar, gr_bar, progress):
    process = decompress_file(temp_input_path, temp_output_path, lpaq8_path)

    monitor(process, temp_input_path, temp_output_path, p_bar, gr_bar, progress)


def back_compress(g_block, quality_block, id_block, lpaq8_path, output_path, save, block_count, p_bar, gr_bar):
    tqdm.write(f"info：开始后端压缩，当前压缩的是第{block_count}块")

    output_dir = os.path.join(os.path.dirname(output_path), "back_compressed")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    temp_input_path = os.path.join(os.path.dirname(output_path), "temp_input")
    temp_output_path = os.path.join(os.path.dirname(output_path), "temp_output")

    with open(output_path, "a+b") as output_file:

        # 对id_regex进行处理
        with open(temp_input_path, "w") as temp_input_file:
            for tokens, regex in id_block:
                temp_input_file.write(regex + '\n')

        compress_with_monitor(temp_input_path, temp_output_path, lpaq8_path, p_bar, gr_bar, 5)

        with open(temp_output_path, "rb") as temp_output_file:
            output_file.write(b"%id_regex%")
            output_file.write(temp_output_file.read())

            if save:
                save_path = os.path.join(output_dir, f"chunk_{block_count}_id_regex.lpaq8")
                temp_output_file.seek(0)
                with open(save_path, "w+b") as save_file:
                    save_file.write(temp_output_file.read())

        # 对id_tokens进行处理
        with open(temp_input_path, "w") as temp_input_file:
            for tokens, regex in id_block:
                temp_input_file.write(' '.join(tokens) + '\n')

        compress_with_monitor(temp_input_path, temp_output_path, lpaq8_path, p_bar, gr_bar, 5)

        with open(temp_output_path, "rb") as temp_output_file:
            output_file.write(b"%id_tokens%")
            output_file.write(temp_output_file.read())

            if save:
                save_path = os.path.join(output_dir, f"chunk_{block_count}_id_tokens.lpaq8")
                temp_output_file.seek(0)
                with open(save_path, "w+b") as save_file:
                    save_file.write(temp_output_file.read())

        # 对g_block进行处理
        g_block.save(temp_input_path, format="tiff")

        compress_with_monitor(temp_input_path, temp_output_path, lpaq8_path, p_bar, gr_bar, 15)

        with open(temp_output_path, "rb") as temp_output_file:
            output_file.write(b"%base_g_prime%")
            output_file.write(temp_output_file.read())

            if save:
                save_path = os.path.join(output_dir, f"chunk_{block_count}_base_g_prime.lpaq8")
                temp_output_file.seek(0)
                with open(save_path, "w+b") as save_file:
                    save_file.write(temp_output_file.read())

        # 对quality进行处理
        quality_block.save(temp_input_path, format="tiff")

        compress_with_monitor(temp_input_path, temp_output_path, lpaq8_path, p_bar, gr_bar, 15)

        with open(temp_output_path, "rb") as temp_output_file:
            output_file.write(b"%quality%")
            output_file.write(temp_output_file.read())

            if save:
                save_path = os.path.join(output_dir, f"chunk_{block_count}_quality.lpaq8")
                temp_output_file.seek(0)
                with open(save_path, "w+b") as save_file:
                    save_file.write(temp_output_file.read())

    p_bar.update(p_bar.total - p_bar.n)

    if gr_bar:
        gr_bar.update(np.ceil(10 * (p_bar.total - p_bar.n)))
    tqdm.write(f"info：后端压缩完成，当前压缩的是第{block_count}块")


def process_block(records, rules_dict, block_count, output_path, lpaq8_path, save, p_bar, gr_bar):
    # 第一步处理：调用函数处理FASTQ文件，并打印处理的reads数量、每块reads数量和图像块数量
    g_block, quality_block, id_block = front_compress(records, rules_dict, block_count, output_path, save, p_bar,
                                                      gr_bar)

    # 第二步处理：将处理后的文件进行后端压缩
    back_compress(g_block, quality_block, id_block, lpaq8_path, output_path, save, block_count, p_bar, gr_bar)


def process_last_block(output_path):
    with open(output_path, "a+b") as output_file:
        output_file.write(b"%eof")


# 对output_path进行处理，允许使用路径或文件作为output_path
def get_output_path(input_path, output_path):
    if input_path is None or not os.path.isfile(input_path):
        tqdm.write("错误：输入路径不存在")
        exit(1)

    if os.path.isfile(output_path):
        return output_path

    basename = os.path.splitext(os.path.basename(input_path))[0]

    if os.path.isdir(output_path):
        return os.path.join(output_path, basename)


def compress(fastq_path, output_path, lpaq8_path, save, gr_progress, block_size = 256*1024*1024):
    output_path = get_output_path(fastq_path, output_path)
    gr_bar = None

    # 初始化规则字典
    rules_dict = init_rules_dict()

    # 获取每块应读取的reads条数
    reads_per_block, total_reads = get_reads_num_per_block(fastq_path, block_size)
    total_block_count = np.ceil(os.path.getsize(fastq_path) / block_size)

    # 将reads_per_block的reads组装成一个records
    records = []
    read_count_per_block = 0
    # read_count = 0
    block_count = 1

    # tqdm.write(f"info：处理的reads数量: {read_count}，每个图像块的reads数量: {reads_per_block}, 一共有: {total_block_count}个图像块。")
    tqdm.write(f"info：每个块的reads数量: {reads_per_block}, 一共有: {total_block_count}个图像块。")

    tqdm.write(f"info：开始读取fastq文件，文件名：{os.path.splitext(os.path.basename(fastq_path))[0]}")
    # 使用 SeqIO.parse() 函数逐条读取 FASTQ 文件
    with open(fastq_path, 'r') as file:
        for record in SeqIO.parse(file, "fastq"):
            records.append(record)

            read_count_per_block += 1
            # read_count += 1

            # 当达到指定的 reads_per_block 时，处理文件块
            if read_count_per_block == reads_per_block:
                with tqdm(total=read_count_per_block, file=sys.stdout, colour='red',
                          desc=f"block:{block_count} / {total_block_count}", dynamic_ncols=True,
                          bar_format='{l_bar}{bar}| {n:.3f}/{total_fmt} [{elapsed}<{remaining}, ' '{rate_fmt}{postfix}]') as p_bar:

                    if gr_progress:
                        gr_bar = gr_progress.tqdm(range(10 * total_reads))
                    process_block(records, rules_dict, block_count, output_path, lpaq8_path, save, p_bar, gr_bar)

                    block_count += 1
                    read_count_per_block = 0
                    records = []
                    # 重置规则字典
                    rules_dict = init_rules_dict()

    # 处理最后一个块（如果有剩余）
    if records:
        with tqdm(total=read_count_per_block, file=sys.stdout, colour='red',
                  desc=f"block:{block_count} / {total_block_count}", dynamic_ncols=True,
                  bar_format='{l_bar}{bar}| {n:.3f}/{total_fmt} [{elapsed}<{remaining}, ' '{rate_fmt}{postfix}]') as p_bar:
            process_block(records, rules_dict, block_count, output_path, lpaq8_path, save, p_bar, gr_bar)
            process_last_block(output_path)
            block_count += 1


def process_compressed_block(output_path, lpaq8_path, id_regex_data, id_tokens_data, g_prime_data, quality_data, save, block_count):
    # 生成保存文件的文件夹
    back_compress_dir = os.path.join(os.path.dirname(output_path), "back_compressed")
    front_compress_dir = os.path.join(os.path.dirname(output_path), "front_compressed")
    if save:
        # 确保输出文件夹存在
        if not os.path.exists(back_compress_dir):
            os.makedirs(back_compress_dir, exist_ok=True)
        if not os.path.exists(front_compress_dir):
            os.makedirs(front_compress_dir, exist_ok=True)

    # 定义不同的临时文件路径，每个数据块使用唯一的路径
    temp_input_path_id_regex = os.path.join(os.path.dirname(output_path), f"temp_input_id_regex_{block_count}")
    temp_output_path_id_regex = os.path.join(os.path.dirname(output_path), f"temp_output_id_regex_{block_count}")

    temp_input_path_id_tokens = os.path.join(os.path.dirname(output_path), f"temp_input_id_tokens_{block_count}")
    temp_output_path_id_tokens = os.path.join(os.path.dirname(output_path), f"temp_output_id_tokens_{block_count}")

    temp_input_path_quality = os.path.join(os.path.dirname(output_path), f"temp_input_quality_{block_count}")
    temp_output_path_quality = os.path.join(os.path.dirname(output_path), f"temp_output_quality_{block_count}")

    temp_input_path_g_prime = os.path.join(os.path.dirname(output_path), f"temp_input_g_prime_{block_count}")
    temp_output_path_g_prime = os.path.join(os.path.dirname(output_path), f"temp_output_g_prime_{block_count}")

    id_regex = None
    id_tokens = None
    g_prime = None
    quality = None
    id_block = None

    try:
        # 处理 id_regex 数据块
        with open(temp_input_path_id_regex, "wb") as temp_input_file:
            with open(temp_output_path_id_regex, "w+") as temp_output_file:
                temp_input_file.write(id_regex_data)
                temp_input_file.flush()
                # 调用已有的解压缩函数
                decompress_with_monitor(temp_input_path_id_regex, temp_output_path_id_regex, lpaq8_path, None, None, None)
                id_regex = [line.strip() for line in temp_output_file.readlines()]

                if save:
                    shutil.copy(temp_input_path_id_regex, os.path.join(back_compress_dir, f"chunk_{block_count}_id_regex.lpaq8"))
                    shutil.copy(temp_output_path_id_regex, os.path.join(front_compress_dir, f"chunk_{block_count}_id_regex.txt"))

                tqdm.write("id_regex处理完毕")

        # 处理 id_tokens 数据块
        with open(temp_input_path_id_tokens, "wb") as temp_input_file:
            with open(temp_output_path_id_tokens, "w+") as temp_output_file:
                temp_input_file.write(id_tokens_data)
                temp_input_file.flush()
                # 调用已有的解压缩函数
                decompress_with_monitor(temp_input_path_id_tokens, temp_output_path_id_tokens, lpaq8_path, None, None, None)
                id_tokens = [line.strip() for line in temp_output_file.readlines()]

                if save:
                    shutil.copy(temp_input_path_id_tokens, os.path.join(back_compress_dir, f"chunk_{block_count}_id_tokens.lpaq8"))
                    shutil.copy(temp_output_path_id_tokens, os.path.join(front_compress_dir, f"chunk_{block_count}_id_tokens.txt"))

                tqdm.write("id_tokens处理完毕")

        id_block = zip(id_tokens, id_regex)

        # 处理 quality 数据块
        with open(temp_input_path_quality, "wb") as temp_input_file:
            temp_input_file.write(quality_data)
            temp_input_file.flush()
            # 调用已有的解压缩函数
            decompress_with_monitor(temp_input_path_quality, temp_output_path_quality, lpaq8_path, None, None, None)
            # 确保文件写入完成
            try:
                with Image.open(temp_output_path_quality) as img:
                    quality = img.copy()
            except UnidentifiedImageError:
                tqdm.write(f"无法识别的图像文件: {temp_output_path_quality}")
                raise

            if save:
                shutil.copy(temp_input_path_quality, os.path.join(back_compress_dir, f'chunk_{block_count}_quality.lpaq8'))
                shutil.copy(temp_output_path_quality, os.path.join(front_compress_dir, f'chunk_{block_count}_quality.tiff'))

            tqdm.write("quality处理完毕")

        # 处理 g_prime 数据块
        with open(temp_input_path_g_prime, "wb") as temp_input_file:
            temp_input_file.write(g_prime_data)
            temp_input_file.flush()
            # 调用已有的解压缩函数
            decompress_with_monitor(temp_input_path_g_prime, temp_output_path_g_prime, lpaq8_path, None, None, None)
            # 确保文件写入完成
            try:
                with Image.open(temp_output_path_g_prime) as img:
                    g_prime = img.copy()
            except UnidentifiedImageError:
                tqdm.write(f"无法识别的图像文件: {temp_output_path_g_prime}")
                raise

            if save:
                shutil.copy(temp_input_path_g_prime, os.path.join(back_compress_dir, f'chunk_{block_count}_base_g_prime.lpaq8'))
                shutil.copy(temp_output_path_g_prime, os.path.join(front_compress_dir, f'chunk_{block_count}_base_g_prime.tiff'))

            tqdm.write("g_prime处理完毕")

    finally:
        pass
    return id_block, g_prime, quality


def decompress(compressed_path, output_path, lpaq8_path, save, gr_progress):
    output_path = get_output_path(compressed_path, output_path)

    id_regex_seperator = b"%id_regex%"
    id_tokens_seperator = b"%id_tokens%"
    base_g_prime_seperator = b"%base_g_prime%"
    quality_seperator = b"%quality%"
    eof_seperator = b"%eof"

    start = 0
    end = 0

    block_count = 1

    id_regex_dat = None
    id_tokens_data = None
    g_prime_data = None
    quality_data = None

    with open(compressed_path, "rb") as input_file:
        content = input_file.read()

        # 取id_regex
        start = content.find(id_regex_seperator, end)
        eof = content.find(eof_seperator, start)

        if start == -1:
            tqdm.write(f"错误：文件开始符缺失")
            exit(1)

        if eof == -1:
            tqdm.write(f"错误：文件结束符缺失")
            exit(1)

        while True:

            # 获取id_regex
            start = content.find(id_regex_seperator, end)
            end = content.find(id_tokens_seperator, start)
            if start == -1:
                break
            if end < start:
                tqdm.write("错误：id_regex_seperator或id_tokens_seperator文件分隔符缺失")
                exit(1)
            id_regex_data = content[start + len(id_regex_seperator): end]

            # 获取id_tokens
            start = content.find(id_tokens_seperator, end)
            end = content.find(base_g_prime_seperator, start)
            if start == -1 or end < start:
                tqdm.write("错误：id_tokens_seperator或base_g_prime_seperator文件分隔符缺失")
                exit(1)
            id_tokens_data = content[start + len(id_tokens_seperator): end]

            # 获取g_prime
            start = content.find(base_g_prime_seperator, end)
            end = content.find(quality_seperator, start)
            if start == -1 or end < start:
                tqdm.write("错误：base_g_prime_seperator或quality_seperator文件分隔符缺失")
                exit(1)
            g_prime_data = content[start + len(base_g_prime_seperator): end]

            # 获取quality
            start = content.find(quality_seperator, end)
            end = content.find(id_regex_seperator, start)

            # 如果end = -1，说明已经到达文件尾
            if end == -1:
                end = eof

            if start == -1 or end < start:
                tqdm.write("错误：quality_seperator或id_regex_seperator文件分隔符缺失")
                exit(1)
            quality_data = content[start + len(quality_seperator): end]

            if id_tokens_data is None or g_prime_data is None or quality_data is None or id_regex_data is None:
                tqdm.write("错误：无法读取到id_regex_data或id_tokens_data或g_prime_data, quality_data")
                exit(1)

            id_block, g_prime, quality = process_compressed_block(output_path, lpaq8_path, id_regex_data,
                                                                  id_tokens_data,
                                                                  g_prime_data, quality_data, save, block_count)

            if id_block is None or g_prime is None or quality is None:
                tqdm.write("错误：无法重建id_block, g_prime, quality")
                exit(1)

            block_count += 1

            # # 判断g_prime_data, quality_data是否相同
            # if g_prime_data == quality_data:
            #     print("错误：g_prime_data, quality_data相同")
            #     exit(1)
            # else:
            #     print("g_prime_data, quality_data不同")
            #
            # # 判断g_prime, quality是否相同
            # if g_prime == quality:
            #     print("错误：g_prime, quality相同")
            #     exit(1)
            # else:
            #     print("g_prime, quality不同")

            reconstruct_fastq(output_path, id_block, g_prime, quality)


def load_id_block(id_block_path, regex_path):
    with open(id_block_path, 'r') as tokens_file, open(regex_path, 'r') as regex_file:
        tokens = [line.strip().split(' ') for line in tokens_file.readlines()]
        regex = [line.strip() for line in regex_file.readlines()]
    return tokens, regex


def reconstruct_id(tokens, regex):
    reconstructed_ids = []
    for t, r in zip(tokens, regex):
        id_str = r
        token_list = t.split()
        for i, token in enumerate(token_list):
            id_str = id_str.replace(f"{{{{T{i+1}}}}}", token)
        reconstructed_ids.append(id_str)
    return reconstructed_ids

def reconstruct_g_from_g_prime(g_prime_array, rules_dict):
    De_g = np.zeros_like(g_prime_array)

    for i in range(g_prime_array.shape[0]):
        for j in range(g_prime_array.shape[1]):
            center = g_prime_array[i, j]
            if i == 0:
                up = 0
            else:
                up = De_g[i - 1, j]
            if j == 0:
                left = 0
            else:
                left = De_g[i, j - 1]
            if i != 0 and j != 0:
                left_up = De_g[i - 1, j - 1]
            else:
                left_up = 0
            matched_rule1 = (up, left_up, left, 32)
            matched_rule2 = (up, left_up, left, 224)
            matched_rule3 = (up, left_up, left, 192)
            matched_rule4 = (up, left_up, left, 64)
            matched_rule5 = (up, left_up, left, 0)

            matched_rules = [matched_rule1, matched_rule2, matched_rule3, matched_rule4, matched_rule5]
            max_freq = -1
            top_rule = None

            for rule in matched_rules:
                freq = rules_dict[rule]
                if freq > max_freq:
                    max_freq = freq
                    top_rule = rule

            De_g[i, j] = top_rule[3] if g_prime_array[i, j] == 1 else center
            matched_rule = (up, left_up, left, De_g[i, j])
            rules_dict[matched_rule] += 1

    return De_g


def reconstruct_base_and_quality(g_prime_img, quality_img):
    bases = []
    qualities = []

    # 创建新的 rules_dict
    rules_dict = defaultdict(int)

    g_prime_array = np.array(g_prime_img)
    bases_array = reconstruct_g_from_g_prime(g_prime_array, rules_dict)

    quality_array = np.array(quality_img)


    for i in range(bases_array.shape[0]):
        base_str = ''.join([gray_to_base[pixel] for pixel in bases_array[i]])
        # 还原质量分数
        quality_scores = [q / 2 for q in quality_array[i]]

        bases.append(base_str)
        qualities.append(quality_scores)
        # print(qualities[0])

    return bases, qualities


# def reconstruct_fastq(id_token_files, id_regex_files, base_files, quality_files, output_fastq_path):
#     ids = []
#     for id_tokens_path, id_regex_path in zip(id_token_files, id_regex_files):
#         tokens, regex = load_id_block(id_tokens_path, id_regex_path)
#         ids.extend(reconstruct_id(tokens, regex))
#
#     rules_dict = defaultdict(int)
#     bases, qualities = reconstruct_base_quality(base_files, quality_files, rules_dict)
#
#     records = []
#     for i in range(len(ids)):
#         seq = Seq(bases[i])
#         record = SeqRecord(seq, id=ids[i], description="")
#         record.letter_annotations["phred_quality"] = qualities[i]
#         records.append(record)
#
#     with open(output_fastq_path, "w") as output_handle:
#         SeqIO.write(records, output_handle, "fastq")

def reconstruct_fastq(output_path, id_block, g_prime_img, quality_img):
    records = []

    # 使用 os.path 确保文件名以 '.fastq' 结尾
    output_path = os.path.splitext(output_path)[0] + '.fastq'

    # 将 id_block 转换为列表
    id_block = list(id_block)

    id_tokens = [item[0] for item in id_block]
    id_regex = [item[1] for item in id_block]

    ids = reconstruct_id(id_tokens, id_regex)

    # 将图像转换为数组
    g_prime_array = np.array(g_prime_img)  # 将 g_prime_img 转换为数组
    quality_array = np.array(quality_img)  # 将 quality_img 转换为数组

    # 使用 reconstruct_base_and_quality 进行重构
    g_prime, quality = reconstruct_base_and_quality(g_prime_array, quality_array)

    # 确保长度一致
    if not (len(ids) == len(g_prime) == len(quality)):
        print("错误：ids、g_prime 和 quality 的长度不一致！")
        return

    for i in range(len(ids)):
        seq = Seq(g_prime[i])
        record = SeqRecord(seq, id=ids[i], description="")

        # 还原质量得分为对应的 ASCII 字符
        record.letter_annotations["phred_quality"] = quality[i]
        records.append(record)

    with open(output_path, 'a+') as output_handle:
        num_written = SeqIO.write(records, output_handle, 'fastq')
        # print(f"写入的记录数量：{num_written}")
        # print(f"FASTQ 文件已保存到：{output_path}")


def main(mode, input_path, output_path, lpaq8_path, save, gr_progress):
    _save = None

    if input_path is None or output_path is None or save is None or mode is None:
        tqdm.write("缺少参数，正在退出进程，，，")
        exit(1)

    if mode != "compress" and mode != "c" and mode != "decompress" and mode != "d":
        tqdm.write("mode不符合格式，正在退出进程，，，")
        exit(1)

    if save != "True" and save != "False" and save is not True and save is not False:
        tqdm.write("save不符合格式，正在退出进程，，，")
        exit(1)
    elif save == "True" or save is True:
        _save = True
    elif save == "False" or save is False:
        _save = False
    else:
        exit(1)

    if mode == "compress" or mode == "c":
        tqdm.write(f"info：开始进行fastq压缩程序，文件路径：{input_path}")
        compress(input_path, output_path, lpaq8_path, _save, gr_progress)
    elif mode == "decompress" or mode == "d":
        tqdm.write(f"info：开始进行fastq解压程序，文件路径：{input_path}")
        decompress(input_path, output_path, lpaq8_path, _save, gr_progress)
    else:
        tqdm.write("错误：指定类型错误")

    # 删除临时文件
    if not _save:
        delete_temp_files(output_path)

def delete_temp_files(output_path):
    temp_dir = output_path  # 使用 output_path 作为临时文件目录
    tqdm.write(f"扫描临时文件所在的目录: {temp_dir}")

    # 列出所有文件，方便调试
    all_files = os.listdir(temp_dir)
    tqdm.write(f"目录下的所有文件: {all_files}")

    # 查找符合条件的临时文件
    temp_files = [
        f for f in all_files
        if f.startswith("temp_input") or f.startswith("temp_output")
    ]
    tqdm.write(f"找到的临时文件: {temp_files}")

    # 删除临时文件
    for temp_file in temp_files:
        temp_path = os.path.join(temp_dir, temp_file)
        if os.path.exists(temp_path):
            try:
                os.chmod(temp_path, 0o777)  # 修改文件权限，确保可以删除
                os.remove(temp_path)
                tqdm.write(f"已删除临时文件: {temp_path}")
            except Exception as e:
                tqdm.write(f"无法删除文件 {temp_path}，错误: {e}")

    # 删除 back_compressed 和 front_compressed 文件夹
    for folder_name in ["back_compressed", "front_compressed"]:
        folder_path = os.path.join(temp_dir, folder_name)
        if os.path.exists(folder_path):
            try:
                shutil.rmtree(folder_path)
                tqdm.write(f"已删除文件夹: {folder_path}")
            except Exception as e:
                tqdm.write(f"无法删除文件夹 {folder_path}，错误: {e}")

if __name__ == '__main__':
    lpaq8_path = f"{os.getcwd()}/lpaq8"

    Debug = True
    if Debug:
        # 手动debug
        fastq_path = f"{os.getcwd()}/input/SRR554369.fastq"
        output_path = f"{os.getcwd()}/output"
        compressed_path = f"{os.getcwd()}/input/SRR21733577_2.fastq"
        main("compress", compressed_path, output_path, lpaq8_path, True, None)
        # main("decompress", compressed_path, output_path, lpaq8_path, False, None)
        exit(0)

    # 命令行解析
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='fastq compress')

    # 添加参数
    parser.add_argument('--input_path', type=str, required=True, help='input_path')
    parser.add_argument('--output_path', type=str, required=True, help='output_path')
    parser.add_argument('--mode', type=str, required=True, help='mode')
    parser.add_argument('--save', type=str, default='False', help='save (default: False)')

    # 解析参数
    args = parser.parse_args()

    # 执行主函数
    main(args.mode, args.input_path, args.output_path, lpaq8_path, args.save, None)