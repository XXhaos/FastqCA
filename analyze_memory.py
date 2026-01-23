#!/usr/bin/env python3
"""
分析内存占用的脚本
"""
import sys

# 估算128MB块的内存占用
block_size_mb = 128
avg_read_bytes = 300  # 150bp * 2 (seq + quality)

reads_per_block = (block_size_mb * 1024 * 1024) // avg_read_bytes
seq_length = 150

print("=" * 60)
print("内存占用详细分析（128MB数据块）")
print("=" * 60)
print(f"\n每个块的reads数量: {reads_per_block:,}")
print(f"每个read的序列长度: {seq_length}bp\n")

# 1. SeqRecord对象的内存
# SeqRecord是Python对象，有巨大的开销
seqrecord_overhead_per_read = 600  # bytes (保守估计)
seqrecord_total_mb = (reads_per_block * seqrecord_overhead_per_read) / (1024 * 1024)

# 2. NumPy数组
base_array_mb = (reads_per_block * seq_length * 1) / (1024 * 1024)  # uint8
quality_array_mb = (reads_per_block * seq_length * 1) / (1024 * 1024)  # uint8
g_prime_mb = (reads_per_block * seq_length * 1) / (1024 * 1024)  # uint8
q_prime_mb = (reads_per_block * seq_length * 1) / (1024 * 1024)  # uint8

# 3. id_block
id_overhead_per_read = 100  # bytes
id_block_mb = (reads_per_block * id_overhead_per_read) / (1024 * 1024)

# 4. rules_dict
rules_dict_mb = 15
rules_dict_q_mb = 15

print("阶段1：读取和解析records")
print("-" * 60)
print(f"  SeqRecord对象（Python开销）:  {seqrecord_total_mb:>8.1f} MB")
print(f"  合计:                        {seqrecord_total_mb:>8.1f} MB")

print("\n阶段2：process_records()中的峰值内存")
print("-" * 60)
print(f"  records（输入）:             {seqrecord_total_mb:>8.1f} MB")
print(f"  base_image_block:            {base_array_mb:>8.1f} MB")
print(f"  quality_block:               {quality_array_mb:>8.1f} MB")
print(f"  g_prime（生成中）:            {g_prime_mb:>8.1f} MB")
print(f"  q_prime（生成中）:            {q_prime_mb:>8.1f} MB")
print(f"  id_block:                    {id_block_mb:>8.1f} MB")
print(f"  rules_dict × 2:              {rules_dict_mb + rules_dict_q_mb:>8.1f} MB")
peak_process_records = seqrecord_total_mb + base_array_mb + quality_array_mb + g_prime_mb + q_prime_mb + id_block_mb + rules_dict_mb + rules_dict_q_mb
print(f"  峰值合计:                    {peak_process_records:>8.1f} MB")

print("\n阶段3：删除records后")
print("-" * 60)
print(f"  g_prime:                     {g_prime_mb:>8.1f} MB")
print(f"  q_prime:                     {q_prime_mb:>8.1f} MB")
print(f"  id_block:                    {id_block_mb:>8.1f} MB")
after_delete = g_prime_mb + q_prime_mb + id_block_mb
print(f"  合计:                        {after_delete:>8.1f} MB")

print("\n阶段4：back_compress_worker()中")
print("-" * 60)
print(f"  g_prime:                     {g_prime_mb:>8.1f} MB")
print(f"  q_prime:                     {q_prime_mb:>8.1f} MB")
print(f"  id_block:                    {id_block_mb:>8.1f} MB")
print(f"  .npy临时文件（磁盘）:        {g_prime_mb:>8.1f} MB")
print(f"  lpaq8进程（估算）:           ~250.0 MB")
npy_file_mb = g_prime_mb + 0.001  # .npy header overhead
back_compress_peak = g_prime_mb + q_prime_mb + id_block_mb + 250
print(f"  峰值合计:                    {back_compress_peak:>8.1f} MB")

print("\n" + "=" * 60)
print("理论最小内存占用:")
print(f"  单个块处理（无lpaq8）:       {after_delete:>8.1f} MB")
print(f"  单个块处理（含lpaq8）:       {back_compress_peak:>8.1f} MB")
print("=" * 60)

print("\n" + "=" * 60)
print("实际观察到的内存: 2900 MB/进程")
print(f"理论最大峰值:     {peak_process_records:.1f} MB/进程")
print(f"差额:             {2900 - peak_process_records:.1f} MB")
print("=" * 60)

print("\n可能的额外内存来源:")
print("1. Python解释器基础开销: ~100-200 MB")
print("2. BioPython库开销: ~50-100 MB")
print("3. 内存碎片: ~100-200 MB")
print("4. lpaq8子进程实际消耗可能更高: ~500-800 MB")
print("5. 多个块重叠处理时的残留内存")
print("6. SeqRecord对象实际开销可能更大（包含字符串、字典等）")

print("\n优化建议:")
print("=" * 60)
print("1. 使用tobytes()代替np.save()，减少.npy格式开销")
print("2. 减少lpaq8的内存消耗（可能需要调整压缩参数）")
print("3. 更激进的内存释放策略")
print("4. 考虑使用memmap或共享内存")
