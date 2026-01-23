# 内存优化方案

## 当前问题分析
- 理论峰值：585 MB/进程
- 实际观察：2900 MB/进程  
- **差距：2315 MB**

## 主要原因
1. **lpaq8内存消耗（最大问题）**：PAQ系列压缩器使用大量内存
   - lpaq8理论需要256MB，实际可能500-800MB
   - 每个worker同时处理base和quality，峰值时有2个lpaq8进程
   - 估算：每个worker 1000-1500MB用于lpaq8

2. **SeqRecord对象开销**：Python对象比原始数据大3-5倍
   - 估算256MB，实际可能400-500MB

3. **内存碎片**：频繁创建/删除大对象导致碎片
   - 估算：200-300MB

4. **Python解释器和库开销**：约200-300MB

## 优化方案（按效果排序）

### 方案1：降低lpaq8压缩级别 ⭐⭐⭐⭐⭐
**预期效果：减少 800-1200 MB/进程**

修改compress_file()调用，将压缩级别从'9'改为'6'或'7'：
- 级别9：~256-500MB内存
- 级别7：~64-128MB内存  
- 级别6：~32-64MB内存

**代价：压缩率降低5-10%，但仍然很好**

### 方案2：使用更快的压缩器替代lpaq8部分场景
**预期效果：减少 500-800 MB/进程**

对于id_regex和id_tokens（文本数据），使用zstd代替lpaq8：
- zstd内存使用：~10-50MB
- lpaq8内存使用：~256-500MB
- 压缩速度更快
- 压缩率损失很小（文本数据）

**只对g_prime和q_prime使用lpaq8（图像数据）**

### 方案3：序列化处理而不是批处理
**预期效果：减少 300-500 MB/进程**

在back_compress_worker()中，不要同时持有g_block和q_block：
1. 处理并压缩g_block
2. 删除g_block，gc.collect()
3. 再处理q_block

### 方案4：减少block_size（需要用户同意）
**预期效果：减少 200-300 MB/进程**

从128MB减少到64MB：
- records开销：256MB → 128MB
- 数组开销：128MB → 64MB
- 但会增加块数量和文件切换开销

### 方案5：使用tobytes()代替np.save()
**预期效果：减少 20-50 MB/进程**

使用array.tobytes()直接获取原始字节：
```python
# 代替 np.save(temp_input_path, g_block)
with open(temp_input_path, 'wb') as f:
    # 写入shape和dtype信息
    f.write(struct.pack('<II', *g_block.shape))
    f.write(g_block.tobytes())
```

## 推荐实施顺序
1. **先尝试方案1**（降低lpaq8级别到7）- 最简单，效果最好
2. **然后方案3**（序列化处理）- 改动小，效果明显  
3. **考虑方案2**（部分使用zstd）- 需要更多测试
4. 如果还不够，考虑方案4或5

