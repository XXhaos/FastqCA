# FastqCA

FastqCA is a FASTQ compression/decompression tool that supports both lossy and lossless modes, with multi-threaded, block-based processing. The entry point is `main_new.py`.

## Features

- **Lossy**: Quantizes quality scores and compresses bases/qualities for smaller archives.
- **LossLess**: Preserves original data for fully lossless compression.
- **Multi-threaded + block-based**: Tune performance via thread count and block size.
- **Read count manifest**: Lossy mode records read counts and validates them on restore.

## Dependencies

Use Python 3.8+ and install the dependencies below:

```bash
pip install biopython numpy pillow tqdm
```

## Quick start

### Lossy compression

```bash
python main_new.py \
  --compressor Lossy \
  --input_path /path/to/input.fastq \
  --output_path /path/to/output \
  --mode compress
```

### Lossy decompression

```bash
python main_new.py \
  --compressor Lossy \
  --input_path /path/to/output \
  --output_path /path/to/restore \
  --mode decompress
```

> Lossy mode generates a `.readcount` file during compression to validate read counts on restore. Use `--manifest` to override the path.

### Lossless compression

```bash
python main_new.py \
  --compressor LossLess \
  --input_path /path/to/input.fastq \
  --output_path /path/to/output \
  --mode compress
```

### Lossless decompression

```bash
python main_new.py \
  --compressor LossLess \
  --input_path /path/to/output \
  --output_path /path/to/restore \
  --mode decompress
```

## Common options

| Option | Description | Default |
| --- | --- | --- |
| `--compressor` | Select compressor: `Lossy` or `LossLess` | `Lossy` |
| `--input_path` | Input FASTQ (compress) or archive path (decompress) | Required |
| `--output_path` | Output path | Required |
| `--mode` | `compress`/`c` or `decompress`/`d` | Required |
| `--save` | Keep intermediate files (`True`/`False`) | `False` |
| `--threads` | Worker threads | `os.cpu_count()` |
| `--block_size` | Block size in bytes | `134217728` |
| `--manifest` | Lossy read-count manifest path | `*.readcount` |

## Notes

- The repository includes `lpaq8` binaries/scripts used by the backend compression steps.
- Lossy decompression checks read counts and raises an error if the restored count differs.
