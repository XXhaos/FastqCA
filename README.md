# FastqCA

FastqCA is a cellular-automaton-inspired compressor for FASTQ files. It turns read identifiers, base calls, and quality scores into image-like tensors, learns local patterns, and feeds the transformed blocks to the lightweight `lpaq8` backend for final entropy coding. Both lossless and lossy (quality-quantized) pipelines are available.

## Features
- **Lossless and lossy modes**: choose strict reconstruction or quantized qualities (Q-score buckets 5/12/18/24) for smaller outputs.
- **Block-wise processing**: splits input FASTQ files into configurable chunks (default 128 MB) to keep memory usage stable.
- **Two-stage compression**: front-end converts reads to predicted base/quality matrices plus ID tokens; back-end packs each part with `lpaq8`.
- **Optional intermediates**: save TIFFs and token text files for inspection or debugging.

## Requirements
- Python 3.9+
- Dependencies: [`biopython`](https://biopython.org/), [`pillow`](https://python-pillow.org/), [`numpy`](https://numpy.org/), [`tqdm`](https://tqdm.github.io/)
- `lpaq8` binary in the repository root (provided) or build it yourself:

  ```bash
  g++ -O2 lpaq8.cpp -o lpaq8
  ```

Install Python dependencies (virtual environment recommended):

```bash
pip install biopython pillow numpy tqdm
```

## Command-line usage
Use `main_new.py` to run either compressor. The `--compressor` flag selects the lossy or lossless pipeline, and `--mode` toggles between compression and decompression.

```bash
python main_new.py \
  --compressor Lossless \
  --input_path path/to/reads.fastq \
  --output_path path/to/output_dir \
  --mode c \
  --save False
```

### Key arguments
- `--compressor`: `Lossy`/`lossy` or `LossLess`/`lossless` (case-insensitive variants supported by the script).
- `--mode`: `c`/`compress` to create a compressed blob; `d`/`decompress` to reconstruct a FASTQ file.
- `--input_path`: FASTQ file to compress, or the previously produced compressed blob when decompressing.
- `--output_path`: directory or file path for results. When compressing to a directory, the output file is named after the input FASTQ stem.
- `--save`: `True` keeps intermediate TIFFs and token text files under `front_compressed/` and `back_compressed/`; `False` cleans them up after finishing.

### Examples
- **Lossless compression**

  ```bash
  python main_new.py --compressor Lossless --input_path input/sample.fastq --output_path output --mode c
  ```

- **Lossy compression (quantized qualities)**

  ```bash
  python main_new.py --compressor Lossy --input_path input/sample.fastq --output_path output --mode compress --save True
  ```

- **Decompression**

  ```bash
  python main_new.py --compressor Lossless --input_path output/sample --output_path output --mode d
  ```

## Outputs
- The compressed blob lives at `OUTPUT_DIR/<input-stem>` (no extension) and contains interleaved `lpaq8` segments for IDs, bases, and qualities.
- When `--save True`, intermediate artifacts are stored alongside the output:
  - `front_compressed/chunk_<n>_base.tiff` and `front_compressed/chunk_<n>_quality.tiff`
  - `front_compressed/chunk_<n>_id_tokens.txt` and `front_compressed/chunk_<n>_id_regex.txt`
  - `back_compressed/*.lpaq8` for each component

## Convenience script
`main.sh` shows a sample lossless invocation you can adapt for quick testing:

```bash
python main_new.py --compressor "Lossless" --input_path "input/SRR554369.fastq" --output_path "output" --mode "c"
```

## Notes
- Default block size is 128 MB; adjust `block_size` in `Lossy.py` or `LossLess.py` if you need different memory characteristics.
- Temporary files prefixed with `temp_input`/`temp_output` are removed automatically unless `--save True`.
