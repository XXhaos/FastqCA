# FastqCA Project Overview

This document summarizes the layout and core logic of the FastqCA repository after an initial review.

## Entry Points and CLI
- `main.py` is the primary CLI wrapper. It accepts `--compressor` (Lossy/LossLess), `--input_path`, `--output_path`, `--mode`, and an optional `--save` flag before dispatching to the corresponding compressor. The path to the bundled `lpaq8` binary is resolved relative to the current working directory.

## Compression Pipelines
### Lossy compression (`Lossy.py`)
- Maps nucleotide bases to grayscale values and quality scores through the `Q4` quantization step.
- Splits FASTQ records into blocks of reads, builds base/quality matrices, and predicts pixels using cellular automaton–style rule dictionaries (`generate_g_prime` and `generate_q_prime`).
- Performs a **front-end** transform that outputs TIFF images for bases/qualities plus identifier tokenization metadata (`front_compress`). Optional saves write intermediate TIFFs and token files under `front_compressed/`.
- A **back-end** step serializes the transformed artifacts into a single stream, compressing each section with the external `lpaq8` binary while tracking progress (`back_compress`).

### Lossless compression (`LossLess.py`)
- Shares the same FASTQ parsing, identifier tokenization, and grayscale mapping as the lossy path, but skips quality quantization and retains full information for exact reconstruction.
- The pipeline is similarly split into front-end block preparation and back-end compression using `lpaq8` for entropy coding.

## lpaq8 Integration
- `lpaq8.py` provides thin wrappers around the bundled `lpaq8` executable for compression/decompression. Functions spawn subprocesses and are reused by both lossy and lossless workflows (`compress_file` / `decompress_file`).

## Notable Dependencies
- BioPython (`Bio.SeqIO`, `SeqRecord`) for FASTQ parsing.
- Pillow for image creation from numpy arrays.
- `tqdm` for progress reporting during both preprocessing and lpaq8 monitoring.

## Repository Contents
- `Lossy.py`, `LossLess.py`, and `LossLess_thread.py` implement compression logic.
- `lpaq8`/`lpaq8.exe` binaries and `lpaq8.cpp` source provide the entropy coder.
- Helper scripts such as `tools.py` and shell wrappers (`main.sh`) support running the compressors.
