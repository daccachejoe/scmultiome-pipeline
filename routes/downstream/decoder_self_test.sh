#!/bin/bash
# Synthetic end-to-end check for the decoder branch: generates
# pseudobulk-shaped synthetic data matching the real sample_metadata
# schema and runs it through 03_train_cv.py -> 04_project_to_peak_space.py
# -> 05_identify_stable_peaks.py, verifying the chain's I/O contracts and
# the CV-leakage fix before real (demultiplexed) data exists. Does not
# touch the R steps (pseudobulk/feature engineering) or a real Seurat
# object -- see scripts/downstream/decoder/self_test_e2e.py.

source "$personal_anaconda_path"
conda activate "$decoder_env_name"

python scripts/downstream/decoder/self_test_e2e.py

conda deactivate
