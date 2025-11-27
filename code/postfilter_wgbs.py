#!/usr/bin/env python3
"""
Post-filter WGBS cleaned datasets across species-specific cleaned_* directories.

For each directory matching data/full-species-24/cleaned_*/wgbs_counts_cleaned.csv:
 1. Load file.
 2. Remove any rows with missing (NaN) values across samples.
 3. Remove rows where no value exceeds a threshold (default > 10).
 4. Write a new file wgbs_counts_cleaned_postfiltered.csv alongside original.
 5. Produce a wgbs_counts_postfilter_summary.txt with stats.
Optionally overwrite original (--overwrite) instead of creating new file.

Usage:
  python code/postfilter_wgbs.py --base-dir data/full-species-24 --threshold 10
  python code/postfilter_wgbs.py --base-dir data/full-species-24 --threshold 10 --overwrite
"""

import argparse
import os
import pandas as pd
from typing import Tuple


def summarize(df: pd.DataFrame) -> str:
    n_features, n_samples = df.shape
    sparsity = (df == 0).sum().sum() / (n_features * n_samples) if n_features and n_samples else 0
    return (
        f"Features: {n_features}\n"
        f"Samples: {n_samples}\n"
        f"Sparsity: {sparsity*100:.2f}%\n"
    )


def process_file(path: str, threshold: float, overwrite: bool) -> Tuple[str, str]:
    df = pd.read_csv(path, index_col=0)
    before = df.shape[0]
    # Drop rows with any NaNs
    df = df.dropna(axis=0, how='any')
    after_dropna = df.shape[0]
    # Keep rows where max value > threshold
    keep_mask = (df.max(axis=1) > threshold)
    df_filtered = df.loc[keep_mask]
    after_threshold = df_filtered.shape[0]

    dir_name = os.path.dirname(path)
    out_name = 'wgbs_counts_cleaned_postfiltered.csv'
    out_path = os.path.join(dir_name, out_name)

    if overwrite:
        out_path = path  # Overwrite original
    df_filtered.to_csv(out_path)

    summary_lines = [
        f"Original rows: {before}",
        f"Rows after dropna: {after_dropna}",
        f"Rows after threshold (> {threshold}): {after_threshold}",
        summarize(df_filtered).strip()
    ]
    summary_text = "\n".join(summary_lines) + "\n"
    with open(os.path.join(dir_name, 'wgbs_counts_postfilter_summary.txt'), 'w') as f:
        f.write(summary_text)
    return out_path, summary_text


def main():
    parser = argparse.ArgumentParser(description="Post-filter WGBS cleaned datasets across species cleaned directories.")
    parser.add_argument('--base-dir', default=os.path.join('data', 'full-species-24'), help='Base directory containing cleaned_* subfolders')
    parser.add_argument('--threshold', type=float, default=10.0, help='Value threshold; require at least one value above this (default 10)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite original wgbs_counts_cleaned.csv instead of writing a new file')
    args = parser.parse_args()

    workspace_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base_dir = os.path.join(workspace_root, args.base_dir)
    assert os.path.isdir(base_dir), f"Base directory not found: {base_dir}"

    targets = []
    for entry in os.listdir(base_dir):
        if entry.startswith('cleaned_'):
            candidate = os.path.join(base_dir, entry, 'wgbs_counts_cleaned.csv')
            if os.path.isfile(candidate):
                targets.append(candidate)

    if not targets:
        print('No wgbs_counts_cleaned.csv files found.')
        return

    print(f"Found {len(targets)} WGBS cleaned files to post-filter.")
    for t in targets:
        print(f"Processing {t} ...")
        out_path, summary = process_file(t, args.threshold, args.overwrite)
        print(f"  -> Written: {out_path}")
        print("  Summary:\n" + summary)

    print("Done post-filtering.")


if __name__ == '__main__':
    main()
