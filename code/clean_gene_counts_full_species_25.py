#!/usr/bin/env python3
"""
Clean/filter gene_count_matrix.csv files in data/full-species-25 and write gene_counts_cleaned.csv
into each species' cleaned_ subdirectory, aligning sample columns to the existing cleaned modalities.

Behavior (matches prior cleaning approach for genes):
- Load {Species}-gene_count_matrix.csv (expects first column 'gene_id')
- Convert counts to numeric, drop rows with all-zero counts
- Align samples to intersection of standardized sample IDs already present in
  cleaned_<species> (lncrna/mirna/wgbs). Columns are sorted alphabetically for consistency.
- Write data/full-species-25/cleaned_<species>/gene_counts_cleaned.csv
- Emit a simple summary file gene_counts_summary.txt and (re)generate combined_summary.txt

Usage:
  python code/clean_gene_counts_full_species_25.py [--cv-threshold 0.0]

Notes:
- CV filtering is available but disabled by default (kept for parity with prior script).
"""

import os
import re
import argparse
from typing import List, Dict

import numpy as np
import pandas as pd


# Standard sample ID like ACR-139-TP1
SAMPLE_STD_PATTERN = re.compile(r"[A-Z]+-\d+-TP\d+")


def workspace_root() -> str:
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def to_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def filter_zero_rows(df: pd.DataFrame) -> pd.DataFrame:
    mask = (df.fillna(0).sum(axis=1) != 0)
    return df.loc[mask]


def cv(series: pd.Series) -> float:
    arr = series.values.astype(float)
    mean = np.nanmean(arr)
    std = np.nanstd(arr)
    if mean == 0:
        return 0.0
    return float(std / mean)


def filter_low_variation(df: pd.DataFrame, cv_threshold: float | None) -> pd.DataFrame:
    if not cv_threshold or cv_threshold <= 0:
        return df
    cvs = df.apply(cv, axis=1)
    return df.loc[cvs >= cv_threshold]


def summarize_df(name: str, df: pd.DataFrame) -> str:
    n_features, n_samples = df.shape
    sparsity = float((df == 0).sum().sum()) / float(n_features * n_samples) if n_features * n_samples else 0.0
    return (
        f"Dataset: {name}\n"
        f"Features: {n_features}\n"
        f"Samples: {n_samples}\n"
        f"Sparsity: {sparsity*100:.2f}%\n"
    )


def load_gene(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Ensure first column is named gene_id
    if df.columns[0].lower() != "gene_id":
        df = df.rename(columns={df.columns[0]: "gene_id"})
    df = df.set_index("gene_id")
    return to_numeric_df(df)


def standardized_columns(cols: List[str]) -> List[str]:
    return [c for c in cols if SAMPLE_STD_PATTERN.fullmatch(c)]


def align_to_cleaned_modalities(species_key: str, df_gene: pd.DataFrame, base_dir: str) -> pd.DataFrame:
    cleaned_dir = os.path.join(base_dir, f"cleaned_{species_key}")
    # Load existing cleaned modalities to determine common standardized samples
    paths = {
        "lncrna": os.path.join(cleaned_dir, "lncrna_counts_cleaned.csv"),
        "mirna": os.path.join(cleaned_dir, "mirna_counts_cleaned.csv"),
        "wgbs": os.path.join(cleaned_dir, "wgbs_counts_cleaned.csv"),
    }

    sample_sets: List[set[str]] = []
    for p in paths.values():
        if os.path.exists(p):
            df = pd.read_csv(p, nrows=1)  # fast header peek
            cols = [c for c in df.columns if c != df.columns[0]]
            sample_sets.append(set(standardized_columns(cols)))

    # If none found, fall back to standardized gene columns only
    gene_std = set(standardized_columns(list(df_gene.columns)))
    common = set.intersection(*sample_sets) if sample_sets else gene_std
    if not common:
        # Nothing matches pattern; keep gene as-is
        common = gene_std if gene_std else set(df_gene.columns)

    ordered = sorted(common)
    return df_gene.loc[:, ordered]


def write_outputs(species_key: str, df_gene: pd.DataFrame, base_dir: str) -> None:
    cleaned_dir = os.path.join(base_dir, f"cleaned_{species_key}")
    os.makedirs(cleaned_dir, exist_ok=True)

    # Write cleaned gene counts
    out_gene = df_gene.copy()
    out_gene.index.name = "gene_id"
    gene_csv_path = os.path.join(cleaned_dir, "gene_counts_cleaned.csv")
    out_gene.to_csv(gene_csv_path)

    # Write summary for genes
    gene_summary = summarize_df("gene", out_gene)
    with open(os.path.join(cleaned_dir, "gene_counts_summary.txt"), "w") as f:
        f.write(gene_summary)

    # Rebuild combined summary including other modalities if present
    summaries: Dict[str, str] = {"gene": gene_summary}

    def maybe_summarize(kind: str, filename: str, index_name: str):
        p = os.path.join(cleaned_dir, filename)
        if os.path.exists(p):
            df = pd.read_csv(p)
            # first column is index
            first = df.columns[0]
            df = df.set_index(first)
            summaries[kind] = summarize_df(kind, df)

    maybe_summarize("lncrna", "lncrna_counts_cleaned.csv", "lncrna_id")
    maybe_summarize("mirna", "mirna_counts_cleaned.csv", "mirna_id")
    maybe_summarize("wgbs", "wgbs_counts_cleaned.csv", "cpg_id")

    with open(os.path.join(cleaned_dir, "combined_summary.txt"), "w") as f:
        f.write("Cleaned datasets summary\n\n")
        for k in ["gene", "lncrna", "mirna", "wgbs"]:
            if k in summaries:
                f.write(summaries[k] + "\n")


def main():
    parser = argparse.ArgumentParser(description="Clean gene_count_matrix files in full-species-25.")
    parser.add_argument("--cv-threshold", type=float, default=None, help="Optional coefficient of variation threshold for filtering")
    args = parser.parse_args()

    root = workspace_root()
    base_dir = os.path.join(root, "data", "full-species-25")
    assert os.path.isdir(base_dir), f"Directory not found: {base_dir}"

    # Discover species by presence of *-gene_count_matrix.csv
    entries = [f for f in os.listdir(base_dir) if f.endswith("-gene_count_matrix.csv")]
    if not entries:
        raise SystemExit("No gene_count_matrix.csv files found in full-species-25")

    for entry in sorted(entries):
        species = entry.split("-", 1)[0]  # e.g., Apul
        species_key = species.lower()
        gene_path = os.path.join(base_dir, entry)

        print(f"Processing {species} -> {gene_path}")
        gene = load_gene(gene_path)
        print(" - drop all-zero rows")
        gene = filter_zero_rows(gene)
        if args.cv_threshold and args.cv_threshold > 0:
            print(f" - CV filter @ {args.cv_threshold}")
            gene = filter_low_variation(gene, args.cv_threshold)

        print(" - align samples to cleaned modalities")
        gene = align_to_cleaned_modalities(species_key, gene, base_dir)
        # Sort columns for consistency
        gene = gene.loc[:, sorted(gene.columns)]

        print(" - write outputs")
        write_outputs(species_key, gene, base_dir)

    print("Done.")


if __name__ == "__main__":
    main()
