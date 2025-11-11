#!/usr/bin/env python3
"""
Clean and standardize raw multi-omics data in data/full-species-24 for a chosen species.

Outputs a new cleaned folder inside the species directory, producing four files:
- gene_counts_cleaned.csv
- lncrna_counts_cleaned.csv
- mirna_counts_cleaned.csv
- wgbs_counts_cleaned.csv
and per-dataset summaries plus a combined summary.

Usage:
  python code/clean_full_species_data.py --species apul --source-dir data/full-species-24 --out-subdir cleaned_apul

Notes:
- Harmonizes sample IDs across modalities to the format ACR-###-TPx
- Drops features with all-zero values
- Optionally drops low-variation features using coefficient of variation threshold (default disabled)
- Ensures consistent column order and alignment across all modalities
"""

import argparse
import os
import re
from typing import Dict, List
import pandas as pd
import numpy as np

# ---------------------------
# Helpers
# ---------------------------

# Standard sample ID like ACR-139-TP1, POR-216-TP2, POC-201-TP3, etc.
SAMPLE_STD_PATTERN = re.compile(r"[A-Z]+-\d+-TP\d+")
# miRNA headers like 1A11_POR-216_TP2 or 1A10_ACR-145_TP4 -> capture group code and TP
MIRNA_COL_PATTERN = re.compile(r".*_([A-Za-z]+-\d+)_TP(\d+)")

SPECIES_KEYS = {
    "apul": {
        "gene": "Component_24_apul_Gene_count_matrix.csv",
        "lncrna": "Apul-lncRNA_counts.clean.filtered.txt",
        "mirna": "Apul_miRNA_counts_formatted.txt",
        "wgbs": "Apul-merged-WGBS-CpG-counts_filtered_n20.csv",
    },
    "peve": {
        "gene": "Component_24_peve_Gene_count_matrix.csv",
        "lncrna": "Peve-lncRNA_counts.clean.filtered.txt",
        "mirna": "Peve_miRNA_counts_formatted.txt",
        "wgbs": "Peve-merged-WGBS-CpG-counts_filtered_n20.csv",
    },
    "ptua": {
        "gene": "Component_24_ptua_Gene_count_matrix.csv",
        "lncrna": "Ptua-lncRNA_counts.clean.filtered.txt",
        "mirna": "Ptua_miRNA_counts_formatted.txt",
        "wgbs": "Ptua-merged-WGBS-CpG-counts_filtered_n20.csv",
    },
}


def _ensure_workspace_root() -> str:
    # code/clean_full_species_data.py -> repo_root
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _to_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _cv(series: pd.Series) -> float:
    arr = series.values.astype(float)
    mean = np.nanmean(arr)
    std = np.nanstd(arr)
    if mean == 0:
        return 0.0
    return float(std / mean)


def _summarize(name: str, df: pd.DataFrame) -> str:
    n_features, n_samples = df.shape
    sparsity = float((df == 0).sum().sum()) / float(n_features * n_samples) if n_features * n_samples else 0.0
    return (
        f"Dataset: {name}\n"
        f"Features: {n_features}\n"
        f"Samples: {n_samples}\n"
        f"Sparsity: {sparsity*100:.2f}%\n"
    )


def _write_summary(out_dir: str, summaries: Dict[str, str]) -> None:
    lines = ["Cleaned datasets summary\n", "\n"]
    for k, v in summaries.items():
        lines.append(v + "\n")
    with open(os.path.join(out_dir, "combined_summary.txt"), "w") as f:
        f.write("".join(lines))


def harmonize_mirna_columns(cols: List[str]) -> List[str]:
    out = []
    for c in cols:
        m = MIRNA_COL_PATTERN.match(c)
        if m:
            out.append(f"{m.group(1)}-TP{m.group(2)}")
        else:
            # Sometimes miRNA files may already be standardized, keep if matches
            if SAMPLE_STD_PATTERN.fullmatch(c):
                out.append(c)
            else:
                out.append(c)
    return out


def load_gene(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.columns[0].lower() != "gene_id":
        df = df.rename(columns={df.columns[0]: "gene_id"})
    df = df.set_index("gene_id")
    df = _to_numeric_df(df)
    return df


def load_lncrna(path: str) -> pd.DataFrame:
    # lncRNA file has meta columns, keep only sample columns
    df = pd.read_csv(path, sep="\t")
    # Expect first column is feature id
    feature_col = df.columns[0]
    # Drop typical meta columns if present
    drop_cols = ["Chr", "Start", "End", "Strand", "Length"]
    keep_cols = [feature_col] + [c for c in df.columns if c not in drop_cols and c != feature_col]
    df = df[keep_cols]
    df = df.rename(columns={feature_col: "lncrna_id"}).set_index("lncrna_id")
    df = _to_numeric_df(df)
    return df


def load_mirna(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # First column is feature id, remainder are samples that need renaming
    feature_col = df.columns[0]
    sample_cols = list(df.columns[1:])
    std_cols = harmonize_mirna_columns(sample_cols)
    col_map = {old: new for old, new in zip(sample_cols, std_cols)}
    df = df.rename(columns=col_map)
    df = df.rename(columns={feature_col: "mirna_id"}).set_index("mirna_id")
    df = _to_numeric_df(df)
    return df


def load_wgbs(path: str) -> pd.DataFrame:
    # Large file; use pandas default csv reader. Assume first column is CpG id (e.g., 'CpG').
    df = pd.read_csv(path)
    first_col = df.columns[0]
    # Set index using the original first column name to avoid case-sensitivity issues
    df = df.set_index(first_col)
    df = _to_numeric_df(df)
    return df


def filter_zero_rows(df: pd.DataFrame) -> pd.DataFrame:
    mask = (df.fillna(0).sum(axis=1) != 0)
    return df.loc[mask]


def filter_low_variation(df: pd.DataFrame, cv_threshold: float = None) -> pd.DataFrame:
    if cv_threshold is None or cv_threshold <= 0:
        return df
    cvs = df.apply(_cv, axis=1)
    return df.loc[cvs >= cv_threshold]


def align_samples(dfs: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    # Keep only standardized sample names across all modalities
    sample_sets = [set([c for c in df.columns if SAMPLE_STD_PATTERN.fullmatch(c)]) for df in dfs.values()]
    common = set.intersection(*sample_sets)
    if not common:
        raise ValueError("No common standardized sample IDs across modalities. Check input headers.")
    ordered = sorted(common)
    return {k: v.loc[:, ordered] for k, v in dfs.items()}


def write_outputs(out_dir: str, dfs: Dict[str, pd.DataFrame]) -> Dict[str, str]:
    os.makedirs(out_dir, exist_ok=True)

    summaries: Dict[str, str] = {}

    for key, df in dfs.items():
        out_name = {
            "gene": "gene_counts_cleaned.csv",
            "lncrna": "lncrna_counts_cleaned.csv",
            "mirna": "mirna_counts_cleaned.csv",
            "wgbs": "wgbs_counts_cleaned.csv",
        }[key]
        df_out = df.copy()
        df_out.index.name = {
            "gene": "gene_id",
            "lncrna": "lncrna_id",
            "mirna": "mirna_id",
            "wgbs": "cpg_id",
        }[key]
        df_out.to_csv(os.path.join(out_dir, out_name))

        # Individual summary
        summary = _summarize(key, df_out)
        summaries[key] = summary
        with open(os.path.join(out_dir, f"{key}_counts_summary.txt"), "w") as f:
            f.write(summary)

    # Combined summary
    _write_summary(out_dir, summaries)
    return summaries


def main():
    parser = argparse.ArgumentParser(description="Clean full-species-24 raw data for a species.")
    parser.add_argument("--species", choices=sorted(SPECIES_KEYS.keys()), default="apul", help="Species key to clean (default: apul)")
    parser.add_argument("--source-dir", default=os.path.join("data", "full-species-24"), help="Source directory containing raw files")
    parser.add_argument("--out-subdir", default=None, help="Name of output subfolder inside source-dir (default: cleaned_<species>)")
    parser.add_argument("--cv-threshold", type=float, default=None, help="Optional coefficient of variation threshold for filtering (e.g., 0.1)")

    args = parser.parse_args()

    workspace_root = _ensure_workspace_root()
    source_dir = os.path.join(workspace_root, args.source_dir)
    assert os.path.isdir(source_dir), f"Source directory not found: {source_dir}"

    spec = SPECIES_KEYS[args.species]
    in_paths = {k: os.path.join(source_dir, v) for k, v in spec.items()}

    # Load
    print(f"Loading raw files for {args.species} from {source_dir} ...")
    gene = load_gene(in_paths["gene"])  # CSV
    lncrna = load_lncrna(in_paths["lncrna"])  # TSV with meta columns
    mirna = load_mirna(in_paths["mirna"])  # TSV with sample name normalization
    wgbs = load_wgbs(in_paths["wgbs"])  # large CSV

    # Basic filters
    print("Applying zero-sum row filtering ...")
    gene = filter_zero_rows(gene)
    lncrna = filter_zero_rows(lncrna)
    mirna = filter_zero_rows(mirna)
    wgbs = filter_zero_rows(wgbs)

    if args.cv_threshold is not None and args.cv_threshold > 0:
        print(f"Applying CV filtering with threshold {args.cv_threshold} ...")
        gene = filter_low_variation(gene, args.cv_threshold)
        lncrna = filter_low_variation(lncrna, args.cv_threshold)
        mirna = filter_low_variation(mirna, args.cv_threshold)
        wgbs = filter_low_variation(wgbs, args.cv_threshold)

    # Align samples across modalities
    print("Aligning common samples across modalities ...")
    dfs = align_samples({"gene": gene, "lncrna": lncrna, "mirna": mirna, "wgbs": wgbs})

    # Order columns consistently
    ordered_cols = sorted(dfs["gene"].columns)
    for k in dfs:
        dfs[k] = dfs[k].loc[:, ordered_cols]

    # Output directory
    out_subdir = args.out_subdir or f"cleaned_{args.species}"
    out_dir = os.path.join(source_dir, out_subdir)

    print(f"Writing cleaned outputs to: {out_dir}")
    summaries = write_outputs(out_dir, dfs)

    print("\nSummary:\n" + "\n".join([f"- {k}: {v.splitlines()[1]}{v.splitlines()[2]}{v.splitlines()[3]}" for k, v in summaries.items()]))
    print("\nDone.")


if __name__ == "__main__":
    main()
