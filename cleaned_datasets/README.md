# Cleaned and Standardized Datasets

This directory contains cleaned and standardized datasets that are ready for correlation analysis and other downstream analyses.

## Overview

All datasets have been processed to ensure:
- **Common sample IDs**: All datasets use the same 40 sample identifiers
- **No zero expression**: Features with zero expression across all samples have been removed
- **Sufficient variation**: Features with limited variation (CV < 0.1) have been removed
- **Consistent structure**: All datasets have the same column order and sample alignment

## Dataset Summary

| Dataset | Features | Samples | Sparsity | Description |
|---------|----------|---------|----------|-------------|
| `gene_counts_cleaned.csv` | 36,084 | 40 | 37.8% | Gene expression counts |
| `lncrna_counts_cleaned.csv` | 15,900 | 40 | 3.8% | Long non-coding RNA expression counts |
| `mirna_counts_cleaned.csv` | 51 | 40 | 7.8% | MicroRNA expression counts |
| `wgbs_counts_cleaned.csv` | 249 | 40 | 41.4% | WGBS CpG methylation counts |

## Sample Structure

All datasets contain the same 40 samples representing different time points (TP1-TP4) across 10 different conditions (ACR-139 through ACR-265):

- **ACR-139**: TP1, TP2, TP3, TP4
- **ACR-145**: TP1, TP2, TP3, TP4
- **ACR-150**: TP1, TP2, TP3, TP4
- **ACR-173**: TP1, TP2, TP3, TP4
- **ACR-186**: TP1, TP2, TP3, TP4
- **ACR-225**: TP1, TP2, TP3, TP4
- **ACR-229**: TP1, TP2, TP3, TP4
- **ACR-237**: TP1, TP2, TP3, TP4
- **ACR-244**: TP1, TP2, TP3, TP4
- **ACR-265**: TP1, TP2, TP3, TP4

## Data Quality Metrics

### Zero Expression Filtering
- **Gene counts**: Removed 8,287 features (18.7% of original)
- **lncRNA counts**: No features removed (already clean)
- **miRNA counts**: No features removed (already clean)
- **WGBS counts**: Removed 2 features (0.8% of original)

### Variation Filtering
- **Minimum CV threshold**: 0.1 (10% coefficient of variation)
- **All datasets**: Passed variation filtering with no features removed

### Data Types
- **Gene counts**: Integer counts
- **lncRNA counts**: Float values (normalized)
- **miRNA counts**: Integer counts
- **WGBS counts**: Float values (methylation percentages)

## File Descriptions

### Main Data Files
- `*_cleaned.csv`: Cleaned and standardized datasets ready for analysis
- `*_summary.txt`: Individual dataset summaries with statistics
- `combined_summary.txt`: Overall summary of all datasets

### Original Data Sources
- **Gene counts**: `data/apul-gene_count_matrix.csv`
- **lncRNA counts**: `data/lncrna_counts_cleaned.csv`
- **miRNA counts**: `data/Apul_miRNA_counts_formatted.txt`
- **WGBS counts**: `data/merged-WGBS-CpG-counts_filtered.csv`

## Usage

These cleaned datasets are ready for:
- Correlation analysis between different data types
- Time series analysis across TP1-TP4 time points
- Multi-omics integration studies
- Statistical modeling and machine learning

## Processing Scripts

The datasets were created using:
- `code/clean_and_standardize_datasets.py`: Main cleaning and standardization script
- `code/verify_cleaned_datasets.py`: Verification script to confirm data quality

## Notes

- All datasets are properly aligned by sample ID
- No missing values in any dataset
- Expression values are preserved (no normalization applied)
- Sample order is consistent across all datasets
- Ready for immediate use in correlation analysis workflows
