# C. elegans Test Dataset for ConTra Validation

## Overview

This dataset contains synthetic multi-omics data modeled after C. elegans developmental 
time series for validating the ConTra analysis pipeline. The data includes known 
regulatory interactions that can be used to assess pipeline performance.

## Dataset Contents

- `gene_counts_cleaned.csv`: Gene expression data (25 genes × 8 samples)
- `lncrna_counts_cleaned.csv`: lncRNA expression data (20 lncRNAs × 8 samples)  
- `mirna_counts_cleaned.csv`: miRNA expression data (6 miRNAs × 8 samples)
- `wgbs_counts_cleaned.csv`: DNA methylation data (50 CpG sites × 8 samples)
- `validation_metadata.json`: Known interactions and validation criteria
- `*_summary.txt`: Individual dataset statistics

## Sample Structure

Samples represent C. elegans developmental stages:
- L1-TP1, L1-TP2: L1 larval stage (early and late)
- L2-TP1, L2-TP2: L2 larval stage (early and late)
- L3-TP1, L3-TP2: L3 larval stage (early and late)  
- L4-TP1, L4-TP2: L4 larval stage (early and late)

## Known Regulatory Interactions

### let-7 Pathway (L3-L4 stages)
- **miRNA**: let-7 (upregulated in L3-L4)
- **Targets**: lin-41, hbl-1 (downregulated in L3-L4)
- **Expected**: Negative correlation, stage-specific

### lin-4 Cascade (L1-L2 stages)  
- **miRNA**: lin-4 (upregulated in L1-L2)
- **Targets**: lin-14, lin-28 (downregulated in L2)
- **Expected**: Negative correlation, stage-specific

## Usage with ConTra

```bash
# Copy test data to main data directory
cp -r data/test_datasets/celegans_modENCODE/* data/cleaned_datasets/

# Run ConTra analysis
cd code
python subset_context_dependent_analysis.py

# Validate results against known interactions
python validate_celegans_results.py
```

## Validation Criteria

A successful validation should detect:
1. let-7 negative regulation of lin-41, hbl-1 in L3-L4 stages
2. lin-4 negative regulation of lin-14, lin-28 in L1-L2 stages  
3. Context-dependent interactions matching developmental timing
4. Appropriate statistical significance (p < 0.05)

## References

- C. elegans developmental timing: Slack & Ruvkun (1997) Cell
- let-7 miRNA function: Reinhart et al. (2000) Nature
- modENCODE project: modENCODE Consortium (2010) Science
