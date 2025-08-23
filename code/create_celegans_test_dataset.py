#!/usr/bin/env python3
"""
Download and format C. elegans modENCODE test dataset for ConTra validation.

This script downloads multi-omics data from the modENCODE project for C. elegans
developmental time series and formats it to match ConTra's input requirements.
"""

import os
import sys
import pandas as pd
import numpy as np
import requests
import gzip
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CelegansTestDataDownloader:
    """Download and format C. elegans test data for ConTra validation."""
    
    def __init__(self, output_dir: str = "data/test_datasets/celegans_modENCODE"):
        """Initialize the downloader."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Sample mapping: developmental stages to sample IDs
        self.sample_mapping = {
            'L1-TP1': 'WBls:0000024',  # L1 larva early
            'L1-TP2': 'WBls:0000025',  # L1 larva late
            'L2-TP1': 'WBls:0000026',  # L2 larva early  
            'L2-TP2': 'WBls:0000027',  # L2 larva late
            'L3-TP1': 'WBls:0000028',  # L3 larva early
            'L3-TP2': 'WBls:0000029',  # L3 larva late
            'L4-TP1': 'WBls:0000030',  # L4 larva early
            'L4-TP2': 'WBls:0000031',  # L4 larva late
        }
        
        # Known regulatory interactions for validation
        self.known_interactions = {
            'let-7_pathway': {
                'miRNA': 'let-7',
                'targets': ['lin-41', 'hbl-1', 'daf-12'],
                'stage_specific': ['L3-TP2', 'L4-TP1', 'L4-TP2'],
                'interaction_type': 'negative_regulation'
            },
            'lin-4_cascade': {
                'miRNA': 'lin-4', 
                'targets': ['lin-14', 'lin-28'],
                'stage_specific': ['L1-TP2', 'L2-TP1'],
                'interaction_type': 'negative_regulation'
            },
            'developmental_timing': {
                'lncRNAs': ['lincRNA-1', 'lincRNA-2'],
                'targets': ['cog-1', 'unc-3'],
                'methylation_sites': ['chr1:1000-2000', 'chr2:5000-6000'],
                'stage_specific': 'all'
            }
        }
    
    def create_synthetic_data(self) -> Dict[str, pd.DataFrame]:
        """
        Create synthetic C. elegans multi-omics data that mimics real patterns
        and includes known regulatory interactions for validation.
        
        Note: This creates synthetic data as a placeholder. In a real implementation,
        this would download actual modENCODE data from ENCODE portal.
        """
        logger.info("Creating synthetic C. elegans test dataset...")
        
        # Create sample names matching ConTra format (condition-timepoint)
        samples = []
        for stage in ['L1', 'L2', 'L3', 'L4']:
            for tp in ['TP1', 'TP2']:
                samples.append(f"{stage}-{tp}")
        
        # Gene expression data (simplified subset of C. elegans genes)
        gene_names = [
            # let-7 pathway targets
            'lin-41', 'hbl-1', 'daf-12', 'lin-28', 'lin-14',
            # Developmental timing genes
            'cog-1', 'unc-3', 'daf-16', 'let-60', 'lin-15',
            # Metabolic genes (background)
            'aco-1', 'aco-2', 'idh-1', 'mdh-1', 'sdha-1',
            # Stress response genes
            'hsp-16.2', 'hsp-70', 'gst-4', 'sod-1', 'ctl-1',
            # Additional developmental genes
            'mab-5', 'egl-1', 'ced-3', 'ced-4', 'ced-9'
        ]
        
        # Create realistic expression patterns
        np.random.seed(42)  # For reproducible synthetic data
        n_genes = len(gene_names)
        n_samples = len(samples)
        
        # Base expression levels
        gene_data = np.random.negative_binomial(100, 0.3, size=(n_genes, n_samples))
        
        # Add stage-specific patterns for known interactions
        for i, gene in enumerate(gene_names):
            if gene in ['lin-41', 'hbl-1']:  # let-7 targets - decrease in late stages
                for j, sample in enumerate(samples):
                    if 'L3' in sample or 'L4' in sample:
                        gene_data[i, j] = int(gene_data[i, j] * 0.3)  # Strong downregulation
            elif gene in ['lin-14', 'lin-28']:  # lin-4 targets - decrease in early stages
                for j, sample in enumerate(samples):
                    if 'L2' in sample:
                        gene_data[i, j] = int(gene_data[i, j] * 0.4)  # Moderate downregulation
        
        gene_df = pd.DataFrame(gene_data, index=gene_names, columns=samples)
        
        # lncRNA expression data
        lncrna_names = [f'lincRNA-{i}' for i in range(1, 21)]  # 20 lncRNAs
        lncrna_data = np.random.negative_binomial(50, 0.4, size=(len(lncrna_names), n_samples))
        lncrna_df = pd.DataFrame(lncrna_data, index=lncrna_names, columns=samples)
        
        # miRNA expression data
        mirna_names = ['let-7', 'lin-4', 'mir-1', 'mir-48', 'mir-84', 'mir-241']
        mirna_data = np.random.negative_binomial(20, 0.5, size=(len(mirna_names), n_samples))
        
        # Add realistic patterns for known miRNAs
        for i, mirna in enumerate(mirna_names):
            if mirna == 'let-7':  # Increases during development
                for j, sample in enumerate(samples):
                    if 'L3' in sample or 'L4' in sample:
                        mirna_data[i, j] = int(mirna_data[i, j] * 3)  # Strong upregulation
            elif mirna == 'lin-4':  # Early expression
                for j, sample in enumerate(samples):
                    if 'L1' in sample or 'L2' in sample:
                        mirna_data[i, j] = int(mirna_data[i, j] * 2)  # Moderate upregulation
        
        mirna_df = pd.DataFrame(mirna_data, index=mirna_names, columns=samples)
        
        # DNA methylation data (CpG sites)
        methylation_sites = [f'chr{chr}_CpG_{site}' 
                           for chr in range(1, 6) 
                           for site in range(1, 11)]  # 50 CpG sites
        
        # Methylation percentages (0-100%)
        meth_data = np.random.beta(2, 3, size=(len(methylation_sites), n_samples)) * 100
        
        # Add patterns related to gene regulation
        # Simulate hypermethylation leading to gene silencing
        for i, site in enumerate(methylation_sites):
            if 'chr1' in site:  # Sites near let-7 targets
                for j, sample in enumerate(samples):
                    if 'L3' in sample or 'L4' in sample:
                        meth_data[i, j] = meth_data[i, j] * 0.5  # Hypomethylation for expression
        
        methylation_df = pd.DataFrame(meth_data, index=methylation_sites, columns=samples)
        
        logger.info(f"Created synthetic dataset with {n_samples} samples:")
        logger.info(f"  - {len(gene_names)} genes")
        logger.info(f"  - {len(lncrna_names)} lncRNAs") 
        logger.info(f"  - {len(mirna_names)} miRNAs")
        logger.info(f"  - {len(methylation_sites)} methylation sites")
        
        return {
            'gene_counts': gene_df,
            'lncrna_counts': lncrna_df, 
            'mirna_counts': mirna_df,
            'wgbs_counts': methylation_df
        }
    
    def save_datasets(self, datasets: Dict[str, pd.DataFrame]) -> None:
        """Save datasets in ConTra format."""
        logger.info("Saving datasets in ConTra format...")
        
        for data_type, df in datasets.items():
            filename = f"{data_type}_cleaned.csv"
            filepath = self.output_dir / filename
            df.to_csv(filepath)
            logger.info(f"Saved {filename}: {df.shape}")
            
            # Create summary file
            summary_file = self.output_dir / f"{data_type}_summary.txt"
            with open(summary_file, 'w') as f:
                f.write(f"Dataset: {data_type}\n")
                f.write(f"Features: {df.shape[0]}\n")
                f.write(f"Samples: {df.shape[1]}\n")
                f.write(f"Sample IDs: {', '.join(df.columns)}\n")
                f.write(f"Total non-zero values: {(df > 0).sum().sum()}\n")
                f.write(f"Mean expression per feature: {df.mean(axis=1).mean():.2f}\n")
                f.write(f"Std expression per feature: {df.mean(axis=1).std():.2f}\n")
    
    def create_validation_metadata(self) -> None:
        """Create metadata file with known interactions for validation."""
        metadata_file = self.output_dir / "validation_metadata.json"
        
        import json
        metadata = {
            'organism': 'Caenorhabditis elegans',
            'data_source': 'Synthetic modENCODE-like dataset',
            'developmental_stages': list(self.sample_mapping.keys()),
            'known_interactions': self.known_interactions,
            'validation_criteria': {
                'sensitivity_threshold': 0.8,
                'specificity_threshold': 0.8,
                'statistical_significance': 0.05
            },
            'expected_results': {
                'let-7_targets_downregulated': ['lin-41', 'hbl-1'],
                'lin-4_targets_downregulated': ['lin-14', 'lin-28'],
                'stage_specific_interactions': {
                    'L1-L2': ['lin-4 -> lin-14', 'lin-4 -> lin-28'],
                    'L3-L4': ['let-7 -> lin-41', 'let-7 -> hbl-1']
                }
            }
        }
        
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Created validation metadata: {metadata_file}")
    
    def create_readme(self) -> None:
        """Create README for the test dataset."""
        readme_content = """# C. elegans Test Dataset for ConTra Validation

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
"""
        
        readme_file = self.output_dir / "README.md"
        with open(readme_file, 'w') as f:
            f.write(readme_content)
        
        logger.info(f"Created README: {readme_file}")

def main():
    """Main function to download and format test dataset."""
    parser = argparse.ArgumentParser(description='Download C. elegans test dataset for ConTra validation')
    parser.add_argument('--output-dir', default='data/test_datasets/celegans_modENCODE',
                       help='Output directory for test dataset')
    
    args = parser.parse_args()
    
    # Create downloader
    downloader = CelegansTestDataDownloader(args.output_dir)
    
    # Create synthetic dataset
    datasets = downloader.create_synthetic_data()
    
    # Save in ConTra format
    downloader.save_datasets(datasets)
    
    # Create validation metadata
    downloader.create_validation_metadata()
    
    # Create README
    downloader.create_readme()
    
    logger.info(f"✅ Test dataset created successfully in {args.output_dir}")
    logger.info("To use with ConTra:")
    logger.info(f"  1. Copy files to data/cleaned_datasets/")
    logger.info(f"  2. Run ConTra analysis")
    logger.info(f"  3. Validate results against known interactions")

if __name__ == '__main__':
    main()