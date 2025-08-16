# ConTra: Context-Dependent Regulation Analysis

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Contributions
Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](CONTRIBUTING.md)
[![Issues](https://img.shields.io/badge/issues-open-orange.svg)](https://github.com/sr320/ConTra/issues)
[![PRs](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/sr320/ConTra/pulls)

**ConTra** is a high-performance Python framework for identifying
context-dependent regulatory interactions in multi-omics data. It
leverages parallel processing, vectorized operations, and
memory-efficient algorithms to analyze complex biological regulatory
networks.

## 🚀 Features

-   **Multi-omics Integration**: Analyzes gene expression, lncRNA,
    miRNA, and DNA methylation data
-   **High-Performance Computing**: Parallel processing across 48+ CPU
    cores with optimized memory usage
-   **Context-Dependent Analysis**: Identifies regulatory interactions
    that vary across different biological contexts
-   **Advanced Statistical Methods**:
    -   Interaction term analysis
    -   Conditional correlation analysis
    -   Multi-variable regression with interaction terms
    -   Context-specific regulatory network inference
-   **Memory Optimization**: Efficient batch processing using up to
    247GB RAM
-   **Comprehensive Output**: Generates plots, tables, and detailed
    reports

## 📋 Requirements

-   Python 3.8+
-   8GB+ RAM (recommended: 16GB+)
-   Multi-core CPU (recommended: 8+ cores)

## 🛠️ Installation

1.  Clone the repository:

``` bash
git clone https://github.com/sr320/ConTra.git
cd ConTra
```

2.  Install dependencies:

``` bash
pip install -r code/requirements.txt
```

or 

``` bash
python3 -m pip install -r code/requirements.txt
```


## 📊 Usage

### Data Format

The repository provides pre-cleaned and standardized multi-omics
datasets ready for immediate analysis. All datasets are in CSV format
with consistent sample alignment:

#### **Available Datasets**

| Dataset                     | Features | Samples | Sparsity | Description                           |
|---------------|---------------|---------------|---------------|---------------|
| `gene_counts_cleaned.csv`   | 36,084   | 40      | 37.8%    | Gene expression counts                |
| `lncrna_counts_cleaned.csv` | 15,900   | 40      | 3.8%     | Long non-coding RNA expression counts |
| `mirna_counts_cleaned.csv`  | 51       | 40      | 7.8%     | MicroRNA expression counts            |
| `wgbs_counts_cleaned.csv`   | 249      | 40      | 41.4%    | WGBS CpG methylation counts           |

#### **Sample Structure**

All datasets contain the same **40 samples** representing different time
points (TP1-TP4) across **10 different conditions**:

-   **ACR-139**: TP1, TP2, TP3, TP4
-   **ACR-145**: TP1, TP2, TP3, TP4\
-   **ACR-150**: TP1, TP2, TP3, TP4
-   **ACR-173**: TP1, TP2, TP3, TP4
-   **ACR-186**: TP1, TP2, TP3, TP4
-   **ACR-225**: TP1, TP2, TP3, TP4
-   **ACR-229**: TP1, TP2, TP3, TP4
-   **ACR-237**: TP1, TP2, TP3, TP4
-   **ACR-244**: TP1, TP2, TP3, TP4
-   **ACR-265**: TP1, TP2, TP3, TP4

#### **Data Quality Features**

-   **Common sample IDs**: All datasets use identical 40 sample
    identifiers
-   **No zero expression**: Features with zero expression across all
    samples removed
-   **Sufficient variation**: Features with limited variation (CV \<
    0.1) filtered out
-   **Consistent structure**: Same column order and sample alignment
    across all datasets
-   **No missing values**: Complete data matrices ready for analysis


#### **File Organization**

```         
data/cleaned_datasets/
├── gene_counts_cleaned.csv      # Main gene expression dataset
├── lncrna_counts_cleaned.csv    # Main lncRNA expression dataset
├── mirna_counts_cleaned.csv     # Main miRNA expression dataset
├── wgbs_counts_cleaned.csv      # Main DNA methylation dataset
├── *_summary.txt                # Individual dataset statistics
├── combined_summary.txt         # Overall dataset summary
└── README.md                    # Detailed data documentation
```

#### **Ready for Analysis**

These datasets are immediately usable for: - Multi-omics correlation
analysis - Time series analysis across TP1-TP4 time points -
Context-dependent regulatory network inference - Statistical modeling
and machine learning workflows

## 🔬 Methodology

ConTra employs several sophisticated approaches to identify
context-dependent regulatory interactions:

1.  **Interaction Term Analysis**: Examines how regulatory relationships
    change across different biological contexts
2.  **Conditional Correlation Analysis**: Identifies correlations that
    are context-specific
3.  **Multi-variable Regression**: Models complex regulatory networks
    with interaction terms
4.  **Network Inference**: Constructs context-specific regulatory
    networks

## 📁 Project Structure

```         
ConTra/
├── code/
│   ├── context_dependent_analysis.py      # Main analysis pipeline
│   ├── subset_context_dependent_analysis.py # Subset analysis tools
│   └── requirements.txt                   # Python dependencies
├── data/
│   └── cleaned_datasets/                 # Input data files
├── output/                               # Generated results (created at runtime)
├── LICENSE                               # MIT License
└── README.md                             # This file
```

### Script Differences

**`context_dependent_analysis.py`** - **Full Analysis Pipeline**
- Analyzes **ALL 36,084 genes** in the dataset
- Comprehensive regulatory interaction analysis across all genes
- Higher computational requirements but complete coverage

**`subset_context_dependent_analysis.py`** - **Subset Analysis Tools**
- Analyzes **500 genes** (randomly sampled from the full dataset)
- Faster execution for testing and development
- Lower computational requirements (~5 min on 48 cores)

## 🤝 Contributing

We welcome contributions from the community! Whether you're a
bioinformatician, data scientist, or developer, there are many ways to
contribute:

### How to Contribute

1.  **Fork** the repository
2.  **Create** a feature branch
    (`git checkout -b feature/amazing-feature`)
3.  **Commit** your changes (`git commit -m 'Add amazing feature'`)
4.  **Push** to the branch (`git push origin feature/amazing-feature`)
5.  **Open** a Pull Request



## 📝 License

This project is licensed under the MIT License - see the
[LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

-   **Steven Roberts** - Project maintainer and primary developer
-   **Open Source Community** - For the excellent libraries that make
    this project possible
-   **Contributors** - Everyone who has helped improve ConTra

## 📞 Contact

-   **Issues**: [GitHub
    Issues](https://github.com/sr320/ConTra/issues)
-   **Discussions**: [GitHub
    Discussions](https://github.com/sr320/ConTra/discussions)
-   


------------------------------------------------------------------------

**⭐ Star this repository if you find it useful!**

**🤝 Contributions are always welcome and appreciated!**
