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

- **Multi-omics Integration**: Gene, lncRNA, miRNA, and DNA methylation
- **Dual Analysis Modes**: Single unified script supports `full` and `subset` modes
- **High-Performance Computing**: Parallel processing across all detected CPU cores (override with `--n-jobs` or `CONTRA_MAX_CORES`)
- **Context-Dependent Analysis**: Interaction + conditional + multi-way modeling
- **Advanced Statistical Methods**:
    - Interaction term regression & F‑tests
    - Vectorized correlation screening
    - Multi-way regulator synergy scoring
    - Context-specific network inference
- **Robust Output Suite**: Publication-ready plots, CSV tables, Markdown + HTML reports
- **Reproducible Subset Mode**: Fixed random seed for quick iteration
- **Memory Efficiency**: Chunked / vectorized operations

## 📋 Requirements

- Python 3.8+
- 8GB+ RAM (recommended: 16GB+)
- Multi-core CPU (recommended: 8+ cores)

## 🛠️ Installation

1. Clone the repository:

```bash
git clone https://github.com/sr320/ConTra.git
cd ConTra
```

2. Install dependencies:

```bash
pip install -r code/requirements.txt
```

or

```bash
python3 -m pip install -r code/requirements.txt
```


## 📊 Usage

### 1. Quick Start

Run (interactive prompt will ask for mode, default = full):

```bash
python3 code/context_dependent_analysis.py
```

Run explicitly in subset (faster dev/test) mode using 8 workers:

```bash
python3 code/context_dependent_analysis.py --mode subset --n-jobs 8
```

Run full analysis using all detected cores:

```bash
python3 code/context_dependent_analysis.py --mode full
```

Arguments:

- `--mode {full,subset}` Select analysis breadth
- `--n-jobs N` Override auto CPU core detection

#### Parallelism / CPU Control

By default the pipeline uses all available CPU cores reported by Python.

Ways to control cores:

| Method | Example | Notes |
|--------|---------|-------|
| Command-line flag | `--n-jobs 32` | Explicitly sets worker count |
| Environment variable | `CONTRA_MAX_CORES=64` | Upper cap when `--n-jobs` not provided |
| Both provided | `CONTRA_MAX_CORES=64 --n-jobs 80` | Flag wins (uses 80 if system has ≥80 cores) |

Examples:

```bash
# Cap to 64 cores via environment variable
CONTRA_MAX_CORES=64 python3 code/context_dependent_analysis.py --mode full

# Explicitly use 32 cores (ignores CONTRA_MAX_CORES)
python3 code/context_dependent_analysis.py --mode subset --n-jobs 32
```

#### Optional Adaptive Permutation Testing

You can obtain empirical p-values for top correlations using on-demand, early-stopping permutations. Disabled by default to keep the default run fast.

Flags:

- `--perm` Enable adaptive permutation testing
- `--perm-min INT` Minimum permutations per tested correlation (default: 1000)
- `--perm-max INT` Maximum permutations (default: 100000)
- `--perm-alpha FLOAT` Target tail probability precision (default: 0.001). Stops early when achievable precision reached or max permutations hit.

How it works (per correlation sign):

1. Generate permutations in batches (inverting the correlation by shuffling one vector)
2. Track exceedances (|r_perm| >= |r_obs|)
3. After each batch, form empirical p = (exceed + 1) / (n_perm + 1)
4. Compute a two-sided binomial CI width; stop if max(n_perm) reached or CI half-width < `perm_alpha/2`

Interpretation:

- Reported p is conservative (add-one smoothing) and bounded below by 1/(n_perm+1)
- Very strong correlations may stop early (e.g. after ~2–5k permutations) saving time
- Weak correlations accumulate permutations until they can be confidently declared non-significant

Example (subset mode, 16 workers, enable permutations):

```bash
python3 code/context_dependent_analysis.py --mode subset --n-jobs 16 --perm --perm-min 2000 --perm-max 50000 --perm-alpha 0.0005
```

Result columns (where applicable) gain `empirical_p` alongside existing statistical metrics.

Performance tips:

- Start with smaller `--perm-min` (e.g. 1000) to gauge runtime
- Increase `--perm-max` only if you need finer p-value resolution (<1e-4)
- Tightening `--perm-alpha` increases runtime; loosening speeds it up
- Use subset mode first to benchmark permutation cost before running full

Outputs (per run) are written to:

```text
output/context_dependent_analysis_<mode>_<YYYYMMDD_HHMMSS>/
    plots/   *.png
    tables/  *.csv
    reports/ *.md, *.html
```

Key tables:

- `methylation_mirna_context.csv`
- `lncrna_mirna_context.csv`
- `multi_way_interactions.csv`
- `*gene_*_correlations.csv` (context-specific networks)

Reports:

- Markdown: `context_dependent_analysis_report.md` (full) or `subset_context_dependent_analysis_report.md`
- HTML copy with embedded images

### 2. Analysis Modes

| Mode | Genes (pairwise) | Genes (multi-way) | Genes (networks) | miRNA top/use | Methylation top/use | lncRNA top/use | Multi-way regulators (miRNA / lncRNA / methylation) | Seed |
|------|------------------|-------------------|------------------|---------------|---------------------|----------------|------------------------------------------------------|------|
| full   | all | all | all | 25 / 10 | 50 / 15 | 50 / 15 | 15 / 30 / 25 | none |
| subset | 500 | 200 | 200 | 10 / 5  | 10 / 5  | 10 / 5  | 5 / 7 / 5    | 42   |

Subset mode greatly reduces runtime and file sizes while preserving pipeline logic (useful for method development / CI tests).

### 3. Data Format

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

#### **Data Quality Features**

- **Common sample IDs**: All datasets use identical 40 sample identifiers
- **No zero expression**: Features with zero expression across all samples removed
- **Sufficient variation**: Features with limited variation (CV < 0.1) filtered out
- **Consistent structure**: Same column order and sample alignment across all datasets
- **No missing values**: Complete data matrices ready for analysis


#### **File Organization**

```text
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

1. **Interaction Term Analysis**: Examines how regulatory relationships change across different biological contexts
2. **Conditional Correlation Analysis**: Identifies correlations that are context-specific
3. **Multi-variable Regression**: Models complex regulatory networks with interaction terms
4. **Network Inference**: Constructs context-specific regulatory networks

## 📁 Project Structure

```text
ConTra/
├── code/
│   ├── context_dependent_analysis.py  # Unified full + subset analysis script
│   └── requirements.txt               # Python dependencies
├── data/
│   └── cleaned_datasets/             # Input data files
├── output/                           # Generated results (created at runtime)
├── LICENSE                           # MIT License
└── README.md                         # This file
```

Deprecated: The previous `subset_context_dependent_analysis.py` has been merged—use the unified script with `--mode subset`.

## 🤝 Contributing

We welcome contributions from the community! Whether you're a
bioinformatician, data scientist, or developer, there are many ways to
contribute:

### How to Contribute

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5. **Open** a Pull Request



## 📝 License

This project is licensed under the MIT License - see the
[LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Steven Roberts** - Project maintainer and primary developer
- **Open Source Community** - For the excellent libraries that make this project possible
- **Contributors** - Everyone who has helped improve ConTra

## 📞 Contact

- **Issues**: [GitHub Issues](https://github.com/sr320/ConTra/issues)
- **Discussions**: [GitHub Discussions](https://github.com/sr320/ConTra/discussions)


------------------------------------------------------------------------

**⭐ Star this repository if you find it useful!**

**🤝 Contributions are always welcome and appreciated!**
