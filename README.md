## ConTra: Context-Dependent Regulation Analysis

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/downloads/) [![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) [![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](CONTRIBUTING.md)

**ConTra** is a high-performance Python pipeline for identifying context-dependent regulatory interactions in multi-omics data. It leverages parallel processing, vectorized operations, and memory-efficient algorithms to analyze complex regulatory networks involving gene expression, lncRNA, miRNA, and DNA methylation.

## 🚀 Features

- **Multi-omics integration**: Jointly analyzes gene, lncRNA, miRNA, and DNA methylation data.
- **Context-dependent analysis**: Detects regulatory relationships that change across biological contexts (e.g. high/low miRNA or methylation states).
- **Advanced modeling**:
  - Interaction term analysis via linear regression.
  - Conditional correlations and context strength metrics.
  - Optional multi-regulator (multi-way) models.
- **High-performance computing**: Parallelized across many cores with attention to memory usage.
- **Rich outputs**: Tables of interactions, context-specific correlation networks, and human-readable reports (Markdown + HTML).

## 📋 Requirements

- **Python**: 3.8+
- **RAM**: 8GB+ (16GB+ recommended for full runs)
- **CPU**: Multi-core (8+ cores recommended)

Install dependencies from the repo root:

```bash
pip install -r code/requirements.txt
# or
python3 -m pip install -r code/requirements.txt
```

## 📊 Running `context_dependent_analysis.py`

### Overview

The core entry point for ConTra is **`code/context_dependent_analysis.py`**, which:

- **Loads four matched data matrices**: gene expression, lncRNA expression, miRNA expression, and DNA methylation.
- **Analyzes all genes** using regression models with interaction terms and conditional correlations.
- **Focuses by default on methylation–miRNA context**, the component that shows robust separation between real and randomized data.
- Optionally runs additional, exploratory modules (lncRNA–miRNA context, multi-way interactions) and context-specific correlation networks.

Each run creates a timestamped output directory:

- `output/context_dependent_analysis_YYYYMMDD_HHMMSS/`
  - `tables/` – CSV result tables (core: `methylation_mirna_context.csv`)
  - `plots/` – PNG figures
  - `reports/` – Markdown and HTML reports summarizing the analysis

### Quick start on cleaned data

If you have a cleaned dataset under `data/cleaned_datasets` (or another cleaned folder with the required files), run:

```bash
python code/context_dependent_analysis.py --data-dir data/cleaned_datasets
```

You can also pass an **absolute path** to `--data-dir`.

### Using multi-species cleaned data

If your repository includes species-specific cleaned folders such as:

- `data/full-species-24/cleaned_apul`
- `data/full-species-24/cleaned_peve`
- `data/full-species-24/cleaned_ptua`

you can run:

```bash
python code/context_dependent_analysis.py --data-dir data/full-species-24/cleaned_apul
python code/context_dependent_analysis.py --data-dir data/full-species-24/cleaned_peve
```

Each invocation produces its own timestamped output folder under `output/`.

### Running on your own data

To use your own multi-omics data, prepare a directory that contains **four CSV files** with **identical sample columns**:

- **Required file names** (must match exactly):
  - `gene_counts_cleaned.csv`
  - `lncrna_counts_cleaned.csv`
  - `mirna_counts_cleaned.csv`
  - `wgbs_counts_cleaned.csv`

- **File layout**:
  - **Rows**: features (genes, lncRNAs, miRNAs, CpGs).
  - **Columns**: samples, with the **same IDs and order** in all four files.
  - First column: feature identifier (used as the index).

- **Recommended preprocessing**:
  - Remove features with zero counts across all samples.
  - Filter out extremely low-variance features.
  - Ensure **no missing values**.

Then point the analysis at your folder:

```bash
python code/context_dependent_analysis.py --data-dir /absolute/path/to/your_cleaned_folder
```

If sample IDs are not aligned or a required file is missing, the script will raise an informative error.

### Key CLI options

- **Core options**:
  - **`--data-dir PATH`**: directory containing the four cleaned CSVs (default: `data/cleaned_datasets`).
  - **`--n-jobs INT`**: number of parallel worker processes (default: up to 48, or all available cores).

- **Empirical FDR (recommended for robustness)**:
  - **`--enable-empirical-fdr`**: generate randomized null datasets and compute empirical FDR thresholds for methylation–miRNA `context_strength`.
  - **`--null-replicates INT`** (default: `1`): number of randomized null datasets to analyze.
  - **`--fdr-alpha FLOAT`** (default: `0.3`): target empirical FDR level (e.g. `0.3` for exploratory, `0.1` for stricter).
  - **`--random-seed INT`**: base random seed for reproducible null generation.

- **Optional exploratory modules**:
  - **`--enable-lncrna-context`**:
    - Runs the lncRNA–miRNA–gene context module.
    - Current evidence suggests its context metrics behave similarly in real and randomized data; treat results as exploratory.
  - **`--enable-multi-way`**:
    - Runs the multi-way interaction module (many regulators per gene).
    - Large improvements can also appear in null data; use for hypothesis generation.

Example usage with empirical FDR and exploratory modules:

```bash
python code/context_dependent_analysis.py \
  --data-dir data/full-species-24/cleaned_apul \
  --enable-empirical-fdr \
  --null-replicates 1 \
  --fdr-alpha 0.3 \
  --enable-lncrna-context \
  --enable-multi-way
```

### Main outputs

- **Core interaction table**: `methylation_mirna_context.csv`
  - One row per gene–CpG–miRNA triplet.
  - Contains:
    - Regression metrics: `r2_regulator1_only`, `r2_regulator1_regulator2`, `r2_with_interaction`.
    - Effect sizes: `improvement_from_regulator2`, `improvement_from_interaction`.
    - Conditional correlations: `corr_high_regulator2`, `corr_low_regulator2`.
    - **`context_strength`** quantifying change in correlation across contexts.
    - When empirical FDR is enabled:
      - `empirical_fdr_threshold`
      - `empirical_fdr_estimated`
      - `empirical_fdr_significant`
  - **Primary table for downstream biological interpretation.**

- **Context correlation networks** (exploratory):
  - `high_methylation_gene_methylation_correlations.csv`
  - `high_methylation_gene_mirna_correlations.csv`
  - `high_methylation_gene_lncrna_correlations.csv`
  - `low_mirna_gene_methylation_correlations.csv`
  - `low_mirna_gene_mirna_correlations.csv`
  - `low_mirna_gene_lncrna_correlations.csv`
  - Each row is a gene–regulator pair with Pearson `correlation` and raw `p_value` in a specific context (high methylation or low miRNA).
  - Best used for **network-style exploration**; p-values are uncorrected and should not be treated as strict significance.

- **Reports**:
  - `reports/context_dependent_analysis_report.md`
  - `reports/context_dependent_analysis_report.html`
  - Summarize:
    - Dataset overview and resources used.
    - Key methylation–miRNA context metrics (including empirical FDR if enabled).
    - Status and top results for exploratory modules.
    - A **“Table Definitions and Statistical Confidence”** section describing each major table and how to interpret its statistics.

## 🔬 Statistical robustness of `context_dependent_analysis.py`

### Empirically strong core: methylation–miRNA context

- **What is robust**:
  - Across multiple real vs randomized comparisons, **methylation–miRNA context** shows:
    - Higher baseline R² in real data (`r2_regulator1_only`, `r2_regulator1_regulator2`).
    - Stronger gene–regulator correlations.
    - Larger **`context_strength`** in real data than in randomized null runs.
- **Empirical FDR**:
  - The script can generate randomized null datasets by **shuffling samples within each feature**, preserving marginal distributions while destroying gene–regulator structure.
  - It then estimates a context_strength threshold such that
    \[
    \widehat{\text{FDR}} \approx \frac{\text{null hits above threshold}}{\text{real hits above threshold}} \le q_{\text{target}}
    \]
  - Interactions with `empirical_fdr_significant == True` are those that pass this **data-driven FDR control** and are the most reliable candidates.
- **Practical guidance**:
  - For production analyses, enable:
    - `--enable-empirical-fdr`
    - A reasonable `--fdr-alpha` (e.g. 0.2–0.3 for exploratory, 0.1 for stricter).
  - Focus interpretation on **FDR-significant rows** and use `context_strength` as the primary effect-size indicator.

### Exploratory components: lncRNA context and multi-way interactions

- **lncRNA–miRNA context**:
  - Baseline lncRNA→gene relationships are strong in real data.
  - However, current context interaction metrics (additional improvements, `context_strength`) often look similar or even stronger in randomized data.
  - Use this module to:
    - Generate hypotheses about lncRNAs associated with genes that already have strong methylation–miRNA context.
    - Explore ceRNA-like patterns that warrant external validation.

- **Multi-way interactions**:
  - Multi-regulator models can show large improvements in both real and randomized datasets.
  - Few genes pass stringent significance thresholds, and null runs can exhibit similar improvement magnitudes.
  - Treat these results as **hypothesis-generating only**, and always cross-check with simpler, better-behaved context metrics (especially methylation–miRNA context with empirical FDR).

### Interpreting p-values vs empirical FDR

- **F-test p-values** in the tables:
  - Used internally for flags like `context_dependent` and `has_significant_interactions`.
  - Useful but can be misleading when considered alone, because some apparent improvements also occur in randomized data.

- **Empirical FDR** (recommended where available):
  - Explicitly measures how often a metric as large as observed in real data appears in randomized null runs.
  - Provides a more realistic sense of **false discovery rate** under the actual analysis pipeline and dataset structure.

- **Correlation-network p-values**:
  - In the high-/low-context correlation tables, p-values are uncorrected across many tests.
  - Use them to **rank** gene–regulator pairs; rely on:
    - Strong |correlation| values,
    - Consistency across contexts,
    - Overlap with FDR-supported interactions from `methylation_mirna_context.csv`.

In summary, **methylation–miRNA context with empirical FDR** is the statistically most robust output of `context_dependent_analysis.py`, while other modules and correlation networks are preserved as clearly labeled exploratory tools for deeper biological investigation.

## 📝 License

This project is licensed under the MIT License – see the `LICENSE` file for details.

## 🙏 Acknowledgments

- **Steven Roberts** – project maintainer and primary developer  
- **Open Source Community** – for the libraries that make this project possible  
- **Contributors** – everyone who has helped improve ConTra

## 📞 Contact

- **Issues**: GitHub Issues  
- **Discussions**: GitHub Discussions

