#!/usr/bin/env python3
"""
Unified Context-Dependent Regulation Analysis (Full/Subet Modes)

This script merges previous full and subset analysis scripts. Features:
- Mode selection (full / subset) via CLI or interactive prompt
- Mode-specific sampling and regulator breadth
- Robust plotting with error handling (from subset script)
- Parallelized and vectorized computations

Configuration summary:
Full mode:
    • All genes for pairwise, multi-way, and context network analyses
    • Regulators: miRNA top 25 (use 10), methylation top 50 (use 15), lncRNA top 50 (use 15)
    • Multi-way regulators: miRNA15 / lncRNA30 / methylation25
Subset mode:
    • Pairwise genes 500, multi-way 200, context networks 200 (random, seed=42)
    • Regulators: miRNA top 10 (use 5), methylation top 10 (use 5), lncRNA top 10 (use 5)
    • Multi-way regulators: miRNA5 / lncRNA7 / methylation5
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import warnings
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import os
import time
import base64
import gc
import psutil
from datetime import datetime
from typing import Dict, List, Tuple, Any

# Attempt to select a more efficient start method on Unix-like systems to reduce
# process spawn overhead and improve CPU saturation. 'fork' avoids re-pickling
# large objects when safe. If already set or unsupported, this is silently ignored.
try:  # noqa: E722
    if os.name == 'posix':
        mp.set_start_method('fork', force=False)
except Exception:  # noqa: E722
    pass

warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class OptimizedContextDependentRegulationAnalysis:
    """Unified context-dependent regulation analysis supporting full and subset modes."""

    CONFIG = {
        'full': {
            'pairwise_gene_limit': None,
            'multi_way_gene_limit': None,
            'methylation_mirna': {'miRNA_top': 25, 'miRNA_use': 10, 'meth_top': 50, 'meth_use': 15},
            'multi_way_regulators': {'miRNA': 15, 'lncRNA': 30, 'methylation': 25},
            'seed': None
        },
        'subset': {
            'pairwise_gene_limit': 500,
            'multi_way_gene_limit': 200,
            'methylation_mirna': {'miRNA_top': 10, 'miRNA_use': 5, 'meth_top': 10, 'meth_use': 5},
            'multi_way_regulators': {'miRNA': 5, 'lncRNA': 7, 'methylation': 5},
            'seed': 42
        }
    }

    def __init__(self, data_dir: str = "data/cleaned_datasets", n_jobs: int = None, mode: str = 'full', output_root: str = 'output'):
        """Initialize analysis object.

        Parameters
        ----------
        data_dir : str
            Relative path to cleaned datasets directory.
        n_jobs : int | None
            Number of parallel workers (default: all cores, optionally capped by CONTRA_MAX_CORES env).
        mode : {'full','subset'}
            Analysis mode controlling sampling/regulator sizes.
        output_root : str
            Base output directory (relative to repo root).
        """
        mode = (mode or 'full').lower()
        if mode not in self.CONFIG:
            raise ValueError("mode must be 'full' or 'subset'")
        self.mode = mode

        # Paths & data containers
        self.workspace_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.data_dir = os.path.join(self.workspace_root, data_dir)
        self.datasets: Dict[str, pd.DataFrame] = {}
        self.results: Dict[str, Any] = {}
        
        # Add caching for correlation calculations
        self.correlation_cache: Dict[str, Tuple[float, float]] = {}

        # Parallel workers (allow env cap)
        if n_jobs is not None:
            self.n_jobs = n_jobs
        else:
            max_env = os.environ.get('CONTRA_MAX_CORES')
            try:
                max_env_val = int(max_env) if max_env else None
            except ValueError:
                max_env_val = None
            cpu_total = mp.cpu_count()
            self.n_jobs = min(cpu_total, max_env_val) if max_env_val else cpu_total

        # Seed
        self.seed = self.CONFIG[self.mode]['seed']
        if self.seed is not None:
            np.random.seed(self.seed)

        # Output directories
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.base_output_root = os.path.join(self.workspace_root, output_root.rstrip('/'))
        self.output_dir = os.path.join(self.base_output_root, f"context_dependent_analysis_{self.mode}_{self.timestamp}")
        os.makedirs(self.output_dir, exist_ok=True)
        self.plots_dir = os.path.join(self.output_dir, "plots")
        self.tables_dir = os.path.join(self.output_dir, "tables")
        self.reports_dir = os.path.join(self.output_dir, "reports")
        for d in (self.plots_dir, self.tables_dir, self.reports_dir):
            os.makedirs(d, exist_ok=True)

        print(f"🚀 Initializing {self.mode.upper()} context-dependent analysis with {self.n_jobs} parallel workers")
        print(f"💾 Available RAM: {self._get_available_ram():.1f} GB")
        print(f"📁 Output directory: {self.output_dir}")

        # Load data
        self.load_datasets()

        # Optional task granularity multiplier (env: CONTRA_CHUNK_MULT)
        try:
            self.chunk_multiplier = max(1, int(os.environ.get('CONTRA_CHUNK_MULT', '1')))
        except ValueError:
            self.chunk_multiplier = 1
        # No persistent process pool: create per stage to avoid pickling locks on macOS spawn method
        self._pool = None  # placeholder (not used directly)

        # Permutation configuration (can be overridden after construction or CLI or post-init)
        self.permutation_enabled = False
        self.permutation_min = 100
        self.permutation_max = 1000
        self.permutation_alpha = 0.05
        self.permutation_batch = 25

    def _get_available_ram(self):
        """Get available RAM in GB."""
        try:
            return psutil.virtual_memory().available / (1024**3)
        except:
            try:
                with open('/proc/meminfo', 'r') as f:
                    for line in f:
                        if line.startswith('MemAvailable:'):
                            return int(line.split()[1]) / (1024**2)
            except:
                pass
        return 8.0  # Conservative fallback
    
    def _monitor_memory_usage(self):
        """Monitor and log current memory usage."""
        try:
            memory = psutil.virtual_memory()
            process = psutil.Process()
            print(f"    Memory: {memory.percent:.1f}% used, Process: {process.memory_info().rss / 1024**2:.1f} MB")
            
            # Trigger garbage collection if memory usage is high
            if memory.percent > 80:
                print("    🧹 High memory usage detected, running garbage collection...")
                gc.collect()
                
        except Exception as e:
            print(f"    ⚠️  Memory monitoring error: {e}")
            
    def _optimize_memory_usage(self):
        """Optimize memory usage by clearing caches and running GC."""
        # Clear correlation cache if it gets too large
        if len(self.correlation_cache) > 10000:
            print("    🧹 Clearing correlation cache...")
            self.correlation_cache.clear()
        
        # Force garbage collection
        gc.collect()
        print(f"    🧹 Memory optimization complete")

    def load_datasets(self):
        """Load all cleaned datasets with parallel processing and memory optimization."""
        print("📂 Loading cleaned datasets with parallel processing...")
        
        # Load datasets in parallel using ThreadPoolExecutor for I/O operations
        with ThreadPoolExecutor(max_workers=4) as executor:
            future_gene = executor.submit(pd.read_csv, f"{self.data_dir}/gene_counts_cleaned.csv", index_col=0)
            future_lncrna = executor.submit(pd.read_csv, f"{self.data_dir}/lncrna_counts_cleaned.csv", index_col=0)
            future_mirna = executor.submit(pd.read_csv, f"{self.data_dir}/mirna_counts_cleaned.csv", index_col=0)
            future_methylation = executor.submit(pd.read_csv, f"{self.data_dir}/wgbs_counts_cleaned.csv", index_col=0)
            
            self.datasets['gene'] = future_gene.result()
            self.datasets['lncrna'] = future_lncrna.result()
            self.datasets['mirna'] = future_mirna.result()
            self.datasets['methylation'] = future_methylation.result()
        
        print(f"✅ Loaded gene expression: {self.datasets['gene'].shape}")
        print(f"✅ Loaded lncRNA expression: {self.datasets['lncrna'].shape}")
        print(f"✅ Loaded miRNA expression: {self.datasets['mirna'].shape}")
        print(f"✅ Loaded DNA methylation: {self.datasets['methylation'].shape}")
        
        # Verify sample alignment
        self.verify_sample_alignment()
        
        # Pre-compute data arrays for faster operations
        self._precompute_data_arrays()
        
    def _precompute_data_arrays(self):
        """Pre-compute numpy arrays for memory efficiency."""
        print("🔢 Pre-computing data arrays...")
        
        # Convert to numpy arrays for faster operations - use float32 to reduce memory
        self.gene_array = self.datasets['gene'].values.astype(np.float32)
        self.lncrna_array = self.datasets['lncrna'].values.astype(np.float32)
        self.mirna_array = self.datasets['mirna'].values.astype(np.float32)
        self.methylation_array = self.datasets['methylation'].values.astype(np.float32)
        
        # Create index mappings for faster lookups
        self.gene_index = {gene: idx for idx, gene in enumerate(self.datasets['gene'].index)}
        self.mirna_index = {mirna: idx for idx, mirna in enumerate(self.datasets['mirna'].index)}
        self.methylation_index = {meth: idx for idx, meth in enumerate(self.datasets['methylation'].index)}
        
        print("✅ Data arrays pre-computed with float32 precision and index mappings")
        
    def verify_sample_alignment(self):
        """Verify that all datasets have the same sample structure."""
        sample_sets = [set(df.columns) for df in self.datasets.values()]
        if not all(sample_sets[0] == s for s in sample_sets):
            raise ValueError("Sample IDs are not aligned across datasets!")
        
        self.samples = sorted(list(sample_sets[0]))
        self.n_samples = len(self.samples)
        print(f"✓ All datasets aligned with {self.n_samples} samples")
        
        # Extract time points and conditions
        self.time_points = sorted(list(set([s.split('-')[-1] for s in self.samples])))
        self.conditions = sorted(list(set([s.split('-')[1] for s in self.samples])))
        print(f"Time points: {self.time_points}")
        print(f"Conditions: {self.conditions}")

    def analyze_context_dependent_regulation(self):
        """Main analysis orchestrator."""
        print("\n" + "="*80)
        print(f"🚀 OPTIMIZED CONTEXT-DEPENDENT REGULATION ANALYSIS (MODE: {self.mode.upper()})")
        print("="*80)
        print(f"Using {self.n_jobs} parallel workers | Mode: {self.mode}")
        print("="*80)
        
        start_time = time.time()
        self._monitor_memory_usage()
        
        # 1. Analyze methylation-gene interactions dependent on miRNA levels (parallelized)
        print("\n1. 🔄 Analyzing methylation-gene interactions dependent on miRNA levels (parallelized)...")
        methylation_mirna_context = self.parallel_analyze_methylation_mirna_context()
        
        # Memory optimization after heavy computation
        self._optimize_memory_usage()
        
        # 2. Analyze multi-way regulatory interactions (parallelized)
        print("\n2. 🔄 Analyzing multi-way regulatory interactions (parallelized)...")
        multi_way_interactions = self.parallel_analyze_multi_way_interactions()
        
        # Final memory cleanup
        self._optimize_memory_usage()
        
        # Store results
        self.results['context_dependent'] = {
            'methylation_mirna_context': methylation_mirna_context,
            'multi_way_interactions': multi_way_interactions
        }
        
        total_time = time.time() - start_time
        print(f"\n⏱️  Total analysis time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
        print("✅ Context-dependent analysis completed with parallel processing!")
        self._monitor_memory_usage()
        
    def parallel_analyze_methylation_mirna_context(self):
        """Parallel methylation/miRNA context interaction analysis."""
        print("  🔄 Parallel processing methylation-miRNA context analysis...")
        cfg = self.CONFIG[self.mode]['methylation_mirna']
        
        # Gene sampling
        if self.CONFIG[self.mode]['pairwise_gene_limit'] is None:
            sampled_genes = self.datasets['gene'].index.to_list()
        else:
            limit = min(self.CONFIG[self.mode]['pairwise_gene_limit'], len(self.datasets['gene']))
            sampled_genes = np.random.choice(self.datasets['gene'].index, limit, replace=False)
        sampled_genes = list(sampled_genes)
        n_genes = len(sampled_genes)
        
        if n_genes == 0:
            return pd.DataFrame()
        
        # Adaptive chunk sizing based on memory and computation complexity
        memory_gb = self._get_available_ram()
        base_chunk_size = max(10, min(100, int(memory_gb / 2)))  # Conservative memory usage
        n_chunks = max(1, min(self.n_jobs * self.chunk_multiplier, n_genes // base_chunk_size + 1))
        gene_chunks = [chunk for chunk in np.array_split(sampled_genes, n_chunks) if len(chunk) > 0]
        
        print(f"  📊 Processing {n_genes} genes in {len(gene_chunks)} chunks (mode={self.mode}, chunk_size≈{n_genes//len(gene_chunks)})")
        
        all_results: List[Dict] = []
        with ProcessPoolExecutor(max_workers=self.n_jobs) as pool:
            futures = {pool.submit(self._process_methylation_mirna_chunk, list(chunk), cfg): i for i, chunk in enumerate(gene_chunks)}
            for fut in as_completed(futures):
                idx = futures[fut]
                try:
                    res = fut.result()
                    all_results.extend(res)
                    print(f"    ✅ Chunk {idx+1}/{len(gene_chunks)}: {len(res)} results")
                except Exception as e:
                    print(f"    ❌ Chunk {idx+1} failed: {e}")
        
        if not all_results:
            return pd.DataFrame()
        
        df = pd.DataFrame(all_results)
        print(f"  🎯 Total methylation-miRNA context interactions: {len(df)}")
        return df

    def _process_methylation_mirna_chunk(self, gene_chunk: List[str], cfg: Dict) -> List[Dict]:
        """Process a chunk of genes for methylation-miRNA context analysis."""
        results: List[Dict] = []
        
        # Pre-calculate correlations to avoid redundant computation
        chunk_start_time = time.time()
        
        for i, gene in enumerate(gene_chunk):
            try:
                # Use pre-computed arrays for faster access
                if gene in self.gene_index:
                    gene_idx = self.gene_index[gene]
                    gene_expression = self.gene_array[gene_idx]
                else:
                    gene_expression = self.datasets['gene'].loc[gene].values.astype(np.float32)
                
                # Get top regulators with caching
                cache_key_mirna = f"mirna_{gene}"
                cache_key_meth = f"methylation_{gene}"
                
                if cache_key_mirna not in self.correlation_cache:
                    mirna_corrs = self._get_top_correlations_vectorized(
                        gene_expression, self.datasets['mirna'], 'mirna', top_n=cfg['miRNA_top']
                    )
                    self.correlation_cache[cache_key_mirna] = mirna_corrs
                else:
                    mirna_corrs = self.correlation_cache[cache_key_mirna]
                
                if cache_key_meth not in self.correlation_cache:
                    meth_corrs = self._get_top_correlations_vectorized(
                        gene_expression, self.datasets['methylation'], 'methylation', top_n=cfg['meth_top']
                    )
                    self.correlation_cache[cache_key_meth] = meth_corrs
                else:
                    meth_corrs = self.correlation_cache[cache_key_meth]
                
                # Process interactions
                for mirna_name, _, _ in mirna_corrs[:cfg['miRNA_use']]:
                    for meth_name, _, _ in meth_corrs[:cfg['meth_use']]:
                        result = self._analyze_methylation_mirna_interaction(
                            gene, gene_expression, mirna_name, meth_name
                        )
                        if result:
                            results.append(result)
                            
                # Progress indicator for long-running chunks
                if i > 0 and i % 10 == 0:
                    elapsed = time.time() - chunk_start_time
                    rate = i / elapsed
                    eta = (len(gene_chunk) - i) / rate if rate > 0 else 0
                    print(f"      Progress: {i}/{len(gene_chunk)} genes ({rate:.1f} genes/sec, ETA: {eta:.0f}s)")
                    
            except Exception as e:
                print(f"    ⚠️  Error processing gene {gene}: {e}")
                continue
                
        return results

    def _analyze_methylation_mirna_interaction(self, gene_name: str, gene_expression: np.ndarray,
                                               mirna_name: str, meth_name: str) -> Dict:
        """Analyze methylation (regulator1) and miRNA (regulator2) interaction on a gene.

        Returns a result dict or None on failure.
        """
        try:
            # Use pre-computed arrays and indices for faster access
            meth_clean_name = meth_name.replace('methylation_', '')
            mirna_clean_name = mirna_name.replace('mirna_', '')
            
            if meth_clean_name in self.methylation_index:
                meth_idx = self.methylation_index[meth_clean_name]
                meth_expr = self.methylation_array[meth_idx]
            else:
                meth_expr = self.datasets['methylation'].loc[meth_clean_name].values.astype(np.float32)
                
            if mirna_clean_name in self.mirna_index:
                mirna_idx = self.mirna_index[mirna_clean_name]
                mirna_expr = self.mirna_array[mirna_idx]
            else:
                mirna_expr = self.datasets['mirna'].loc[mirna_clean_name].values.astype(np.float32)
            
            # Use numpy arrays directly for faster computation
            data_array = np.column_stack([gene_expression.astype(np.float32), meth_expr, mirna_expr])
            
            # Fast standardization
            data_scaled = (data_array - data_array.mean(axis=0)) / data_array.std(axis=0)
            
            # Add interaction term
            interaction = data_scaled[:, 1] * data_scaled[:, 2]
            
            # Sequential models using sklearn's efficient implementations
            target = data_scaled[:, 0]
            regulator1 = data_scaled[:, 1:2]
            regulator1_regulator2 = data_scaled[:, 1:3]
            full_features = np.column_stack([data_scaled[:, 1:3], interaction])
            
            # Fit models
            model1 = LinearRegression(fit_intercept=False).fit(regulator1, target)
            r2_1 = model1.score(regulator1, target)
            
            model2 = LinearRegression(fit_intercept=False).fit(regulator1_regulator2, target)
            r2_2 = model2.score(regulator1_regulator2, target)
            
            model3 = LinearRegression(fit_intercept=False).fit(full_features, target)
            r2_3 = model3.score(full_features, target)
            
            improvement_from_regulator2 = r2_2 - r2_1
            improvement_from_interaction = r2_3 - r2_2
            
            # Optimized F-test for interaction term
            n_samples = len(data_scaled)
            df1 = 1
            df2 = n_samples - 4  # parameters: intercept + 3 predictors
            
            if (df2 > 0 and improvement_from_interaction > 0 and r2_3 < 1.0 and (1 - r2_3) > 1e-10):
                try:
                    f_stat = (improvement_from_interaction / df1) / ((1 - r2_3) / df2)
                    p_value = 1 - stats.f.cdf(f_stat, df1, df2)
                    context_dependent = p_value < 0.05
                except Exception:
                    context_dependent = False
                    p_value = 1.0
            else:
                context_dependent = False
                p_value = 1.0
            
            # Faster conditional correlations using boolean indexing
            regulator2_scaled = data_scaled[:, 2]
            high_mask = regulator2_scaled > 0.5
            low_mask = regulator2_scaled < -0.5
            
            if np.sum(high_mask) > 2 and np.sum(low_mask) > 2:
                # Fast correlation using numpy
                target_high = target[high_mask]
                regulator1_high = data_scaled[high_mask, 1]
                target_low = target[low_mask]
                regulator1_low = data_scaled[low_mask, 1]
                
                corr_high = np.corrcoef(target_high, regulator1_high)[0, 1]
                corr_low = np.corrcoef(target_low, regulator1_low)[0, 1]
                
                context_strength = abs(corr_high - corr_low)
                context_direction = 'positive' if corr_high > corr_low else 'negative'
            else:
                corr_high = corr_low = context_strength = np.nan
                context_direction = 'NA'
            
            return {
                'interaction_type': 'methylation_mirna',
                'target': gene_name,
                'regulator1': meth_name,
                'regulator2': mirna_name,
                'r2_regulator1_only': r2_1,
                'r2_regulator1_regulator2': r2_2,
                'r2_with_interaction': r2_3,
                'improvement_from_regulator2': improvement_from_regulator2,
                'improvement_from_interaction': improvement_from_interaction,
                'context_dependent': context_dependent,
                'corr_high_regulator2': corr_high,
                'corr_low_regulator2': corr_low,
                'context_strength': context_strength,
                'context_direction': context_direction,
                'interaction_p_value': p_value
            }
        except Exception as e:  # pragma: no cover
            print(f"    ⚠️  Error analyzing {gene_name}-{mirna_name}-{meth_name}: {e}")
            return None

    def _get_top_correlations_vectorized(self, target_expression: np.ndarray, regulator_df: pd.DataFrame,
                                         prefix: str, top_n: int = 25) -> List[Tuple[str, float, float]]:
        """Return top correlated regulators with vectorized correlation computation.

        Returns list of (name_with_prefix, corr, p_value) sorted by |corr| desc.
        """
        # Vectorized correlation computation for all regulators at once
        target = target_expression.astype(np.float32)
        regulator_values = regulator_df.values.astype(np.float32)
        
        # Normalize target
        target_norm = (target - target.mean()) / target.std()
        
        # Normalize all regulators
        regulator_means = regulator_values.mean(axis=1)[:, np.newaxis]
        regulator_stds = regulator_values.std(axis=1)[:, np.newaxis]
        regulator_norm = (regulator_values - regulator_means) / regulator_stds
        
        # Compute correlations in batch
        correlations = np.dot(regulator_norm, target_norm) / (len(target) - 1)
        
        # Calculate p-values using t-distribution approximation (faster than pearsonr)
        n = len(target)
        t_stats = correlations * np.sqrt((n - 2) / (1 - correlations**2 + 1e-10))
        p_values = 2 * (1 - stats.t.cdf(np.abs(t_stats), n - 2))
        
        # Create results list
        results = []
        for i, (reg_name, corr, pval) in enumerate(zip(regulator_df.index, correlations, p_values)):
            if not np.isnan(corr):
                if self.permutation_enabled and pval < 0.2:
                    # Only use expensive permutation for promising correlations
                    pval = self._adaptive_permutation_pvalue(target, regulator_values[i], observed_abs=abs(corr))
                results.append((f"{prefix}_{reg_name}", corr, pval))
        
        # Sort by absolute correlation magnitude
        results.sort(key=lambda t: abs(t[1]), reverse=True)
        return results[:top_n]

    def _adaptive_permutation_pvalue(self, x: np.ndarray, y: np.ndarray, observed_abs: float) -> float:
        """Adaptive permutation test for correlation significance.

        Stops early when lower confidence bound of p exceeds alpha or max perms reached.
        """
        rng = np.random.default_rng()
        min_perms = self.permutation_min
        max_perms = self.permutation_max
        batch = self.permutation_batch
        alpha = self.permutation_alpha
        n = 0
        extreme = 0
        while n < max_perms:
            to_run = min(batch, max_perms - n)
            perms = rng.permutation(len(y))
            # Generate permutations one-by-one (len small) for clarity
            for _ in range(to_run):
                y_perm = rng.permutation(y)
                try:
                    corr, _ = pearsonr(x, y_perm)
                    if abs(corr) >= observed_abs - 1e-12:
                        extreme += 1
                except Exception:
                    pass
            n += to_run
            if n >= min_perms:
                # +1 smoothing (conservative)
                p_hat = (extreme + 1) / (n + 1)
                # Early stop if clearly > alpha
                if p_hat > alpha * 2:  # heuristic margin
                    return p_hat
        return (extreme + 1) / (n + 1)

            
    def parallel_analyze_multi_way_interactions(self):
        """Parallel multi-way regulator interaction analysis."""
        print("  🔄 Parallel processing multi-way regulatory interactions...")
        if self.CONFIG[self.mode]['multi_way_gene_limit'] is None:
            sampled_genes = self.datasets['gene'].index.to_list()
        else:
            limit = min(self.CONFIG[self.mode]['multi_way_gene_limit'], len(self.datasets['gene']))
            sampled_genes = np.random.choice(self.datasets['gene'].index, limit, replace=False)
        sampled_genes = list(sampled_genes)
        n_genes = len(sampled_genes)
        if n_genes == 0:
            return pd.DataFrame()
        n_chunks = max(1, min(self.n_jobs * self.chunk_multiplier, n_genes))
        gene_chunks = [chunk for chunk in np.array_split(sampled_genes, n_chunks) if len(chunk) > 0]
        print(f"  📊 Processing {n_genes} genes in {len(gene_chunks)} chunks (mode={self.mode}, mult={self.chunk_multiplier})")
        all_results: List[Dict] = []
        with ProcessPoolExecutor(max_workers=self.n_jobs) as pool:
            futures = {pool.submit(self._process_multi_way_chunk, list(chunk)): i for i, chunk in enumerate(gene_chunks)}
            for fut in as_completed(futures):
                idx = futures[fut]
                try:
                    res = fut.result()
                    all_results.extend(res)
                    print(f"    ✅ Chunk {idx+1}/{len(gene_chunks)}: {len(res)} results")
                except Exception as e:
                    print(f"    ❌ Chunk {idx+1} failed: {e}")
        if not all_results:
            return pd.DataFrame()
        df = pd.DataFrame(all_results)
        print(f"  🎯 Total multi-way interactions: {len(df)}")
        return df

    def _process_multi_way_chunk(self, gene_chunk: List[str]) -> List[Dict]:
        """Process a chunk of genes for multi-way interaction analysis."""
        results = []
        
        for gene in gene_chunk:
            gene_expression = self.datasets['gene'].loc[gene].values
            
            regulators = self._get_top_regulators_vectorized(gene)
            
            if len(regulators) >= 3:
                # Analyze multi-way interactions
                multi_way_result = self._analyze_multi_regulator_interaction(
                    gene_expression, regulators, gene
                )
                
                if multi_way_result:
                    results.append(multi_way_result)
        
        return results
        
    def _get_top_regulators_vectorized(self, gene: str) -> Dict[str, np.ndarray]:
        cfg = self.CONFIG[self.mode]['multi_way_regulators']
        regulators: Dict[str, np.ndarray] = {}
        
        # Use pre-computed gene array if available
        if gene in self.gene_index:
            gene_idx = self.gene_index[gene]
            gene_expression = self.gene_array[gene_idx].astype(np.float32)
        else:
            gene_expression = self.datasets['gene'].loc[gene].values.astype(np.float32)
        
        # Get top correlations for each regulator type
        mirna_corrs = self._get_top_correlations_vectorized(
            gene_expression, self.datasets['mirna'], 'mirna', top_n=cfg['miRNA']
        )
        for name, _, _ in mirna_corrs:
            clean_name = name.replace('mirna_', '')
            if clean_name in self.mirna_index:
                mirna_idx = self.mirna_index[clean_name]
                regulators[name] = self.mirna_array[mirna_idx].astype(np.float32)
            else:
                regulators[name] = self.datasets['mirna'].loc[clean_name].values.astype(np.float32)
        
        lncrna_corrs = self._get_top_correlations_vectorized(
            gene_expression, self.datasets['lncrna'], 'lncrna', top_n=cfg['lncRNA']
        )
        for name, _, _ in lncrna_corrs:
            clean_name = name.replace('lncrna_', '')
            regulators[name] = self.datasets['lncrna'].loc[clean_name].values.astype(np.float32)
        
        meth_corrs = self._get_top_correlations_vectorized(
            gene_expression, self.datasets['methylation'], 'methylation', top_n=cfg['methylation']
        )
        for name, _, _ in meth_corrs:
            clean_name = name.replace('methylation_', '')
            if clean_name in self.methylation_index:
                meth_idx = self.methylation_index[clean_name]
                regulators[name] = self.methylation_array[meth_idx].astype(np.float32)
            else:
                regulators[name] = self.datasets['methylation'].loc[clean_name].values.astype(np.float32)
        
        return regulators

    def _analyze_multi_regulator_interaction(self, gene_expression: np.ndarray, regulators: Dict[str, np.ndarray], 
                                           gene_name: str) -> Dict:
        """Analyze multi-regulator interactions for a specific gene."""
        try:
            # Create feature matrix
            feature_data = {}
            for reg_name, reg_values in regulators.items():
                feature_data[reg_name] = reg_values
            
            feature_data['target'] = gene_expression
            data = pd.DataFrame(feature_data)
            
            # Scale data
            scaler = StandardScaler()
            data_scaled = pd.DataFrame(
                scaler.fit_transform(data),
                columns=data.columns
            )
            
            # Fit models with increasing complexity
            X = data_scaled.drop('target', axis=1)
            y = data_scaled['target']
            
            # Base model (first regulator only)
            base_model = LinearRegression()
            base_model.fit(X.iloc[:, :1], y)
            r2_base = base_model.score(X.iloc[:, :1], y)
            
            # Full model (all regulators)
            full_model = LinearRegression()
            full_model.fit(X, y)
            r2_full = full_model.score(X, y)
            
            # Calculate improvement
            improvement_from_regulators = r2_full - r2_base
            
            # FIXED: Use proper statistical testing instead of arbitrary R^2 threshold
            # Perform nested F-test to compare base vs full model
            from sklearn.metrics import r2_score
            from scipy import stats
            
            # Calculate F-statistic for nested model comparison
            n_samples = len(data_scaled)
            df1 = len(regulators) - 1  # Additional parameters in full model
            df2 = n_samples - len(regulators) - 1  # Degrees of freedom for full model
            
            # FIXED: Add proper error handling for edge cases
            if (df2 > 0 and 
                improvement_from_regulators > 0 and 
                r2_full < 1.0 and  # Prevent division by zero when R² = 1
                (1 - r2_full) > 1e-10):  # Prevent division by very small numbers
                try:
                    f_stat = (improvement_from_regulators / df1) / ((1 - r2_full) / df2)
                    p_value = 1 - stats.f.cdf(f_stat, df1, df2)
                    has_significant_interactions = p_value < 0.05  # Use p-value threshold
                except (ZeroDivisionError, ValueError, RuntimeWarning):
                    # Handle any numerical errors gracefully
                    has_significant_interactions = False
                    p_value = 1.0
            else:
                has_significant_interactions = False
                p_value = 1.0
            
            return {
                'gene': gene_name,
                'n_regulators': len(regulators),
                'r2_base_model': r2_base,
                'r2_full_model': r2_full,
                'improvement_from_regulators': improvement_from_regulators,
                'has_significant_interactions': has_significant_interactions,
                'interaction_p_value': p_value,
                'regulator_types': list(regulators.keys())
            }
            
        except Exception as e:
            print(f"    ⚠️  Error analyzing multi-regulator interaction for {gene_name}: {e}")
            return None
            
    # Context-specific network analysis removed per scope reduction.
        
    def generate_context_visualizations(self):
        print("\n" + "="*60)
        print("GENERATING CONTEXT-DEPENDENT VISUALIZATIONS")
        print("="*60)
        with ThreadPoolExecutor(max_workers=2) as ex:
            futures = [
                ex.submit(self.plot_context_dependent_interactions),
                ex.submit(self.plot_interaction_improvements)
            ]
            for fut in futures: fut.result()
        print(f"✅ All visualizations saved to {self.plots_dir}/")
        
    def plot_context_dependent_interactions(self):
        try:
            if 'context_dependent' not in self.results:
                print("    ⚠️  No context-dependent results available for plotting")
                return
            fig, axes = plt.subplots(1,2, figsize=(15,6))
            meth_mirna = self.results['context_dependent'].get('methylation_mirna_context', pd.DataFrame())
            if not meth_mirna.empty:
                try:
                    axes[0].hist(meth_mirna['improvement_from_interaction'].dropna(), bins=30, alpha=0.7, edgecolor='black')
                    axes[0].set_title('Methylation-miRNA Interaction Improvement')
                    axes[0].set_xlabel('Improvement from Interaction')
                    axes[0].set_ylabel('Frequency')
                    axes[0].axvline(x=0.1, color='red', linestyle='--', alpha=0.7, label='Threshold')
                    axes[0].legend()
                except Exception as e:
                    print(f"    ⚠️  Error plotting methylation-miRNA: {e}")
                    axes[0].text(0.5,0.5,'Plot Error', ha='center', va='center', transform=axes[0].transAxes)
            if not meth_mirna.empty:
                try:
                    strength = meth_mirna['context_strength'].dropna()
                    if len(strength) > 0:
                        axes[1].hist(strength, bins=30, alpha=0.7, edgecolor='black')
                    else:
                        axes[1].text(0.5,0.5,'No Data', ha='center', va='center', transform=axes[1].transAxes)
                    axes[1].set_title('Methylation-miRNA Context Strength')
                except Exception as e:
                    print(f"    ⚠️  Error plotting methylation strength: {e}")
                    axes[1].text(0.5,0.5,'Plot Error', ha='center', va='center', transform=axes[1].transAxes)
            plt.tight_layout()
            try:
                path = os.path.join(self.plots_dir, "context_dependent_interactions.png")
                plt.savefig(path, dpi=300, bbox_inches='tight')
                print(f"    ✅ Saved context-dependent interactions plot: {path}")
            except Exception as e:
                print(f"    ❌ Error saving plot: {e}")
            finally:
                plt.close()
        except Exception as e:
            print(f"    ❌ Error in plot_context_dependent_interactions: {e}")
            try:
                plt.close('all')
            except:  # noqa
                pass
        
        
    def plot_interaction_improvements(self):
        try:
            if 'context_dependent' not in self.results:
                print("    ⚠️  No context-dependent results available for plotting")
                return
            fig, axes = plt.subplots(1,2, figsize=(12,5))
            meth_mirna = self.results['context_dependent'].get('methylation_mirna_context', pd.DataFrame())
            if meth_mirna.empty:
                for ax, title in zip(axes, ['Interaction Improvement','Context Strength']):
                    ax.text(0.5,0.5,'No Data', ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(title)
            else:
                try:
                    axes[0].boxplot([meth_mirna['improvement_from_interaction'].dropna()], labels=['Methylation-miRNA'])
                    axes[0].set_title('Interaction Improvement')
                    axes[0].set_ylabel('ΔR² (Interaction)')
                    axes[0].grid(True, alpha=0.3)
                    axes[1].boxplot([meth_mirna['context_strength'].dropna()], labels=['Methylation-miRNA'])
                    axes[1].set_title('Context Strength')
                    axes[1].set_ylabel('|Corr High - Corr Low|')
                    axes[1].grid(True, alpha=0.3)
                except Exception as e:
                    print(f"    ⚠️  Error plotting improvements: {e}")
                    for ax in axes:
                        ax.text(0.5,0.5,'Plot Error', ha='center', va='center', transform=ax.transAxes)
            plt.tight_layout()
            try:
                path = os.path.join(self.plots_dir, "interaction_improvements.png")
                plt.savefig(path, dpi=300, bbox_inches='tight')
                print(f"    ✅ Saved interaction improvements plot: {path}")
            except Exception as e:
                print(f"    ❌ Error saving plot: {e}")
            finally:
                plt.close()
        except Exception as e:
            print(f"    ❌ Error in plot_interaction_improvements: {e}")
            try:
                plt.close('all')
            except:  # noqa
                pass
        
    def save_context_results(self):
        """Save all context-dependent analysis results."""
        print("\n" + "="*60)
        print("SAVING CONTEXT-DEPENDENT RESULTS")
        print("="*60)
        
        # Save context-dependent results
        if 'context_dependent' in self.results:
            for analysis_type, results in self.results['context_dependent'].items():
                if isinstance(results, pd.DataFrame) and not results.empty:
                    results.to_csv(os.path.join(self.tables_dir, f"{analysis_type}.csv"), index=False)
                    print(f"Saved {analysis_type}: {len(results)} results")
                elif isinstance(results, dict):
                    # Save context networks
                    for context_name, context_data in results.items():
                        if isinstance(context_data, dict):
                            # Save each correlation type
                            for corr_type, corr_data in context_data.items():
                                if isinstance(corr_data, list) and corr_data:
                                    corr_df = pd.DataFrame(corr_data)
                                    corr_df.to_csv(os.path.join(self.tables_dir, f"{context_name}_{corr_type}.csv"), index=False)
                                    print(f"Saved {context_name}_{corr_type}: {len(corr_data)} correlations")
        
        print(f"✅ All context-dependent results saved to {self.tables_dir}/")
        
    def generate_markdown_report(self):
        """Generate comprehensive markdown report (mode-aware)."""
        print("\n" + "="*60)
        print("GENERATING MARKDOWN REPORT")
        print("="*60)
        if self.mode == 'full':
            report_filename = "context_dependent_analysis_report.md"
            title = "Context-Dependent Regulation Analysis Report"
        else:
            report_filename = "subset_context_dependent_analysis_report.md"
            title = "Subset Context-Dependent Regulation Analysis Report"
        report_path = os.path.join(self.reports_dir, report_filename)
        
        with open(report_path, 'w') as f:
            # Header
            f.write(f"# {title}\n\n")
            f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Analysis ID:** {self.timestamp}\n\n")
            
            # Executive Summary
            f.write("## Executive Summary\n\n")
            if self.mode == 'full':
                f.write("This report presents the results of context-dependent regulatory interaction analysis between:\n")
            else:
                f.write("This report presents the results of subset context-dependent regulatory interaction analysis between:\n")
            f.write("- Gene expression\n")
            f.write("- miRNA expression\n")
            f.write("- lncRNA expression\n")
            f.write("- DNA methylation\n\n")
            
            # Analysis Overview
            f.write("## Analysis Overview\n\n")
            f.write(f"- **Parallel Processing:** {self.n_jobs} CPU cores\n")
            f.write(f"- **Available RAM:** {self._get_available_ram():.1f} GB\n")
            f.write(f"- **Datasets Loaded:** {len(self.datasets)} data types\n")
            for data_type, data in self.datasets.items():
                f.write(f"  - {data_type.capitalize()}: {data.shape[1]} samples × {data.shape[0]} features\n")
            f.write("\n")
            
            # Results Summary
            if 'context_dependent' in self.results:
                f.write("## Results Summary\n\n")
                
                # Methylation-miRNA context
                meth_mirna = self.results['context_dependent'].get('methylation_mirna_context', pd.DataFrame())
                if not meth_mirna.empty:
                    f.write("### Methylation-miRNA Context Analysis\n\n")
                    f.write(f"- **Total interactions analyzed:** {len(meth_mirna)}\n")
                    f.write(f"- **Context-dependent interactions:** {meth_mirna['context_dependent'].sum()}\n")
                    f.write(f"- **Mean improvement from interaction:** {meth_mirna['improvement_from_interaction'].mean():.3f}\n")
                    f.write(f"- **Mean context strength:** {meth_mirna['context_strength'].mean():.3f}\n\n")
                    
                    # Top interactions table
                    if len(meth_mirna) > 0:
                        f.write("#### Top 10 Methylation-miRNA Interactions\n\n")
                        f.write("*This table shows genes whose regulation by DNA methylation is context-dependent on miRNA expression levels, indicating epigenetic-regulatory crosstalk.*\n\n")
                        top_interactions = meth_mirna.nlargest(10, 'improvement_from_interaction')
                        f.write("| Rank | Gene | Methylation | miRNA | Improvement | Context Strength |\n")
                        f.write("|------|------|-------------|-------|-------------|------------------|\n")
                        for i, (_, row) in enumerate(top_interactions.iterrows(), 1):
                            f.write(f"| {i} | {row.get('target', 'N/A')} | {row.get('regulator1', 'N/A')} | {row.get('regulator2', 'N/A')} | {row.get('improvement_from_interaction', 0):.3f} | {row.get('context_strength', 0):.3f} |\n")
                        f.write("\n")
                else:
                    f.write("### Methylation-miRNA Context Analysis\n\n")
                    f.write("⚠️ **No methylation-miRNA interactions found.**\n\n")
                
                # lncRNA-miRNA context section removed
                
                # Multi-way interactions
                multi_way = self.results['context_dependent'].get('multi_way_interactions', pd.DataFrame())
                if not multi_way.empty:
                    f.write("### Multi-Way Interaction Analysis\n\n")
                    f.write(f"- **Total genes analyzed:** {len(multi_way)}\n")
                    f.write(f"- **Genes with significant interactions:** {multi_way['has_significant_interactions'].sum()}\n")
                    f.write(f"- **Mean improvement from interactions:** {multi_way['improvement_from_regulators'].mean():.3f}\n\n")
                    
                    # Top multi-way interactions table
                    if len(multi_way) > 0:
                        f.write("#### Top 10 Multi-Way Interactions\n\n")
                        f.write("*This table shows genes with the most complex regulatory patterns, where multiple RNA and epigenetic factors coordinately influence gene expression.*\n\n")
                        top_interactions = multi_way.nlargest(10, 'improvement_from_regulators')
                        f.write("| Rank | Gene | Improvement | Significant Interactions |\n")
                        f.write("|------|------|-------------|------------------------|\n")
                        for i, (_, row) in enumerate(top_interactions.iterrows(), 1):
                            f.write(f"| {i} | {row.get('gene', 'N/A')} | {row.get('improvement_from_regulators', 0):.3f} | {row.get('has_significant_interactions', False)} |\n")
                        f.write("\n")
                else:
                    f.write("### Multi-Way Interaction Analysis\n\n")
                    f.write("⚠️ **No multi-way interactions found.**\n\n")
            else:
                f.write("## Results Summary\n\n")
                f.write("⚠️ **No context-dependent results available yet.**\n\n")
                f.write("This may occur if:\n")
                f.write("- The analysis is still running\n")
                f.write("- No significant interactions were found\n")
                f.write("- An error occurred during analysis\n\n")
                f.write("Please check the console output for any error messages.\n\n")
            
            
            # Data Files
            f.write("## Data Files\n\n")
            f.write("The following data files were generated:\n\n")
            f.write("These CSV files contain the detailed results of context-dependent regulatory analysis, including interaction statistics, p-values, context strengths, and regulatory relationships for each analyzed gene-regulator pair.\n\n")
            
            # List CSV files in tables directory
            if os.path.exists(self.tables_dir):
                csv_files = [f for f in os.listdir(self.tables_dir) if f.endswith('.csv')]
                for csv_file in sorted(csv_files):
                    file_path = os.path.join(self.tables_dir, csv_file)
                    if os.path.exists(file_path):
                        file_size = os.path.getsize(file_path) / 1024  # KB
                        f.write(f"- **{csv_file}** ({file_size:.1f} KB)\n")
                f.write("\n")
            
            # Analysis Parameters
            f.write("## Analysis Parameters\n\n")
            f.write(f"- **Parallel workers:** {self.n_jobs}\n")
            f.write(f"- **Data directory:** {self.data_dir}\n")
            f.write(f"- **Output directory:** {self.output_dir}\n")
            f.write(f"- **Analysis timestamp:** {self.timestamp}\n\n")
            
            # Conclusion
            f.write("## Conclusion\n\n")
            if self.mode == 'full':
                f.write("This analysis successfully identified context-dependent regulatory interactions using optimized parallel processing.\n")
            else:
                f.write("This subset analysis successfully identified context-dependent regulatory interactions using optimized parallel processing.\n")
            f.write("The results provide insights into how different regulatory layers interact in a context-specific manner.\n\n")
            
            # Footer
            f.write("---\n")
            f.write(f"*Report generated by OptimizedContextDependentRegulationAnalysis on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n")
        
        print(f"✅ Markdown report generated: {report_path}")
        
    def generate_html_report(self):
        """Generate HTML report with embedded images (mode-aware)."""
        print("\n" + "="*60)
        print("GENERATING HTML REPORT")
        print("="*60)
        if self.mode == 'full':
            html_report_path = os.path.join(self.reports_dir, "context_dependent_analysis_report.html")
            md_report_path = os.path.join(self.reports_dir, "context_dependent_analysis_report.md")
        else:
            html_report_path = os.path.join(self.reports_dir, "subset_context_dependent_analysis_report.html")
            md_report_path = os.path.join(self.reports_dir, "subset_context_dependent_analysis_report.md")
        
        # Read the markdown content
        if not os.path.exists(md_report_path):
            print("⚠️  Markdown report not found. Generate markdown report first.")
            return
            
        with open(md_report_path, 'r') as f:
            md_content = f.read()
        
        # Convert images to base64 and embed them
        html_content = self._convert_md_to_html_with_embedded_images(md_content)
        
        # Write HTML report
        with open(html_report_path, 'w') as f:
            f.write(html_content)
        
        print(f"✅ HTML report generated: {html_report_path}")
        
    def _convert_md_to_html_with_embedded_images(self, md_content):
        """Convert markdown to HTML with embedded base64 images."""
        # Basic HTML structure
        html_head = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Context-Dependent Regulation Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; line-height: 1.6; max-width: 1200px; margin: 0 auto; padding: 20px; }
        h1, h2, h3, h4 { color: #333; }
        h1 { border-bottom: 2px solid #333; padding-bottom: 10px; }
        h2 { border-bottom: 1px solid #ccc; padding-bottom: 5px; margin-top: 30px; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f5f5f5; font-weight: bold; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        img { max-width: 100%; height: auto; display: block; margin: 20px auto; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        .warning { background-color: #fff3cd; border: 1px solid #ffeaa7; border-radius: 4px; padding: 10px; margin: 10px 0; }
        code { background-color: #f4f4f4; padding: 2px 4px; border-radius: 3px; }
        pre { background-color: #f4f4f4; padding: 10px; border-radius: 5px; overflow-x: auto; }
        ul { margin: 10px 0; }
        li { margin: 5px 0; }
    </style>
</head>
<body>"""
        
        html_footer = """</body>
</html>"""
        
        # Convert markdown to HTML (simple conversion)
        html_body = md_content
        
        # Convert headers
        html_body = html_body.replace('#### ', '<h4>').replace('\n', '</h4>\n', html_body.count('#### '))
        html_body = html_body.replace('### ', '<h3>').replace('\n', '</h3>\n', html_body.count('### '))
        html_body = html_body.replace('## ', '<h2>').replace('\n', '</h2>\n', html_body.count('## '))
        html_body = html_body.replace('# ', '<h1>').replace('\n', '</h1>\n', html_body.count('# '))
        
        # Convert bold text
        import re
        html_body = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', html_body)
        
        # Convert lists
        lines = html_body.split('\n')
        in_list = False
        result_lines = []
        
        for line in lines:
            if line.strip().startswith('- '):
                if not in_list:
                    result_lines.append('<ul>')
                    in_list = True
                result_lines.append(f'<li>{line.strip()[2:]}</li>')
            else:
                if in_list:
                    result_lines.append('</ul>')
                    in_list = False
                result_lines.append(line)
        
        if in_list:
            result_lines.append('</ul>')
            
        html_body = '\n'.join(result_lines)
        
        # Convert paragraphs
        paragraphs = html_body.split('\n\n')
        html_paragraphs = []
        for para in paragraphs:
            para = para.strip()
            if para and not para.startswith('<') and not para.startswith('|'):
                html_paragraphs.append(f'<p>{para}</p>')
            else:
                html_paragraphs.append(para)
        
        html_body = '\n\n'.join(html_paragraphs)
        
        # Convert tables
        table_lines = html_body.split('\n')
        result_lines = []
        in_table = False
        
        for i, line in enumerate(table_lines):
            if '|' in line and line.strip().startswith('|'):
                if not in_table:
                    result_lines.append('<table>')
                    in_table = True
                
                # Check if this is a header separator line
                if '---' in line:
                    continue
                
                # Process table row
                cells = [cell.strip() for cell in line.split('|')[1:-1]]  # Remove empty first/last
                if i < len(table_lines) - 1 and '---' in table_lines[i + 1]:
                    # This is a header row
                    row_html = '<tr>' + ''.join(f'<th>{cell}</th>' for cell in cells) + '</tr>'
                else:
                    # This is a data row
                    row_html = '<tr>' + ''.join(f'<td>{cell}</td>' for cell in cells) + '</tr>'
                
                result_lines.append(row_html)
            else:
                if in_table:
                    result_lines.append('</table>')
                    in_table = False
                result_lines.append(line)
        
        if in_table:
            result_lines.append('</table>')
            
        html_body = '\n'.join(result_lines)
        
        # Convert images to embedded base64
        image_pattern = r'!\[(.*?)\]\((.*?)\)'
        
        def replace_image(match):
            alt_text = match.group(1)
            image_path = match.group(2)
            
            # Convert relative path to absolute
            if image_path.startswith('plots/'):
                full_image_path = os.path.join(self.plots_dir, image_path.replace('plots/', ''))
            else:
                full_image_path = image_path
            
            if os.path.exists(full_image_path):
                try:
                    with open(full_image_path, 'rb') as img_file:
                        img_data = img_file.read()
                        img_base64 = base64.b64encode(img_data).decode('utf-8')
                        return f'<img src="data:image/png;base64,{img_base64}" alt="{alt_text}" title="{alt_text}">'
                except Exception as e:
                    print(f"⚠️  Could not embed image {full_image_path}: {e}")
                    return f'<p class="warning">⚠️ Image not found: {alt_text}</p>'
            else:
                return f'<p class="warning">⚠️ Image not found: {alt_text} ({image_path})</p>'
        
        html_body = re.sub(image_pattern, replace_image, html_body)
        
        # Handle warning messages
        html_body = html_body.replace('⚠️', '<span class="warning">⚠️')
        html_body = html_body.replace('**\n\n', '**</span>\n\n')
        
        return html_head + '\n' + html_body + '\n' + html_footer
        
    def print_context_summary(self):
        """Print summary of context-dependent analysis."""
        print("\n" + "="*60)
        print("CONTEXT-DEPENDENT ANALYSIS SUMMARY")
        print("="*60)
        
        if 'context_dependent' not in self.results:
            print("No context-dependent results available.")
            return
        
        # Summary of methylation-miRNA context
        meth_mirna = self.results['context_dependent']['methylation_mirna_context']
        if not meth_mirna.empty:
            print(f"\nMETHYLATION-MIRNA CONTEXT ANALYSIS:")
            print(f"  Total interactions analyzed: {len(meth_mirna)}")
            print(f"  Context-dependent interactions: {meth_mirna['context_dependent'].sum()}")
            print(f"  Mean improvement from interaction: {meth_mirna['improvement_from_interaction'].mean():.3f}")
            print(f"  Mean context strength: {meth_mirna['context_strength'].mean():.3f}")
        
    # Summary of multi-way interactions
        multi_way = self.results['context_dependent']['multi_way_interactions']
        if not multi_way.empty:
            print(f"\nMULTI-WAY INTERACTION ANALYSIS:")
            print(f"  Total genes analyzed: {len(multi_way)}")
            print(f"  Genes with significant interactions: {multi_way['has_significant_interactions'].sum()}")
            print(f"  Mean improvement from interactions: {multi_way['improvement_from_regulators'].mean():.3f}")
    # lncRNA-miRNA and context-specific network summaries removed
    
    def run_complete_context_analysis(self):
        """Run the complete optimized context-dependent analysis (mode-aware)."""
        start_time = time.time()
        
        print("="*80)
        print(f"🚀 OPTIMIZED CONTEXT-DEPENDENT REGULATION ANALYSIS ({self.mode.upper()} MODE)")
        print("="*80)
        print("This optimized analysis utilizes:")
        print(f"  • {self.n_jobs} parallel CPU cores")
        print(f"  • Vectorized numpy/pandas operations")
        print(f"  • Memory-efficient batch processing")
        print(f"  • Concurrent data loading and processing")
        print("="*80)
        
        # Run context-dependent analysis
        self.analyze_context_dependent_regulation()
        
        # Generate visualizations
        self.generate_context_visualizations()
        
        # Save results
        self.save_context_results()
        
        # Generate markdown report
        self.generate_markdown_report()
        
        # Generate HTML report (without images)
        self.generate_html_report()
        
        # Print summary
        self.print_context_summary()
        
        total_time = time.time() - start_time
        print(f"\n⏱️  Total analysis time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
        print("\n" + "="*80)
        print("🎉 OPTIMIZED CONTEXT-DEPENDENT ANALYSIS COMPLETED SUCCESSFULLY!")
        print("="*80)
        # Gracefully shutdown shared pool
    # No persistent pool to shutdown

def main(mode: str = None):
    """Main entry point with optional mode parameter and CLI prompt."""
    if mode is None:
        import argparse
        parser = argparse.ArgumentParser(description="Unified context-dependent regulation analysis")
        parser.add_argument('--mode', choices=['full','subset'], help="Run mode: full or subset", default=None)
        parser.add_argument('--n-jobs', type=int, default=None, help="Number of parallel workers")
        parser.add_argument('--output-root', type=str, default='output', help="Base output directory relative to repo root (default: output)")
        parser.add_argument('--perm', action='store_true', help='Enable adaptive permutation p-values for correlations')
        parser.add_argument('--perm-min', type=int, default=100, help='Minimum permutations (default 100)')
        parser.add_argument('--perm-max', type=int, default=1000, help='Maximum permutations (default 1000)')
        parser.add_argument('--perm-alpha', type=float, default=0.05, help='Alpha for adaptive early stopping (default 0.05)')
        args = parser.parse_args()
        mode = args.mode
        n_jobs = args.n_jobs
        output_root = args.output_root
    else:
        n_jobs = None
        output_root = 'output'
    if mode is None:
        try:
            user_input = input("Run mode? (full/subset) [full]: ").strip().lower()
            mode = user_input if user_input in ('full','subset') else 'full'
        except EOFError:
            mode = 'full'
    analysis = OptimizedContextDependentRegulationAnalysis(mode=mode, n_jobs=n_jobs, output_root=output_root)
    if mode is None:
        mode = analysis.mode  # ensure mode set
    # Apply permutation settings if provided
    if 'args' in locals() and args.perm:
        analysis.permutation_enabled = True
        analysis.permutation_min = max(10, args.perm_min)
        analysis.permutation_max = max(analysis.permutation_min, args.perm_max)
        analysis.permutation_alpha = args.perm_alpha
        analysis.permutation_batch = max(10, analysis.permutation_min // 10)
        print(f"🔁 Permutation testing enabled: min={analysis.permutation_min}, max={analysis.permutation_max}, alpha={analysis.permutation_alpha}, batch={analysis.permutation_batch}")
    analysis.run_complete_context_analysis()

if __name__ == "__main__":
    main()
