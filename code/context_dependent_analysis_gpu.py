#!/usr/bin/env python3
"""
GPU-OPTIMIZED Context-Dependent Regulation Analysis

This GPU-optimized script identifies context-dependent regulatory interactions using:
- NVIDIA RAPIDS cuDF for GPU-accelerated DataFrames
- CuPy for GPU-accelerated NumPy operations
- cuML for GPU-accelerated machine learning
- Batch GPU processing with memory management

Requirements:
- NVIDIA GPU with CUDA support
- RAPIDS suite (cudf, cuml, cupy)
- Fallback to CPU if GPU is not available

Methods used:
- Interaction term analysis (GPU-accelerated)
- Conditional correlation analysis (GPU-vectorized)
- Multi-variable regression with interaction terms (cuML)
- Context-specific regulatory network inference (GPU-optimized)
"""

import os
import gc
import time
import base64
import warnings
from datetime import datetime
from typing import Dict, List, Tuple, Any, Optional

# Standard imports (always available)
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless environments
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr

warnings.filterwarnings('ignore')

# Try to import GPU libraries
GPU_AVAILABLE = False
try:
    import cupy as cp
    import cudf
    from cuml.linear_model import LinearRegression as cuLinearRegression
    from cuml.preprocessing import StandardScaler as cuStandardScaler
    GPU_AVAILABLE = True
    print("🚀 GPU acceleration available via RAPIDS")
except ImportError:
    print("⚠️  RAPIDS not available. Install with: conda install -c rapidsai -c conda-forge -c nvidia rapids")
    print("   Falling back to CPU-optimized implementation...")
    # CPU fallbacks
    cp = np  # Use numpy as fallback
    cudf = pd  # Use pandas as fallback
    from sklearn.linear_model import LinearRegression as cuLinearRegression
    from sklearn.preprocessing import StandardScaler as cuStandardScaler

# CPU fallback imports (always needed for some operations)
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing as mp

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")


class GPUContextDependentRegulationAnalysis:
    """GPU-optimized context-dependent regulatory analysis."""
    
    def __init__(self, data_dir: str = "data/cleaned_datasets", 
                 gpu_memory_fraction: float = 0.8,
                 batch_size: int = 1000):
        """
        Initialize the GPU-optimized context-dependent analysis.
        
        Args:
            data_dir: Directory containing cleaned datasets
            gpu_memory_fraction: Fraction of GPU memory to use (0.0-1.0)
            batch_size: Number of genes to process per GPU batch
        """
        # Get workspace root directory (parent of code directory)
        self.workspace_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        # Set absolute paths relative to workspace root
        self.data_dir = os.path.join(self.workspace_root, data_dir)
        self.datasets = {}
        self.gpu_datasets = {}  # GPU versions of datasets
        self.results = {}
        
        # GPU configuration
        self.gpu_available = GPU_AVAILABLE
        self.gpu_memory_fraction = gpu_memory_fraction
        self.batch_size = batch_size
        
        # Configure GPU if available
        if self.gpu_available:
            self._configure_gpu()
        
        # CPU fallback workers
        self.n_jobs = mp.cpu_count()
        
        # Create timestamped output directory
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir = os.path.join(
            self.workspace_root, 
            f"output/context_dependent_analysis_gpu_{self.timestamp}"
        )
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Create subdirectories for different output types
        self.plots_dir = os.path.join(self.output_dir, "plots")
        self.tables_dir = os.path.join(self.output_dir, "tables")
        self.reports_dir = os.path.join(self.output_dir, "reports")
        
        os.makedirs(self.plots_dir, exist_ok=True)
        os.makedirs(self.tables_dir, exist_ok=True)
        os.makedirs(self.reports_dir, exist_ok=True)
        
        print(f"{'🚀 GPU' if self.gpu_available else '💻 CPU'} Context-Dependent Analysis")
        print(f"📁 Output directory: {self.output_dir}")
        
        if self.gpu_available:
            print(f"🎮 GPU Memory Fraction: {self.gpu_memory_fraction}")
            print(f"📦 Batch Size: {self.batch_size} genes per batch")
        else:
            print(f"🔧 CPU Workers: {self.n_jobs}")
        
        self.load_datasets()
    
    def _configure_gpu(self):
        """Configure GPU settings and memory pool."""
        if not self.gpu_available:
            return
            
        try:
            # Get GPU device info
            device = cp.cuda.Device()
            self.gpu_device_id = device.id
            self.gpu_name = cp.cuda.runtime.getDeviceProperties(device.id)['name'].decode()
            self.gpu_memory_total = device.mem_info[1] / (1024**3)  # GB
            self.gpu_memory_available = device.mem_info[0] / (1024**3)  # GB
            
            print(f"🎮 GPU Device: {self.gpu_name}")
            print(f"💾 GPU Memory: {self.gpu_memory_available:.1f} GB available / {self.gpu_memory_total:.1f} GB total")
            
            # Configure memory pool for better memory management
            mempool = cp.get_default_memory_pool()
            mempool.set_limit(size=int(self.gpu_memory_total * self.gpu_memory_fraction * 1024**3))
            
        except Exception as e:
            print(f"⚠️  GPU configuration warning: {e}")
            self.gpu_memory_total = 8.0  # Default
            self.gpu_memory_available = 6.0
    
    def _get_available_ram(self) -> float:
        """Get available RAM in GB."""
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemAvailable:'):
                        return int(line.split()[1]) / (1024**2)
        except:
            pass
        
        try:
            import subprocess
            result = subprocess.run(['sysctl', '-n', 'hw.memsize'], 
                                  capture_output=True, text=True, check=True)
            return int(result.stdout.strip()) / (1024**3)
        except:
            pass
        
        try:
            import psutil
            return psutil.virtual_memory().total / (1024**3)
        except ImportError:
            pass
        
        return 16.0
    
    def load_datasets(self):
        """Load all cleaned datasets with parallel processing."""
        print("📂 Loading cleaned datasets...")
        
        # Load datasets in parallel using ThreadPoolExecutor for I/O
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
        
        # Transfer data to GPU if available
        if self.gpu_available:
            self._transfer_to_gpu()
        else:
            self._precompute_cpu_arrays()
    
    def _transfer_to_gpu(self):
        """Transfer datasets to GPU memory using cuDF."""
        print("🚀 Transferring data to GPU...")
        
        try:
            for data_type, df in self.datasets.items():
                # Convert pandas DataFrame to cuDF DataFrame
                self.gpu_datasets[data_type] = cudf.DataFrame.from_pandas(df)
                print(f"  ✅ {data_type}: transferred to GPU")
            
            # Also create CuPy arrays for fast numerical operations
            self.gene_gpu = cp.asarray(self.datasets['gene'].values.astype(np.float32))
            self.lncrna_gpu = cp.asarray(self.datasets['lncrna'].values.astype(np.float32))
            self.mirna_gpu = cp.asarray(self.datasets['mirna'].values.astype(np.float32))
            self.methylation_gpu = cp.asarray(self.datasets['methylation'].values.astype(np.float32))
            
            print("✅ GPU data transfer complete")
            
        except Exception as e:
            print(f"⚠️  GPU transfer failed: {e}")
            print("   Falling back to CPU arrays...")
            self._precompute_cpu_arrays()
    
    def _precompute_cpu_arrays(self):
        """Pre-compute numpy arrays for CPU processing."""
        print("🔢 Pre-computing CPU data arrays...")
        
        self.gene_array = self.datasets['gene'].values.astype(np.float32)
        self.lncrna_array = self.datasets['lncrna'].values.astype(np.float32)
        self.mirna_array = self.datasets['mirna'].values.astype(np.float32)
        self.methylation_array = self.datasets['methylation'].values.astype(np.float32)
        
        print("✅ CPU data arrays pre-computed")
    
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

    # =========================================================================
    # GPU-ACCELERATED CORRELATION FUNCTIONS
    # =========================================================================
    
    def gpu_pearson_correlation_batch(self, x_batch: cp.ndarray, y_batch: cp.ndarray) -> Tuple[cp.ndarray, cp.ndarray]:
        """
        Compute Pearson correlations for batches of vectors using GPU.
        
        Args:
            x_batch: Shape (n_features_x, n_samples)
            y_batch: Shape (n_features_y, n_samples)
        
        Returns:
            correlations: Shape (n_features_x, n_features_y)
            p_values: Shape (n_features_x, n_features_y)
        """
        if not self.gpu_available:
            return self._cpu_pearson_correlation_batch(x_batch, y_batch)
        
        n_samples = x_batch.shape[1]
        
        # Center the data (subtract mean)
        x_centered = x_batch - cp.mean(x_batch, axis=1, keepdims=True)
        y_centered = y_batch - cp.mean(y_batch, axis=1, keepdims=True)
        
        # Compute standard deviations
        x_std = cp.std(x_batch, axis=1, keepdims=True)
        y_std = cp.std(y_batch, axis=1, keepdims=True)
        
        # Avoid division by zero
        x_std = cp.where(x_std == 0, 1e-10, x_std)
        y_std = cp.where(y_std == 0, 1e-10, y_std)
        
        # Normalize
        x_norm = x_centered / x_std
        y_norm = y_centered / y_std
        
        # Compute correlation matrix via matrix multiplication
        correlations = cp.dot(x_norm, y_norm.T) / n_samples
        
        # Compute p-values using t-distribution approximation
        # t = r * sqrt((n-2)/(1-r^2))
        n = n_samples
        r_squared = correlations ** 2
        r_squared = cp.clip(r_squared, 0, 0.9999)  # Avoid division by zero
        
        t_stat = correlations * cp.sqrt((n - 2) / (1 - r_squared + 1e-10))
        
        # Two-tailed p-value using scipy on CPU (no GPU equivalent)
        t_stat_cpu = cp.asnumpy(cp.abs(t_stat))
        p_values = 2 * stats.t.sf(t_stat_cpu, n - 2)
        p_values = cp.asarray(p_values)
        
        return correlations, p_values
    
    def _cpu_pearson_correlation_batch(self, x_batch: np.ndarray, y_batch: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """CPU fallback for batch Pearson correlation."""
        n_samples = x_batch.shape[1]
        
        # Center the data
        x_centered = x_batch - np.mean(x_batch, axis=1, keepdims=True)
        y_centered = y_batch - np.mean(y_batch, axis=1, keepdims=True)
        
        # Compute standard deviations
        x_std = np.std(x_batch, axis=1, keepdims=True)
        y_std = np.std(y_batch, axis=1, keepdims=True)
        
        # Avoid division by zero
        x_std = np.where(x_std == 0, 1e-10, x_std)
        y_std = np.where(y_std == 0, 1e-10, y_std)
        
        # Normalize
        x_norm = x_centered / x_std
        y_norm = y_centered / y_std
        
        # Compute correlation matrix
        correlations = np.dot(x_norm, y_norm.T) / n_samples
        
        # Compute p-values
        n = n_samples
        r_squared = np.clip(correlations ** 2, 0, 0.9999)
        t_stat = correlations * np.sqrt((n - 2) / (1 - r_squared + 1e-10))
        p_values = 2 * stats.t.sf(np.abs(t_stat), n - 2)
        
        return correlations, p_values
    
    def gpu_linear_regression_batch(self, X: cp.ndarray, y: cp.ndarray) -> Tuple[cp.ndarray, float]:
        """
        Perform linear regression using GPU.
        
        Args:
            X: Feature matrix (n_samples, n_features)
            y: Target vector (n_samples,)
        
        Returns:
            coefficients: Regression coefficients
            r2_score: R-squared score
        """
        if not self.gpu_available:
            return self._cpu_linear_regression(X, y)
        
        try:
            # Use cuML LinearRegression
            model = cuLinearRegression()
            
            # Convert to cuDF if needed
            if isinstance(X, cp.ndarray):
                X_df = cudf.DataFrame(cp.asnumpy(X))
                y_series = cudf.Series(cp.asnumpy(y))
            else:
                X_df = cudf.DataFrame(X)
                y_series = cudf.Series(y)
            
            model.fit(X_df, y_series)
            predictions = model.predict(X_df)
            
            # Calculate R²
            ss_res = ((y_series - predictions) ** 2).sum()
            ss_tot = ((y_series - y_series.mean()) ** 2).sum()
            r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
            
            return cp.asarray(model.coef_), float(r2)
            
        except Exception as e:
            # Fallback to CPU
            return self._cpu_linear_regression(X, y)
    
    def _cpu_linear_regression(self, X: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, float]:
        """CPU fallback for linear regression."""
        if hasattr(X, 'get'):  # CuPy array
            X = X.get()
        if hasattr(y, 'get'):
            y = y.get()
        
        model = LinearRegression()
        model.fit(X, y)
        r2 = model.score(X, y)
        
        return np.array(model.coef_), r2

    # =========================================================================
    # GPU-ACCELERATED ANALYSIS METHODS
    # =========================================================================
    
    def analyze_context_dependent_regulation(self):
        """Main analysis for context-dependent regulation using GPU acceleration."""
        print("\n" + "="*80)
        print(f"{'🚀 GPU' if self.gpu_available else '💻 CPU'} CONTEXT-DEPENDENT REGULATION ANALYSIS")
        print("="*80)
        
        start_time = time.time()
        
        # 1. Analyze methylation-gene interactions dependent on miRNA levels
        print("\n1. 🔄 Analyzing methylation-gene interactions dependent on miRNA levels...")
        methylation_mirna_context = self.gpu_analyze_methylation_mirna_context()
        
        # 2. Analyze lncRNA-gene interactions dependent on miRNA levels
        print("\n2. 🔄 Analyzing lncRNA-gene interactions dependent on miRNA levels...")
        lncrna_mirna_context = self.gpu_analyze_lncrna_mirna_context()
        
        # 3. Analyze multi-way regulatory interactions
        print("\n3. 🔄 Analyzing multi-way regulatory interactions...")
        multi_way_interactions = self.gpu_analyze_multi_way_interactions()
        
        # 4. Context-specific regulatory network inference
        print("\n4. 🔄 Inferring context-specific regulatory networks...")
        context_networks = self.gpu_infer_context_specific_networks()
        
        # Store results
        self.results['context_dependent'] = {
            'methylation_mirna_context': methylation_mirna_context,
            'lncrna_mirna_context': lncrna_mirna_context,
            'multi_way_interactions': multi_way_interactions,
            'context_networks': context_networks
        }
        
        total_time = time.time() - start_time
        print(f"\n⏱️  Total analysis time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
        print("✅ Context-dependent analysis completed!")
        
        # Clear GPU memory
        if self.gpu_available:
            self._clear_gpu_memory()
    
    def _clear_gpu_memory(self):
        """Clear GPU memory pool."""
        if self.gpu_available:
            try:
                mempool = cp.get_default_memory_pool()
                mempool.free_all_blocks()
                gc.collect()
            except:
                pass
    
    def gpu_analyze_methylation_mirna_context(self) -> pd.DataFrame:
        """Analyze methylation-gene interactions dependent on miRNA levels using GPU."""
        all_results = []
        
        gene_index = self.datasets['gene'].index.tolist()
        n_genes = len(gene_index)
        
        # Process in batches to manage GPU memory
        n_batches = (n_genes + self.batch_size - 1) // self.batch_size
        print(f"  📊 Processing {n_genes} genes in {n_batches} batches...")
        
        for batch_idx in range(n_batches):
            start_idx = batch_idx * self.batch_size
            end_idx = min((batch_idx + 1) * self.batch_size, n_genes)
            batch_genes = gene_index[start_idx:end_idx]
            
            batch_results = self._process_methylation_mirna_batch_gpu(batch_genes)
            all_results.extend(batch_results)
            
            # Progress update
            if (batch_idx + 1) % 10 == 0 or batch_idx == n_batches - 1:
                print(f"    ✅ Batch {batch_idx + 1}/{n_batches} completed: {len(all_results)} results so far")
            
            # Clear GPU memory periodically
            if self.gpu_available and (batch_idx + 1) % 20 == 0:
                self._clear_gpu_memory()
        
        if all_results:
            results_df = pd.DataFrame(all_results)
            print(f"  🎯 Total methylation-miRNA context interactions: {len(results_df)}")
            return results_df
        else:
            return pd.DataFrame()
    
    def _process_methylation_mirna_batch_gpu(self, gene_batch: List[str]) -> List[Dict]:
        """Process a batch of genes for methylation-miRNA context analysis using GPU."""
        results = []
        
        for gene in gene_batch:
            gene_expression = self.datasets['gene'].loc[gene].values.astype(np.float32)
            
            if self.gpu_available:
                gene_gpu = cp.asarray(gene_expression)
            else:
                gene_gpu = gene_expression
            
            # Get top miRNAs for this gene
            mirna_corrs = self._get_top_correlations_gpu(
                gene_gpu, 'mirna', top_n=25
            )
            
            # Get top methylation sites for this gene
            meth_corrs = self._get_top_correlations_gpu(
                gene_gpu, 'methylation', top_n=50
            )
            
            # Analyze interactions for top regulators
            for mirna_name, mirna_corr, mirna_pval in mirna_corrs[:10]:
                for meth_name, meth_corr, meth_pval in meth_corrs[:15]:
                    interaction_result = self._analyze_interaction_gpu(
                        gene, gene_expression, 
                        meth_name, mirna_name,
                        'methylation', 'mirna',
                        'methylation_mirna'
                    )
                    if interaction_result:
                        results.append(interaction_result)
        
        return results
    
    def _get_top_correlations_gpu(self, target_data, regulator_type: str, 
                                   top_n: int = 10) -> List[Tuple[str, float, float]]:
        """Get top correlations using GPU-accelerated operations."""
        correlations = []
        
        # Get regulator data
        regulator_df = self.datasets[regulator_type]
        
        if self.gpu_available:
            # GPU-accelerated batch correlation
            target_gpu = target_data if isinstance(target_data, cp.ndarray) else cp.asarray(target_data)
            
            if regulator_type == 'mirna':
                regulator_gpu = self.mirna_gpu
            elif regulator_type == 'methylation':
                regulator_gpu = self.methylation_gpu
            elif regulator_type == 'lncrna':
                regulator_gpu = self.lncrna_gpu
            else:
                regulator_gpu = cp.asarray(regulator_df.values.astype(np.float32))
            
            # Compute all correlations at once
            target_batch = target_gpu.reshape(1, -1)
            corr_matrix, pval_matrix = self.gpu_pearson_correlation_batch(target_batch, regulator_gpu)
            
            # Extract correlations and p-values
            corrs = cp.asnumpy(corr_matrix[0])
            pvals = cp.asnumpy(pval_matrix[0])
            
            # Filter by significance and collect results
            for i, regulator_name in enumerate(regulator_df.index):
                if pvals[i] < 0.1:
                    correlations.append((
                        f"{regulator_type}_{regulator_name}",
                        float(corrs[i]),
                        float(pvals[i])
                    ))
        else:
            # CPU fallback - vectorized numpy
            target_array = target_data if isinstance(target_data, np.ndarray) else np.asarray(target_data)
            regulator_array = regulator_df.values.astype(np.float32)
            
            target_batch = target_array.reshape(1, -1)
            corr_matrix, pval_matrix = self._cpu_pearson_correlation_batch(target_batch, regulator_array)
            
            corrs = corr_matrix[0]
            pvals = pval_matrix[0]
            
            for i, regulator_name in enumerate(regulator_df.index):
                if pvals[i] < 0.1:
                    correlations.append((
                        f"{regulator_type}_{regulator_name}",
                        float(corrs[i]),
                        float(pvals[i])
                    ))
        
        # Sort by absolute correlation and return top N
        correlations.sort(key=lambda x: abs(x[1]), reverse=True)
        return correlations[:top_n]
    
    def _analyze_interaction_gpu(self, gene_name: str, gene_expression: np.ndarray,
                                  reg1_name: str, reg2_name: str,
                                  reg1_type: str, reg2_type: str,
                                  interaction_type: str) -> Optional[Dict]:
        """Analyze interaction between two regulators for a specific gene using GPU."""
        try:
            # Get regulator data
            reg1_data = self.datasets[reg1_type].loc[reg1_name.replace(f'{reg1_type}_', '')].values.astype(np.float32)
            reg2_data = self.datasets[reg2_type].loc[reg2_name.replace(f'{reg2_type}_', '')].values.astype(np.float32)
            
            if self.gpu_available:
                # GPU-accelerated processing
                return self._analyze_interaction_gpu_core(
                    gene_name, gene_expression, reg1_data, reg2_data,
                    reg1_name, reg2_name, interaction_type
                )
            else:
                # CPU fallback
                return self._analyze_interaction_cpu(
                    gene_name, gene_expression, reg1_data, reg2_data,
                    reg1_name, reg2_name, interaction_type
                )
                
        except Exception as e:
            return None
    
    def _analyze_interaction_gpu_core(self, gene_name: str, gene_expression: np.ndarray,
                                       reg1_data: np.ndarray, reg2_data: np.ndarray,
                                       reg1_name: str, reg2_name: str,
                                       interaction_type: str) -> Optional[Dict]:
        """GPU core implementation of interaction analysis."""
        try:
            # Transfer to GPU
            target_gpu = cp.asarray(gene_expression)
            reg1_gpu = cp.asarray(reg1_data)
            reg2_gpu = cp.asarray(reg2_data)
            interaction_gpu = reg1_gpu * reg2_gpu
            
            # Stack features
            X1 = cp.column_stack([reg1_gpu])
            X2 = cp.column_stack([reg1_gpu, reg2_gpu])
            X3 = cp.column_stack([reg1_gpu, reg2_gpu, interaction_gpu])
            
            # Standardize using GPU
            def gpu_standardize(X):
                mean = cp.mean(X, axis=0)
                std = cp.std(X, axis=0)
                std = cp.where(std == 0, 1, std)
                return (X - mean) / std
            
            target_scaled = gpu_standardize(target_gpu.reshape(-1, 1)).ravel()
            X1_scaled = gpu_standardize(X1)
            X2_scaled = gpu_standardize(X2)
            X3_scaled = gpu_standardize(X3)
            
            # Fit models using cuML or numpy fallback
            _, r2_1 = self.gpu_linear_regression_batch(cp.asnumpy(X1_scaled), cp.asnumpy(target_scaled))
            _, r2_2 = self.gpu_linear_regression_batch(cp.asnumpy(X2_scaled), cp.asnumpy(target_scaled))
            _, r2_3 = self.gpu_linear_regression_batch(cp.asnumpy(X3_scaled), cp.asnumpy(target_scaled))
            
            # Calculate improvements
            improvement_from_reg2 = r2_2 - r2_1
            improvement_from_interaction = r2_3 - r2_2
            
            # Statistical testing
            n_samples = len(gene_expression)
            df1 = 1
            df2 = n_samples - 4
            
            if (df2 > 0 and improvement_from_interaction > 0 and 
                r2_3 < 1.0 and (1 - r2_3) > 1e-10):
                try:
                    f_stat = (improvement_from_interaction / df1) / ((1 - r2_3) / df2)
                    p_value = 1 - stats.f.cdf(f_stat, df1, df2)
                    context_dependent = p_value < 0.05
                except:
                    context_dependent = False
                    p_value = 1.0
            else:
                context_dependent = False
                p_value = 1.0
            
            # Conditional correlations on GPU
            reg2_scaled = gpu_standardize(reg2_gpu.reshape(-1, 1)).ravel()
            high_mask = reg2_scaled > 0.5
            low_mask = reg2_scaled < -0.5
            
            high_count = int(cp.sum(high_mask))
            low_count = int(cp.sum(low_mask))
            
            if high_count > 2 and low_count > 2:
                # Transfer to CPU for correlation calculation
                target_np = cp.asnumpy(target_scaled)
                reg1_scaled_np = cp.asnumpy(gpu_standardize(reg1_gpu.reshape(-1, 1)).ravel())
                high_mask_np = cp.asnumpy(high_mask)
                low_mask_np = cp.asnumpy(low_mask)
                
                corr_high, _ = pearsonr(target_np[high_mask_np], reg1_scaled_np[high_mask_np])
                corr_low, _ = pearsonr(target_np[low_mask_np], reg1_scaled_np[low_mask_np])
                context_strength = abs(corr_high - corr_low)
            else:
                corr_high = corr_low = context_strength = np.nan
            
            if pd.isna(corr_high) or pd.isna(corr_low):
                context_direction = 'NA'
            else:
                context_direction = 'positive' if corr_high > corr_low else 'negative'
            
            return {
                'interaction_type': interaction_type,
                'target': gene_name,
                'regulator1': reg1_name,
                'regulator2': reg2_name,
                'r2_regulator1_only': r2_1,
                'r2_regulator1_regulator2': r2_2,
                'r2_with_interaction': r2_3,
                'improvement_from_regulator2': improvement_from_reg2,
                'improvement_from_interaction': improvement_from_interaction,
                'context_dependent': context_dependent,
                'corr_high_regulator2': corr_high,
                'corr_low_regulator2': corr_low,
                'context_strength': context_strength,
                'context_direction': context_direction,
                'interaction_p_value': p_value
            }
            
        except Exception as e:
            return None
    
    def _analyze_interaction_cpu(self, gene_name: str, gene_expression: np.ndarray,
                                  reg1_data: np.ndarray, reg2_data: np.ndarray,
                                  reg1_name: str, reg2_name: str,
                                  interaction_type: str) -> Optional[Dict]:
        """CPU fallback for interaction analysis."""
        try:
            # Create interaction dataset
            data = pd.DataFrame({
                'target': gene_expression,
                'regulator1': reg1_data,
                'regulator2': reg2_data,
                'interaction': reg1_data * reg2_data
            })
            
            # Scale data
            scaler = StandardScaler()
            data_scaled = pd.DataFrame(
                scaler.fit_transform(data),
                columns=data.columns
            )
            
            # Fit models
            model1 = LinearRegression()
            model1.fit(data_scaled[['regulator1']], data_scaled['target'])
            r2_1 = model1.score(data_scaled[['regulator1']], data_scaled['target'])
            
            model2 = LinearRegression()
            model2.fit(data_scaled[['regulator1', 'regulator2']], data_scaled['target'])
            r2_2 = model2.score(data_scaled[['regulator1', 'regulator2']], data_scaled['target'])
            
            model3 = LinearRegression()
            model3.fit(data_scaled[['regulator1', 'regulator2', 'interaction']], data_scaled['target'])
            r2_3 = model3.score(data_scaled[['regulator1', 'regulator2', 'interaction']], data_scaled['target'])
            
            # Calculate improvements
            improvement_from_reg2 = r2_2 - r2_1
            improvement_from_interaction = r2_3 - r2_2
            
            # Statistical testing
            n_samples = len(data_scaled)
            df1 = 1
            df2 = n_samples - 4
            
            if (df2 > 0 and improvement_from_interaction > 0 and 
                r2_3 < 1.0 and (1 - r2_3) > 1e-10):
                try:
                    f_stat = (improvement_from_interaction / df1) / ((1 - r2_3) / df2)
                    p_value = 1 - stats.f.cdf(f_stat, df1, df2)
                    context_dependent = p_value < 0.05
                except:
                    context_dependent = False
                    p_value = 1.0
            else:
                context_dependent = False
                p_value = 1.0
            
            # Conditional correlations
            high_mask = data_scaled['regulator2'] > 0.5
            low_mask = data_scaled['regulator2'] < -0.5
            
            if high_mask.sum() > 2 and low_mask.sum() > 2:
                corr_high, _ = pearsonr(
                    data_scaled.loc[high_mask, 'target'],
                    data_scaled.loc[high_mask, 'regulator1']
                )
                corr_low, _ = pearsonr(
                    data_scaled.loc[low_mask, 'target'],
                    data_scaled.loc[low_mask, 'regulator1']
                )
                context_strength = abs(corr_high - corr_low)
            else:
                corr_high = corr_low = context_strength = np.nan
            
            if pd.isna(corr_high) or pd.isna(corr_low):
                context_direction = 'NA'
            else:
                context_direction = 'positive' if corr_high > corr_low else 'negative'
            
            return {
                'interaction_type': interaction_type,
                'target': gene_name,
                'regulator1': reg1_name,
                'regulator2': reg2_name,
                'r2_regulator1_only': r2_1,
                'r2_regulator1_regulator2': r2_2,
                'r2_with_interaction': r2_3,
                'improvement_from_regulator2': improvement_from_reg2,
                'improvement_from_interaction': improvement_from_interaction,
                'context_dependent': context_dependent,
                'corr_high_regulator2': corr_high,
                'corr_low_regulator2': corr_low,
                'context_strength': context_strength,
                'context_direction': context_direction,
                'interaction_p_value': p_value
            }
            
        except Exception as e:
            return None
    
    def gpu_analyze_lncrna_mirna_context(self) -> pd.DataFrame:
        """Analyze lncRNA-gene interactions dependent on miRNA levels using GPU."""
        all_results = []
        
        gene_index = self.datasets['gene'].index.tolist()
        n_genes = len(gene_index)
        
        n_batches = (n_genes + self.batch_size - 1) // self.batch_size
        print(f"  📊 Processing {n_genes} genes in {n_batches} batches...")
        
        for batch_idx in range(n_batches):
            start_idx = batch_idx * self.batch_size
            end_idx = min((batch_idx + 1) * self.batch_size, n_genes)
            batch_genes = gene_index[start_idx:end_idx]
            
            batch_results = self._process_lncrna_mirna_batch_gpu(batch_genes)
            all_results.extend(batch_results)
            
            if (batch_idx + 1) % 10 == 0 or batch_idx == n_batches - 1:
                print(f"    ✅ Batch {batch_idx + 1}/{n_batches} completed: {len(all_results)} results so far")
            
            if self.gpu_available and (batch_idx + 1) % 20 == 0:
                self._clear_gpu_memory()
        
        if all_results:
            results_df = pd.DataFrame(all_results)
            print(f"  🎯 Total lncRNA-miRNA context interactions: {len(results_df)}")
            return results_df
        else:
            return pd.DataFrame()
    
    def _process_lncrna_mirna_batch_gpu(self, gene_batch: List[str]) -> List[Dict]:
        """Process a batch of genes for lncRNA-miRNA context analysis using GPU."""
        results = []
        
        for gene in gene_batch:
            gene_expression = self.datasets['gene'].loc[gene].values.astype(np.float32)
            
            if self.gpu_available:
                gene_gpu = cp.asarray(gene_expression)
            else:
                gene_gpu = gene_expression
            
            # Get top lncRNAs for this gene
            lncrna_corrs = self._get_top_correlations_gpu(
                gene_gpu, 'lncrna', top_n=50
            )
            
            # Get top miRNAs for this gene
            mirna_corrs = self._get_top_correlations_gpu(
                gene_gpu, 'mirna', top_n=25
            )
            
            # Analyze interactions for top regulators
            for lncrna_name, lncrna_corr, lncrna_pval in lncrna_corrs[:15]:
                for mirna_name, mirna_corr, mirna_pval in mirna_corrs[:10]:
                    interaction_result = self._analyze_interaction_gpu(
                        gene, gene_expression,
                        lncrna_name, mirna_name,
                        'lncrna', 'mirna',
                        'lncrna_mirna'
                    )
                    if interaction_result:
                        results.append(interaction_result)
        
        return results
    
    def gpu_analyze_multi_way_interactions(self) -> pd.DataFrame:
        """Analyze complex multi-way regulatory interactions using GPU."""
        all_results = []
        
        gene_index = self.datasets['gene'].index.tolist()
        n_genes = len(gene_index)
        
        n_batches = (n_genes + self.batch_size - 1) // self.batch_size
        print(f"  📊 Processing {n_genes} genes in {n_batches} batches...")
        
        for batch_idx in range(n_batches):
            start_idx = batch_idx * self.batch_size
            end_idx = min((batch_idx + 1) * self.batch_size, n_genes)
            batch_genes = gene_index[start_idx:end_idx]
            
            batch_results = self._process_multi_way_batch_gpu(batch_genes)
            all_results.extend(batch_results)
            
            if (batch_idx + 1) % 10 == 0 or batch_idx == n_batches - 1:
                print(f"    ✅ Batch {batch_idx + 1}/{n_batches} completed: {len(all_results)} results so far")
            
            if self.gpu_available and (batch_idx + 1) % 20 == 0:
                self._clear_gpu_memory()
        
        if all_results:
            results_df = pd.DataFrame(all_results)
            print(f"  🎯 Total multi-way interactions: {len(results_df)}")
            return results_df
        else:
            return pd.DataFrame()
    
    def _process_multi_way_batch_gpu(self, gene_batch: List[str]) -> List[Dict]:
        """Process a batch of genes for multi-way interaction analysis using GPU."""
        results = []
        
        for gene in gene_batch:
            gene_expression = self.datasets['gene'].loc[gene].values.astype(np.float32)
            
            if self.gpu_available:
                gene_gpu = cp.asarray(gene_expression)
            else:
                gene_gpu = gene_expression
            
            # Get top regulators for this gene
            regulators = self._get_top_regulators_gpu(gene, gene_gpu, 20)
            
            if len(regulators) >= 3:
                multi_way_result = self._analyze_multi_regulator_gpu(
                    gene_expression, regulators, gene
                )
                
                if multi_way_result:
                    results.append(multi_way_result)
        
        return results
    
    def _get_top_regulators_gpu(self, gene: str, gene_gpu, n_regulators: int) -> Dict[str, np.ndarray]:
        """Get top regulators for a specific gene using GPU."""
        regulators = {}
        
        # Get top miRNA regulators
        mirna_corrs = self._get_top_correlations_gpu(gene_gpu, 'mirna', top_n=15)
        for mirna_name, corr, pval in mirna_corrs:
            regulators[mirna_name] = self.datasets['mirna'].loc[mirna_name.replace('mirna_', '')].values
        
        # Get top lncRNA regulators
        lncrna_corrs = self._get_top_correlations_gpu(gene_gpu, 'lncrna', top_n=30)
        for lncrna_name, corr, pval in lncrna_corrs:
            regulators[lncrna_name] = self.datasets['lncrna'].loc[lncrna_name.replace('lncrna_', '')].values
        
        # Get top methylation regulators
        meth_corrs = self._get_top_correlations_gpu(gene_gpu, 'methylation', top_n=25)
        for meth_name, corr, pval in meth_corrs:
            regulators[meth_name] = self.datasets['methylation'].loc[meth_name.replace('methylation_', '')].values
        
        return regulators
    
    def _analyze_multi_regulator_gpu(self, gene_expression: np.ndarray, 
                                      regulators: Dict[str, np.ndarray],
                                      gene_name: str) -> Optional[Dict]:
        """Analyze multi-regulator interactions using GPU."""
        try:
            # Create feature matrix
            feature_data = {reg_name: reg_values for reg_name, reg_values in regulators.items()}
            feature_data['target'] = gene_expression
            data = pd.DataFrame(feature_data)
            
            # Scale data
            scaler = StandardScaler()
            data_scaled = pd.DataFrame(
                scaler.fit_transform(data),
                columns=data.columns
            )
            
            X = data_scaled.drop('target', axis=1).values
            y = data_scaled['target'].values
            
            if self.gpu_available:
                X_gpu = cp.asarray(X.astype(np.float32))
                y_gpu = cp.asarray(y.astype(np.float32))
                
                # Base model (first regulator only)
                _, r2_base = self.gpu_linear_regression_batch(cp.asnumpy(X_gpu[:, :1]), cp.asnumpy(y_gpu))
                
                # Full model (all regulators)
                _, r2_full = self.gpu_linear_regression_batch(cp.asnumpy(X_gpu), cp.asnumpy(y_gpu))
            else:
                model_base = LinearRegression()
                model_base.fit(X[:, :1], y)
                r2_base = model_base.score(X[:, :1], y)
                
                model_full = LinearRegression()
                model_full.fit(X, y)
                r2_full = model_full.score(X, y)
            
            improvement = r2_full - r2_base
            
            # Statistical testing
            n_samples = len(data_scaled)
            df1 = len(regulators) - 1
            df2 = n_samples - len(regulators) - 1
            
            if (df2 > 0 and improvement > 0 and r2_full < 1.0 and (1 - r2_full) > 1e-10):
                try:
                    f_stat = (improvement / df1) / ((1 - r2_full) / df2)
                    p_value = 1 - stats.f.cdf(f_stat, df1, df2)
                    has_significant = p_value < 0.05
                except:
                    has_significant = False
                    p_value = 1.0
            else:
                has_significant = False
                p_value = 1.0
            
            return {
                'gene': gene_name,
                'n_regulators': len(regulators),
                'r2_base_model': r2_base,
                'r2_full_model': r2_full,
                'improvement_from_regulators': improvement,
                'has_significant_interactions': has_significant,
                'interaction_p_value': p_value,
                'regulator_types': list(regulators.keys())
            }
            
        except Exception as e:
            return None
    
    def gpu_infer_context_specific_networks(self) -> Dict:
        """Infer context-specific regulatory networks using GPU."""
        print("  🔄 Inferring context-specific regulatory networks...")
        
        context_networks = {}
        
        contexts = [
            ('high_mirna', 0.5),
            ('low_mirna', -0.5),
            ('high_methylation', 0.5)
        ]
        
        for context_name, threshold in contexts:
            print(f"    📊 Processing {context_name} context...")
            context_result = self._analyze_context_network_gpu(context_name, threshold)
            context_networks[context_name] = context_result
            print(f"    ✅ {context_name} completed")
        
        return context_networks
    
    def _analyze_context_network_gpu(self, context_name: str, threshold: float) -> Dict:
        """Analyze regulatory network for a specific context using GPU."""
        context_networks = {
            'gene_mirna_correlations': [],
            'gene_lncrna_correlations': [],
            'gene_methylation_correlations': []
        }
        
        gene_index = self.datasets['gene'].index.tolist()
        
        # Get context mask
        if 'high_mirna' in context_name or 'low_mirna' in context_name:
            sentinel_mirna = self.datasets['mirna'].index[0]
            mirna_values = self.datasets['mirna'].loc[sentinel_mirna].values
            scaler = StandardScaler()
            z_scores = scaler.fit_transform(mirna_values.reshape(-1, 1)).ravel()
            
            if threshold > 0:
                context_mask = z_scores > threshold
            else:
                context_mask = z_scores < abs(threshold)
                
        elif 'high_methylation' in context_name:
            sentinel_cpg = self.datasets['methylation'].index[0]
            meth_values = self.datasets['methylation'].loc[sentinel_cpg].values
            scaler = StandardScaler()
            z_scores = scaler.fit_transform(meth_values.reshape(-1, 1)).ravel()
            context_mask = z_scores > threshold
        else:
            context_mask = np.ones(len(self.samples), dtype=bool)
        
        context_samples = [col for i, col in enumerate(self.datasets['gene'].columns) if context_mask[i]]
        
        if len(context_samples) < 10:
            return context_networks
        
        # Process in batches
        n_batches = (len(gene_index) + self.batch_size - 1) // self.batch_size
        
        for batch_idx in range(n_batches):
            start_idx = batch_idx * self.batch_size
            end_idx = min((batch_idx + 1) * self.batch_size, len(gene_index))
            batch_genes = gene_index[start_idx:end_idx]
            
            batch_results = self._process_context_network_batch_gpu(batch_genes, context_samples)
            
            for corr_type in context_networks:
                context_networks[corr_type].extend(batch_results[corr_type])
            
            if self.gpu_available and (batch_idx + 1) % 20 == 0:
                self._clear_gpu_memory()
        
        return context_networks
    
    def _process_context_network_batch_gpu(self, gene_batch: List[str], 
                                            context_samples: List[str]) -> Dict:
        """Process a batch of genes for context network analysis using GPU."""
        chunk_results = {
            'gene_mirna_correlations': [],
            'gene_lncrna_correlations': [],
            'gene_methylation_correlations': []
        }
        
        for gene in gene_batch:
            gene_expression = self.datasets['gene'].loc[gene, context_samples].values.astype(np.float32)
            
            if self.gpu_available:
                gene_gpu = cp.asarray(gene_expression).reshape(1, -1)
                
                # miRNA correlations
                mirna_data = self.datasets['mirna'][context_samples].values.astype(np.float32)
                mirna_gpu = cp.asarray(mirna_data)
                
                corrs, pvals = self.gpu_pearson_correlation_batch(gene_gpu, mirna_gpu)
                corrs_np = cp.asnumpy(corrs[0])
                pvals_np = cp.asnumpy(pvals[0])
                
                for i, mirna in enumerate(self.datasets['mirna'].index):
                    if pvals_np[i] < 0.1:
                        chunk_results['gene_mirna_correlations'].append({
                            'gene': gene, 'mirna': mirna, 
                            'correlation': float(corrs_np[i]), 
                            'p_value': float(pvals_np[i])
                        })
                
                # lncRNA correlations
                lncrna_data = self.datasets['lncrna'][context_samples].values.astype(np.float32)
                lncrna_gpu = cp.asarray(lncrna_data)
                
                corrs, pvals = self.gpu_pearson_correlation_batch(gene_gpu, lncrna_gpu)
                corrs_np = cp.asnumpy(corrs[0])
                pvals_np = cp.asnumpy(pvals[0])
                
                for i, lncrna in enumerate(self.datasets['lncrna'].index):
                    if pvals_np[i] < 0.1:
                        chunk_results['gene_lncrna_correlations'].append({
                            'gene': gene, 'lncrna': lncrna,
                            'correlation': float(corrs_np[i]),
                            'p_value': float(pvals_np[i])
                        })
                
                # Methylation correlations
                meth_data = self.datasets['methylation'][context_samples].values.astype(np.float32)
                meth_gpu = cp.asarray(meth_data)
                
                corrs, pvals = self.gpu_pearson_correlation_batch(gene_gpu, meth_gpu)
                corrs_np = cp.asnumpy(corrs[0])
                pvals_np = cp.asnumpy(pvals[0])
                
                for i, cpg in enumerate(self.datasets['methylation'].index):
                    if pvals_np[i] < 0.1:
                        chunk_results['gene_methylation_correlations'].append({
                            'gene': gene, 'cpg': cpg,
                            'correlation': float(corrs_np[i]),
                            'p_value': float(pvals_np[i])
                        })
            else:
                # CPU fallback
                for mirna in self.datasets['mirna'].index:
                    mirna_expression = self.datasets['mirna'].loc[mirna, context_samples].values
                    if len(mirna_expression) > 5:
                        corr, pval = pearsonr(gene_expression, mirna_expression)
                        if pval < 0.1:
                            chunk_results['gene_mirna_correlations'].append({
                                'gene': gene, 'mirna': mirna, 'correlation': corr, 'p_value': pval
                            })
                
                for lncrna in self.datasets['lncrna'].index:
                    lncrna_expression = self.datasets['lncrna'].loc[lncrna, context_samples].values
                    if len(lncrna_expression) > 5:
                        corr, pval = pearsonr(gene_expression, lncrna_expression)
                        if pval < 0.1:
                            chunk_results['gene_lncrna_correlations'].append({
                                'gene': gene, 'lncrna': lncrna, 'correlation': corr, 'p_value': pval
                            })
                
                for cpg in self.datasets['methylation'].index:
                    meth_expression = self.datasets['methylation'].loc[cpg, context_samples].values
                    if len(meth_expression) > 5:
                        corr, pval = pearsonr(gene_expression, meth_expression)
                        if pval < 0.1:
                            chunk_results['gene_methylation_correlations'].append({
                                'gene': gene, 'cpg': cpg, 'correlation': corr, 'p_value': pval
                            })
        
        return chunk_results

    # =========================================================================
    # VISUALIZATION METHODS (Same as CPU version)
    # =========================================================================
    
    def generate_context_visualizations(self):
        """Generate context-dependent analysis visualizations."""
        print("\n" + "="*60)
        print("GENERATING CONTEXT-DEPENDENT VISUALIZATIONS")
        print("="*60)
        
        self.plot_context_dependent_interactions()
        self.plot_context_networks()
        self.plot_interaction_improvements()
        
        print(f"✅ All visualizations saved to {self.plots_dir}/")
    
    def plot_context_dependent_interactions(self):
        """Plot context-dependent interaction analysis."""
        if 'context_dependent' not in self.results:
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        meth_mirna = self.results['context_dependent']['methylation_mirna_context']
        if not meth_mirna.empty:
            axes[0, 0].hist(meth_mirna['improvement_from_interaction'], bins=30, alpha=0.7, edgecolor='black')
            axes[0, 0].set_title('Methylation-miRNA Context Interactions')
            axes[0, 0].set_xlabel('Improvement from Interaction')
            axes[0, 0].set_ylabel('Frequency')
            axes[0, 0].axvline(x=0.1, color='red', linestyle='--', alpha=0.7, label='Significance threshold')
            axes[0, 0].legend()
        
        lncrna_mirna = self.results['context_dependent']['lncrna_mirna_context']
        if not lncrna_mirna.empty:
            axes[0, 1].hist(lncrna_mirna['improvement_from_interaction'], bins=30, alpha=0.7, edgecolor='black')
            axes[0, 1].set_title('lncRNA-miRNA Context Interactions')
            axes[0, 1].set_xlabel('Improvement from Interaction')
            axes[0, 1].set_ylabel('Frequency')
            axes[0, 1].axvline(x=0.1, color='red', linestyle='--', alpha=0.7, label='Significance threshold')
            axes[0, 1].legend()
        
        if not meth_mirna.empty:
            axes[1, 0].hist(meth_mirna['context_strength'].dropna(), bins=30, alpha=0.7, edgecolor='black')
            axes[1, 0].set_title('Methylation-miRNA Context Strength')
            axes[1, 0].set_xlabel('Context Strength')
            axes[1, 0].set_ylabel('Frequency')
        
        if not lncrna_mirna.empty:
            axes[1, 1].hist(lncrna_mirna['context_strength'].dropna(), bins=30, alpha=0.7, edgecolor='black')
            axes[1, 1].set_title('lncRNA-miRNA Context Strength')
            axes[1, 1].set_xlabel('Context Strength')
            axes[1, 1].set_ylabel('Frequency')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.plots_dir, "context_dependent_interactions.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_context_networks(self):
        """Plot context-specific regulatory networks."""
        if 'context_dependent' not in self.results:
            return
        
        context_networks = self.results['context_dependent']['context_networks']
        if not context_networks:
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        contexts = list(context_networks.keys())
        for i, context in enumerate(contexts[:4]):
            if context in context_networks:
                network = context_networks[context]
                
                mirna_count = len(network.get('gene_mirna_correlations', []))
                lncrna_count = len(network.get('gene_lncrna_correlations', []))
                meth_count = len(network.get('gene_methylation_correlations', []))
                
                categories = ['miRNA', 'lncRNA', 'Methylation']
                counts = [mirna_count, lncrna_count, meth_count]
                
                axes[i//2, i%2].bar(categories, counts, alpha=0.7, color=['blue', 'green', 'red'])
                axes[i//2, i%2].set_title(f'{context.replace("_", " ").title()} Network')
                axes[i//2, i%2].set_ylabel('Number of Regulatory Relationships')
                
                for j, count in enumerate(counts):
                    axes[i//2, i%2].text(j, count + 1, str(count), ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.plots_dir, "context_networks.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_interaction_improvements(self):
        """Plot interaction improvement distributions."""
        if 'context_dependent' not in self.results:
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        meth_mirna = self.results['context_dependent']['methylation_mirna_context']
        lncrna_mirna = self.results['context_dependent']['lncrna_mirna_context']
        
        if not meth_mirna.empty and not lncrna_mirna.empty:
            data_to_plot = [
                meth_mirna['improvement_from_interaction'].dropna(),
                lncrna_mirna['improvement_from_interaction'].dropna()
            ]
            
            axes[0].boxplot(data_to_plot, labels=['Methylation-miRNA', 'lncRNA-miRNA'])
            axes[0].set_title('Interaction Improvement Comparison')
            axes[0].set_ylabel('Improvement from Interaction')
            axes[0].grid(True, alpha=0.3)
            
            data_to_plot = [
                meth_mirna['context_strength'].dropna(),
                lncrna_mirna['context_strength'].dropna()
            ]
            
            axes[1].boxplot(data_to_plot, labels=['Methylation-miRNA', 'lncRNA-miRNA'])
            axes[1].set_title('Context Strength Comparison')
            axes[1].set_ylabel('Context Strength')
            axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.plots_dir, "interaction_improvements.png"), dpi=300, bbox_inches='tight')
        plt.close()

    # =========================================================================
    # RESULTS SAVING AND REPORTING
    # =========================================================================
    
    def save_context_results(self):
        """Save all context-dependent analysis results."""
        print("\n" + "="*60)
        print("SAVING CONTEXT-DEPENDENT RESULTS")
        print("="*60)
        
        if 'context_dependent' in self.results:
            for analysis_type, results in self.results['context_dependent'].items():
                if isinstance(results, pd.DataFrame) and not results.empty:
                    results.to_csv(os.path.join(self.tables_dir, f"{analysis_type}.csv"), index=False)
                    print(f"Saved {analysis_type}: {len(results)} results")
                elif isinstance(results, dict):
                    for context_name, context_data in results.items():
                        if isinstance(context_data, dict):
                            for corr_type, corr_data in context_data.items():
                                if isinstance(corr_data, list) and corr_data:
                                    corr_df = pd.DataFrame(corr_data)
                                    corr_df.to_csv(os.path.join(self.tables_dir, f"{context_name}_{corr_type}.csv"), index=False)
                                    print(f"Saved {context_name}_{corr_type}: {len(corr_data)} correlations")
        
        print(f"✅ All results saved to {self.tables_dir}/")
    
    def generate_markdown_report(self):
        """Generate comprehensive markdown report."""
        print("\n" + "="*60)
        print("GENERATING MARKDOWN REPORT")
        print("="*60)
        
        report_path = os.path.join(self.reports_dir, "context_dependent_analysis_report.md")
        
        with open(report_path, 'w') as f:
            f.write("# GPU-Optimized Context-Dependent Regulation Analysis Report\n\n")
            f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Analysis ID:** {self.timestamp}\n")
            f.write(f"**Compute Mode:** {'GPU (RAPIDS)' if self.gpu_available else 'CPU (Fallback)'}\n\n")
            
            f.write("## Executive Summary\n\n")
            f.write("This report presents the results of GPU-accelerated context-dependent regulatory interaction analysis.\n\n")
            
            f.write("## Analysis Configuration\n\n")
            if self.gpu_available:
                f.write(f"- **GPU Device:** {getattr(self, 'gpu_name', 'Unknown')}\n")
                f.write(f"- **GPU Memory:** {getattr(self, 'gpu_memory_total', 'N/A'):.1f} GB\n")
                f.write(f"- **Batch Size:** {self.batch_size} genes\n")
            else:
                f.write(f"- **CPU Workers:** {self.n_jobs}\n")
                f.write(f"- **Available RAM:** {self._get_available_ram():.1f} GB\n")
            f.write("\n")
            
            f.write("## Datasets\n\n")
            for data_type, data in self.datasets.items():
                f.write(f"- **{data_type.capitalize()}:** {data.shape[1]} samples × {data.shape[0]} features\n")
            f.write("\n")
            
            if 'context_dependent' in self.results:
                f.write("## Results Summary\n\n")
                
                meth_mirna = self.results['context_dependent'].get('methylation_mirna_context', pd.DataFrame())
                if not meth_mirna.empty:
                    f.write("### Methylation-miRNA Context Analysis\n\n")
                    f.write(f"- **Total interactions:** {len(meth_mirna)}\n")
                    f.write(f"- **Context-dependent:** {meth_mirna['context_dependent'].sum()}\n")
                    f.write(f"- **Mean improvement:** {meth_mirna['improvement_from_interaction'].mean():.3f}\n\n")
                
                lncrna_mirna = self.results['context_dependent'].get('lncrna_mirna_context', pd.DataFrame())
                if not lncrna_mirna.empty:
                    f.write("### lncRNA-miRNA Context Analysis\n\n")
                    f.write(f"- **Total interactions:** {len(lncrna_mirna)}\n")
                    f.write(f"- **Context-dependent:** {lncrna_mirna['context_dependent'].sum()}\n")
                    f.write(f"- **Mean improvement:** {lncrna_mirna['improvement_from_interaction'].mean():.3f}\n\n")
                
                multi_way = self.results['context_dependent'].get('multi_way_interactions', pd.DataFrame())
                if not multi_way.empty:
                    f.write("### Multi-Way Interactions\n\n")
                    f.write(f"- **Total genes:** {len(multi_way)}\n")
                    f.write(f"- **Significant interactions:** {multi_way['has_significant_interactions'].sum()}\n\n")
            
            f.write("---\n")
            f.write(f"*Generated by GPUContextDependentRegulationAnalysis*\n")
        
        print(f"✅ Report saved: {report_path}")
    
    def print_context_summary(self):
        """Print summary of context-dependent analysis."""
        print("\n" + "="*60)
        print("CONTEXT-DEPENDENT ANALYSIS SUMMARY")
        print("="*60)
        
        if 'context_dependent' not in self.results:
            print("No context-dependent results available.")
            return
        
        meth_mirna = self.results['context_dependent']['methylation_mirna_context']
        if not meth_mirna.empty:
            print(f"\nMETHYLATION-MIRNA CONTEXT ANALYSIS:")
            print(f"  Total interactions analyzed: {len(meth_mirna)}")
            print(f"  Context-dependent interactions: {meth_mirna['context_dependent'].sum()}")
            print(f"  Mean improvement: {meth_mirna['improvement_from_interaction'].mean():.3f}")
        
        lncrna_mirna = self.results['context_dependent']['lncrna_mirna_context']
        if not lncrna_mirna.empty:
            print(f"\nLNCRNA-MIRNA CONTEXT ANALYSIS:")
            print(f"  Total interactions analyzed: {len(lncrna_mirna)}")
            print(f"  Context-dependent interactions: {lncrna_mirna['context_dependent'].sum()}")
            print(f"  Mean improvement: {lncrna_mirna['improvement_from_interaction'].mean():.3f}")
        
        multi_way = self.results['context_dependent']['multi_way_interactions']
        if not multi_way.empty:
            print(f"\nMULTI-WAY INTERACTION ANALYSIS:")
            print(f"  Total genes analyzed: {len(multi_way)}")
            print(f"  Significant interactions: {multi_way['has_significant_interactions'].sum()}")
        
        context_networks = self.results['context_dependent']['context_networks']
        if context_networks:
            print(f"\nCONTEXT-SPECIFIC NETWORKS:")
            for context_name, network in context_networks.items():
                total = (
                    len(network.get('gene_mirna_correlations', [])) +
                    len(network.get('gene_lncrna_correlations', [])) +
                    len(network.get('gene_methylation_correlations', []))
                )
                print(f"  {context_name}: {total} regulatory relationships")
    
    def run_complete_context_analysis(self):
        """Run the complete GPU-optimized context-dependent analysis."""
        start_time = time.time()
        
        print("="*80)
        print(f"{'🚀 GPU' if self.gpu_available else '💻 CPU'} CONTEXT-DEPENDENT REGULATION ANALYSIS")
        print("="*80)
        if self.gpu_available:
            print(f"  • GPU: {getattr(self, 'gpu_name', 'Unknown')}")
            print(f"  • Memory: {getattr(self, 'gpu_memory_available', 'N/A'):.1f} GB available")
            print(f"  • Batch size: {self.batch_size} genes per batch")
        else:
            print(f"  • CPU workers: {self.n_jobs}")
            print(f"  • RAM: {self._get_available_ram():.1f} GB")
        print("="*80)
        
        # Run analysis
        self.analyze_context_dependent_regulation()
        
        # Generate outputs
        self.generate_context_visualizations()
        self.save_context_results()
        self.generate_markdown_report()
        self.print_context_summary()
        
        total_time = time.time() - start_time
        print(f"\n⏱️  Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
        print("\n" + "="*80)
        print("🎉 ANALYSIS COMPLETED SUCCESSFULLY!")
        print("="*80)


def main():
    """Main function to run the GPU-optimized context-dependent analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(description='GPU-optimized Context-Dependent Regulation Analysis')
    parser.add_argument('--data-dir', type=str, default='data/cleaned_datasets',
                        help='Directory containing cleaned datasets')
    parser.add_argument('--batch-size', type=int, default=1000,
                        help='Number of genes per GPU batch')
    parser.add_argument('--gpu-memory', type=float, default=0.8,
                        help='Fraction of GPU memory to use (0.0-1.0)')
    
    args = parser.parse_args()
    
    # Initialize and run analysis
    analysis = GPUContextDependentRegulationAnalysis(
        data_dir=args.data_dir,
        batch_size=args.batch_size,
        gpu_memory_fraction=args.gpu_memory
    )
    
    analysis.run_complete_context_analysis()


if __name__ == "__main__":
    main()

