# Critical Bug Fixes Applied to Context-Dependent Analysis

## Overview
This document summarizes the critical bugs that were identified and fixed in the context-dependent analysis Python scripts to ensure proper statistical analysis and correct sample selection.

## Files Fixed
- ✅ **`code/context_dependent_analysis.py`** - Main analysis script with all critical bug fixes
- ✅ **`code/subset_context_dependent_analysis.py`** - Subset analysis script with all critical bug fixes

## Bug Fixes Applied

### 1. ✅ Context Mask Uses Features, Not Samples (CRITICAL)
**Problem**: In `_analyze_context_network_parallel`, context was computed from `mirna_means` (indexed by miRNA features) and then applied to samples via `i % len(context_mask)`. This mismatched axes and yielded quasi-random sample subsets.

**Fix Applied**: 
- Implemented proper per-sample context variables using sentinel approach
- For miRNA contexts: Use first miRNA as context proxy, compute z-scores per sample
- For methylation contexts: Use first CpG as context proxy, compute z-scores per sample
- Context mask now properly indexes samples, not features

**Code Changes**:
```python
# FIXED: Get context mask based on per-sample context variables, not features
if 'high_mirna' in context_name or 'low_mirna' in context_name:
    # Sentinel miRNA approach: use first miRNA as context proxy
    sentinel_mirna = self.datasets['mirna'].index[0]
    mirna_values = self.datasets['mirna'].loc[sentinel_mirna].values
    z_scores = StandardScaler().fit_transform(mirna_values.reshape(-1, 1)).ravel()
    
    if threshold > 0:  # high_mirna
        context_mask = z_scores > threshold
    else:  # low_mirna
        context_mask = z_scores < abs(threshold)
        
elif 'high_methylation' in context_name:
    # FIXED: Implement methylation-based context
    sentinel_cpg = self.datasets['methylation'].index[0]
    meth_values = self.datasets['methylation'].loc[sentinel_cpg].values
    z_scores = StandardScaler().fit_transform(meth_values.reshape(-1, 1)).ravel()
    context_mask = z_scores > threshold
```

### 2. ✅ "High Methylation" Context Not Implemented (CRITICAL)
**Problem**: `contexts_to_analyze` included 'high_methylation' but mask was derived only from miRNAs.

**Fix Applied**: 
- Implemented methylation-based context using sentinel CpG approach
- Uses first CpG site as context proxy, computes z-scores per sample
- Applies threshold to identify high methylation samples

### 3. ✅ In-Sample R² Used to Declare Significance (CRITICAL)
**Problem**: Using training R² improvements (threshold 0.1) inflates false positives.

**Fix Applied**: 
- Replaced arbitrary R² threshold with proper nested F-tests
- Compares models with vs without interaction terms using F-statistics
- Uses p-value threshold (0.05) for statistical significance
- Applied to both methylation-miRNA and lncRNA-miRNA interactions

**Code Changes**:
```python
# FIXED: Use proper statistical testing instead of arbitrary R^2 threshold
# Perform nested F-test to compare models with and without interaction
from sklearn.metrics import r2_score
from scipy import stats

# Calculate F-statistic for nested model comparison
n_samples = len(data_scaled)
df1 = 1  # Additional parameter in full model
df2 = n_samples - 4  # Degrees of freedom for full model (4 parameters)

if df2 > 0 and improvement_from_interaction > 0:
    f_stat = (improvement_from_interaction / df1) / ((1 - r2_3) / df2)
    p_value = 1 - stats.f.cdf(f_stat, df1, df2)
    context_dependent = p_value < 0.05  # Use p-value threshold
else:
    context_dependent = False
    p_value = 1.0
```

### 4. ✅ NaN Handling in Context Direction (CRITICAL)
**Problem**: 'positive'/'negative' was assigned even when `corr_high` or `corr_low` was NaN.

**Fix Applied**: 
- Check for NaN values before assigning context direction
- If any correlation is NaN, set `context_direction='NA'` and `context_strength=NaN`
- Applied to both interaction analysis methods

**Code Changes**:
```python
# FIXED: Handle NaN values in context direction
if pd.isna(corr_high_mirna) or pd.isna(corr_low_mirna):
    context_direction = 'NA'
    context_strength = np.nan
else:
    context_direction = 'positive' if corr_high_mirna > corr_low_mirna else 'negative'
```

### 5. ✅ Unused/Incorrect Precomputations and Reporting (CRITICAL)
**Problem**: Huge sample correlation matrices were computed but unused; MD report labels "samples × features" reversed.

**Fix Applied**: 
- Removed unused correlation matrix precomputations to save memory
- Fixed dimension reporting in markdown report: `{data.shape[1]} samples × {data.shape[0]} features`

**Code Changes**:
```python
# FIXED: Remove unused huge correlation matrices to save memory
# These were computed but never used in the analysis

# Fixed dimension reporting
for data_type, data in self.datasets.items():
    f.write(f"  - {data_type.capitalize()}: {data.shape[1]} samples × {data.shape[0]} features\n")
```

### 6. ✅ Division by Zero Errors in F-Test Calculations (CRITICAL)
**Problem**: "float division by zero" errors occurred when calculating F-statistics for nested model comparisons, particularly when R² values were 1.0 or when degrees of freedom were insufficient.

**Fix Applied**: 
- Added comprehensive error handling for edge cases in F-test calculations
- Prevent division by zero when R² = 1.0 or when (1 - R²) is very small
- Added try-catch blocks to handle numerical errors gracefully
- Applied to all interaction analysis methods (methylation-miRNA, lncRNA-miRNA, multi-way)

**Code Changes**:
```python
# FIXED: Add proper error handling for edge cases
if (df2 > 0 and 
    improvement_from_interaction > 0 and 
    r2_3 < 1.0 and  # Prevent division by zero when R² = 1
    (1 - r2_3) > 1e-10):  # Prevent division by very small numbers
    try:
        f_stat = (improvement_from_interaction / df1) / ((1 - r2_3) / df2)
        p_value = 1 - stats.f.cdf(f_stat, df1, df2)
        context_dependent = p_value < 0.05  # Use p-value threshold
    except (ZeroDivisionError, ValueError, RuntimeWarning):
        # Handle any numerical errors gracefully
        context_dependent = False
        p_value = 1.0
else:
    context_dependent = False
    p_value = 1.0
```

### 7. ✅ Matplotlib Backend Errors in Headless Environments (CRITICAL)
**Problem**: "matplotlib.backend_bases._get_renderer.<locals>.Done" errors occurred when trying to save plots in headless environments (servers without display), causing the visualization step to crash.

**Fix Applied**: 
- Configured matplotlib to use 'Agg' backend for headless environments
- Added comprehensive error handling to all plotting functions
- Implemented graceful degradation when plotting fails
- Added proper cleanup of matplotlib figures to prevent memory leaks

**Code Changes**:
```python
# FIXED: Configure matplotlib for headless environments
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless environments

# FIXED: Add error handling for plot saving
try:
    plot_path = os.path.join(self.plots_dir, "context_dependent_interactions.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"    ✅ Saved context-dependent interactions plot: {plot_path}")
except Exception as e:
    print(f"    ❌ Error saving plot: {e}")
finally:
    plt.close()
```

## Additional Improvements

### Statistical Rigor
- Added `interaction_p_value` field to all interaction results
- Implemented proper F-test calculations for nested model comparisons
- Applied consistent p-value threshold (0.05) across all analyses
- Added robust error handling for numerical edge cases

### Visualization Stability
- Configured matplotlib for headless server environments
- Added comprehensive error handling to all plotting functions
- Implemented graceful degradation when plotting fails
- Added proper cleanup to prevent memory leaks

### Memory Efficiency
- Removed unused correlation matrix precomputations
- Maintained numpy array conversions for vectorized operations
- Optimized context mask generation

### Code Quality
- Added comprehensive documentation of bug fixes
- Ensured consistent error handling
- Maintained backward compatibility with existing result structures
- Added graceful degradation for problematic statistical calculations

## Impact of Fixes

These fixes ensure:
1. **Correct sample selection** for context-dependent analysis
2. **Proper statistical significance testing** using F-tests instead of arbitrary thresholds
3. **Robust handling of edge cases** (NaN values, insufficient samples, division by zero)
4. **Memory efficiency** by removing unused computations
5. **Accurate reporting** of dataset dimensions and analysis results
6. **Stable execution** without crashes from numerical errors
7. **Reliable visualization** in headless server environments

## Testing Recommendations

After applying these fixes, verify:
1. Context masks properly select samples (not features)
2. High methylation context works correctly
3. Statistical significance is properly calculated
4. NaN values are handled gracefully
5. Memory usage is reasonable
6. Report dimensions are correct
7. No "division by zero" errors occur during analysis
8. Analysis completes successfully for all gene subsets
9. Visualization step completes without matplotlib backend errors
10. Plots are generated successfully in headless environments

## Files Modified
- `code/context_dependent_analysis.py` - Main analysis code with all critical bug fixes
- `code/subset_context_dependent_analysis.py` - Subset analysis code with all critical bug fixes
- `code/CRITICAL_BUG_FIXES_SUMMARY.md` - This summary document

All critical bugs have been addressed in both Python scripts, and the code should now perform statistically sound context-dependent regulatory analysis with proper sample selection, significance testing, and robust error handling.
