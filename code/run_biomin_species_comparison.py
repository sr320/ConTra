#!/usr/bin/env python3
"""
Cross-Species Biomin Analysis and Comparison

This script:
1. Runs context_dependent_analysis.py on all three species in full-species-biomin
2. Extracts OG_XXXXX orthologous group IDs from results
3. Compares context-dependent interactions across species to identify:
   - Conserved regulatory patterns (same OG in multiple species)
   - Species-specific regulatory patterns
   - Cross-species correlation of effect sizes

Usage:
    python code/run_biomin_species_comparison.py [--skip-analysis] [--output-dir PATH]
    
Options:
    --skip-analysis: Skip running the per-species analysis (use existing output)
    --output-dir: Directory for comparison results (default: output/biomin_comparison_TIMESTAMP)
"""

import os
import sys
import subprocess
import argparse
import re
from datetime import datetime
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

# Species configuration
SPECIES_CONFIG = {
    'apul': {
        'name': 'A. pulchra',
        'data_dir': 'data/full-species-biomin/cleaned_apul',
        'color': '#E74C3C'  # Red
    },
    'peve': {
        'name': 'P. evermanni', 
        'data_dir': 'data/full-species-biomin/cleaned_peve',
        'color': '#3498DB'  # Blue
    },
    'ptua': {
        'name': 'P. tuahiniensis',
        'data_dir': 'data/full-species-biomin/cleaned_ptua', 
        'color': '#2ECC71'  # Green
    }
}


def get_workspace_root() -> str:
    """Get workspace root directory."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def extract_og_id(gene_id: str) -> Optional[str]:
    """Extract OG_XXXXX ID from gene identifier.
    
    Handles formats like:
    - "('OG_13910', 'FUN_002435')"
    - "OG_13910"
    - "OG_13910_FUN_002435"
    """
    if pd.isna(gene_id):
        return None
    
    gene_id = str(gene_id)
    
    # Pattern to match OG_XXXXX (digits after OG_)
    match = re.search(r'OG_(\d+)', gene_id)
    if match:
        return f"OG_{match.group(1)}"
    return None


def run_species_analysis(species_key: str, workspace_root: str) -> Optional[str]:
    """Run context_dependent_analysis.py for a single species.
    
    Returns the output directory path if successful, None otherwise.
    """
    config = SPECIES_CONFIG[species_key]
    data_dir = os.path.join(workspace_root, config['data_dir'])
    
    if not os.path.exists(data_dir):
        print(f"  ❌ Data directory not found: {data_dir}")
        return None
    
    print(f"\n{'='*60}")
    print(f"Running analysis for {config['name']} ({species_key})")
    print(f"{'='*60}")
    
    # Run the analysis script
    script_path = os.path.join(workspace_root, 'code', 'context_dependent_analysis.py')
    cmd = [
        sys.executable, script_path,
        '--data-dir', data_dir,
        '--n-jobs', '8'  # Use reasonable parallelism
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=workspace_root)
        
        if result.returncode != 0:
            print(f"  ❌ Analysis failed for {species_key}")
            print(f"  stderr: {result.stderr[:500]}")
            return None
        
        # Parse output directory from stdout
        for line in result.stdout.split('\n'):
            if 'Output directory:' in line:
                output_dir = line.split('Output directory:')[-1].strip()
                print(f"  ✅ Analysis complete: {output_dir}")
                return output_dir
        
        # Fallback: find most recent output directory
        output_base = os.path.join(workspace_root, 'output')
        dirs = [d for d in os.listdir(output_base) 
                if d.startswith('context_dependent_analysis_')]
        if dirs:
            latest = sorted(dirs)[-1]
            output_dir = os.path.join(output_base, latest)
            print(f"  ✅ Analysis complete: {output_dir}")
            return output_dir
            
    except Exception as e:
        print(f"  ❌ Error running analysis for {species_key}: {e}")
        return None
    
    return None


def find_latest_output_dir(species_key: str, workspace_root: str) -> Optional[str]:
    """Find the most recent output directory for a species.
    
    Looks for output directories and tries to match by data patterns.
    """
    output_base = os.path.join(workspace_root, 'output')
    
    if not os.path.exists(output_base):
        return None
    
    # Get all context_dependent_analysis directories
    dirs = [d for d in os.listdir(output_base) 
            if d.startswith('context_dependent_analysis_') and 
            os.path.isdir(os.path.join(output_base, d))]
    
    if not dirs:
        return None
    
    # Sort by timestamp (newest first)
    dirs = sorted(dirs, reverse=True)
    
    # Return the latest one for now
    # In a production setting, you might want to track which dir corresponds to which species
    return os.path.join(output_base, dirs[0])


def load_species_results(output_dir: str, species_key: str) -> Dict[str, pd.DataFrame]:
    """Load analysis results for a species."""
    results = {}
    
    tables_dir = os.path.join(output_dir, 'tables')
    
    if not os.path.exists(tables_dir):
        print(f"  ⚠️ Tables directory not found: {tables_dir}")
        return results
    
    # Load main context results
    files_to_load = [
        'methylation_mirna_context.csv',
        'lncrna_mirna_context.csv',
        'multi_way_interactions.csv'
    ]
    
    for filename in files_to_load:
        filepath = os.path.join(tables_dir, filename)
        if os.path.exists(filepath):
            try:
                df = pd.read_csv(filepath)
                # Extract OG IDs from target column
                if 'target' in df.columns:
                    df['og_id'] = df['target'].apply(extract_og_id)
                elif 'gene' in df.columns:
                    df['og_id'] = df['gene'].apply(extract_og_id)
                df['species'] = species_key
                results[filename.replace('.csv', '')] = df
                print(f"  ✅ Loaded {filename}: {len(df)} rows")
            except Exception as e:
                print(f"  ⚠️ Error loading {filename}: {e}")
    
    return results


def compare_methylation_mirna_context(
    species_results: Dict[str, Dict[str, pd.DataFrame]],
    output_dir: str
) -> pd.DataFrame:
    """Compare methylation-miRNA context results across species.
    
    Identifies OGs with context-dependent regulation in multiple species.
    """
    print("\n" + "="*60)
    print("COMPARING METHYLATION-MIRNA CONTEXT ACROSS SPECIES")
    print("="*60)
    
    # Collect all results
    all_results = []
    og_by_species: Dict[str, Dict[str, pd.DataFrame]] = defaultdict(dict)
    
    for species_key, results in species_results.items():
        if 'methylation_mirna_context' in results:
            df = results['methylation_mirna_context']
            all_results.append(df)
            
            # Group by OG ID for this species
            for og_id, group in df.groupby('og_id'):
                if og_id and pd.notna(og_id):
                    # Filter to rows with valid context_strength
                    valid_group = group.dropna(subset=['context_strength'])
                    if not valid_group.empty:
                        # Take the interaction with highest context_strength for this OG
                        best_idx = valid_group['context_strength'].idxmax()
                        og_by_species[og_id][species_key] = valid_group.loc[best_idx]
    
    if not all_results:
        print("  ⚠️ No methylation-miRNA context results found")
        return pd.DataFrame()
    
    # Find OGs present in multiple species
    multi_species_ogs = {og: species for og, species in og_by_species.items() 
                         if len(species) >= 2}
    
    print(f"\n📊 Summary:")
    print(f"  Total unique OGs across all species: {len(og_by_species)}")
    print(f"  OGs in 2+ species: {len([og for og, sp in og_by_species.items() if len(sp) >= 2])}")
    print(f"  OGs in all 3 species: {len([og for og, sp in og_by_species.items() if len(sp) == 3])}")
    
    # Create comparison dataframe
    comparison_rows = []
    
    for og_id, species_data in og_by_species.items():
        row = {'og_id': og_id, 'n_species': len(species_data)}
        
        for species_key in SPECIES_CONFIG.keys():
            if species_key in species_data:
                data = species_data[species_key]
                row[f'{species_key}_context_strength'] = data.get('context_strength', np.nan)
                row[f'{species_key}_context_dependent'] = data.get('context_dependent', False)
                row[f'{species_key}_improvement'] = data.get('improvement_from_interaction', np.nan)
                row[f'{species_key}_r2'] = data.get('r2_with_interaction', np.nan)
            else:
                row[f'{species_key}_context_strength'] = np.nan
                row[f'{species_key}_context_dependent'] = np.nan
                row[f'{species_key}_improvement'] = np.nan
                row[f'{species_key}_r2'] = np.nan
        
        comparison_rows.append(row)
    
    comparison_df = pd.DataFrame(comparison_rows)
    
    # Sort by number of species and average context strength
    if not comparison_df.empty:
        strength_cols = [f'{sp}_context_strength' for sp in SPECIES_CONFIG.keys()]
        comparison_df['mean_context_strength'] = comparison_df[strength_cols].mean(axis=1)
        comparison_df = comparison_df.sort_values(
            ['n_species', 'mean_context_strength'], 
            ascending=[False, False]
        )
        
        # Save comparison table
        comparison_path = os.path.join(output_dir, 'cross_species_methylation_mirna_comparison.csv')
        comparison_df.to_csv(comparison_path, index=False)
        print(f"\n✅ Saved comparison table: {comparison_path}")
        
        # Print top conserved OGs
        print("\n🏆 Top 15 Conserved Context-Dependent OGs (in 2+ species):")
        top_conserved = comparison_df[comparison_df['n_species'] >= 2].head(15)
        print(top_conserved[['og_id', 'n_species', 'mean_context_strength'] + 
                           [f'{sp}_context_strength' for sp in SPECIES_CONFIG.keys()]].to_string(index=False))
    
    return comparison_df


def plot_species_comparison(
    comparison_df: pd.DataFrame,
    species_results: Dict[str, Dict[str, pd.DataFrame]],
    output_dir: str
) -> None:
    """Generate comparison visualizations."""
    
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # 1. Venn-style bar chart of OG overlaps
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # OG counts by species
    species_keys = list(SPECIES_CONFIG.keys())
    og_counts = {sp: set() for sp in species_keys}
    
    for species_key, results in species_results.items():
        if 'methylation_mirna_context' in results:
            df = results['methylation_mirna_context']
            og_counts[species_key] = set(df['og_id'].dropna().unique())
    
    # Bar chart of total OGs per species
    counts = [len(og_counts[sp]) for sp in species_keys]
    colors = [SPECIES_CONFIG[sp]['color'] for sp in species_keys]
    names = [SPECIES_CONFIG[sp]['name'] for sp in species_keys]
    
    axes[0, 0].bar(names, counts, color=colors, alpha=0.7, edgecolor='black')
    axes[0, 0].set_ylabel('Number of OGs with Context-Dependent Regulation')
    axes[0, 0].set_title('OGs per Species')
    for i, (name, count) in enumerate(zip(names, counts)):
        axes[0, 0].text(i, count + 5, str(count), ha='center', fontweight='bold')
    
    # Overlap statistics
    all_ogs = set.union(*og_counts.values()) if og_counts.values() else set()
    ogs_in_2plus = set()
    ogs_in_3 = set()
    
    for og in all_ogs:
        n_species = sum(1 for sp in species_keys if og in og_counts[sp])
        if n_species >= 2:
            ogs_in_2plus.add(og)
        if n_species == 3:
            ogs_in_3.add(og)
    
    overlap_labels = ['In 1 species only', 'In 2+ species', 'In all 3 species']
    overlap_counts = [
        len(all_ogs) - len(ogs_in_2plus),
        len(ogs_in_2plus) - len(ogs_in_3),
        len(ogs_in_3)
    ]
    
    axes[0, 1].bar(overlap_labels, overlap_counts, color=['#95A5A6', '#F39C12', '#27AE60'], 
                   alpha=0.7, edgecolor='black')
    axes[0, 1].set_ylabel('Number of OGs')
    axes[0, 1].set_title('OG Conservation Across Species')
    for i, count in enumerate(overlap_counts):
        axes[0, 1].text(i, count + 2, str(count), ha='center', fontweight='bold')
    
    # 2. Scatter plot comparing context strength between species pairs
    if not comparison_df.empty and len(species_keys) >= 2:
        # Compare first two species
        sp1, sp2 = species_keys[0], species_keys[1]
        x = comparison_df[f'{sp1}_context_strength'].dropna()
        y = comparison_df[f'{sp2}_context_strength'].dropna()
        
        common_idx = x.index.intersection(y.index)
        if len(common_idx) > 5:
            x_common = x.loc[common_idx]
            y_common = y.loc[common_idx]
            
            axes[1, 0].scatter(x_common, y_common, alpha=0.5, s=20, c='#3498DB')
            axes[1, 0].set_xlabel(f'{SPECIES_CONFIG[sp1]["name"]} Context Strength')
            axes[1, 0].set_ylabel(f'{SPECIES_CONFIG[sp2]["name"]} Context Strength')
            axes[1, 0].set_title(f'Context Strength: {SPECIES_CONFIG[sp1]["name"]} vs {SPECIES_CONFIG[sp2]["name"]}')
            
            # Add correlation
            corr, pval = spearmanr(x_common, y_common)
            axes[1, 0].text(0.05, 0.95, f'ρ = {corr:.3f}\np = {pval:.2e}',
                           transform=axes[1, 0].transAxes, fontsize=10,
                           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat'))
            
            # Add diagonal line
            lims = [min(x_common.min(), y_common.min()), max(x_common.max(), y_common.max())]
            axes[1, 0].plot(lims, lims, 'r--', alpha=0.5, label='y=x')
        else:
            axes[1, 0].text(0.5, 0.5, 'Insufficient overlapping data', 
                           ha='center', va='center', transform=axes[1, 0].transAxes)
    
    # 3. Distribution of context strength by species
    strength_data = []
    for species_key, results in species_results.items():
        if 'methylation_mirna_context' in results:
            df = results['methylation_mirna_context']
            if 'context_strength' in df.columns:
                for val in df['context_strength'].dropna():
                    strength_data.append({
                        'species': SPECIES_CONFIG[species_key]['name'],
                        'context_strength': val
                    })
    
    if strength_data:
        strength_df = pd.DataFrame(strength_data)
        species_order = [SPECIES_CONFIG[sp]['name'] for sp in species_keys]
        palette = {SPECIES_CONFIG[sp]['name']: SPECIES_CONFIG[sp]['color'] for sp in species_keys}
        
        sns.violinplot(data=strength_df, x='species', y='context_strength', 
                      order=species_order, palette=palette, ax=axes[1, 1])
        axes[1, 1].set_xlabel('Species')
        axes[1, 1].set_ylabel('Context Strength')
        axes[1, 1].set_title('Distribution of Context Strength by Species')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'species_comparison_overview.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved comparison plots: {plots_dir}/species_comparison_overview.png")


def generate_comparison_report(
    comparison_df: pd.DataFrame,
    species_results: Dict[str, Dict[str, pd.DataFrame]],
    output_dir: str
) -> None:
    """Generate markdown report summarizing cross-species comparison."""
    
    report_path = os.path.join(output_dir, 'cross_species_comparison_report.md')
    
    with open(report_path, 'w') as f:
        f.write("# Cross-Species Biomin Context-Dependent Analysis Comparison\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Overview\n\n")
        f.write("This report compares context-dependent regulatory interactions across three coral species:\n\n")
        for species_key, config in SPECIES_CONFIG.items():
            f.write(f"- **{config['name']}** (`{species_key}`)\n")
        f.write("\n")
        
        f.write("## Per-Species Summary\n\n")
        for species_key, results in species_results.items():
            config = SPECIES_CONFIG[species_key]
            f.write(f"### {config['name']} ({species_key})\n\n")
            
            if 'methylation_mirna_context' in results:
                df = results['methylation_mirna_context']
                f.write(f"- Total methylation-miRNA interactions: {len(df)}\n")
                if 'context_dependent' in df.columns:
                    f.write(f"- Context-dependent interactions: {df['context_dependent'].sum()}\n")
                if 'context_strength' in df.columns:
                    f.write(f"- Mean context strength: {df['context_strength'].mean():.3f}\n")
                f.write(f"- Unique OGs: {df['og_id'].nunique()}\n")
            else:
                f.write("- No methylation-miRNA context results\n")
            f.write("\n")
        
        f.write("## Cross-Species Conservation\n\n")
        
        if not comparison_df.empty:
            n_total = len(comparison_df)
            n_2plus = len(comparison_df[comparison_df['n_species'] >= 2])
            n_all3 = len(comparison_df[comparison_df['n_species'] == 3])
            
            f.write(f"- **Total unique OGs:** {n_total}\n")
            f.write(f"- **OGs in 2+ species:** {n_2plus} ({100*n_2plus/n_total:.1f}%)\n")
            f.write(f"- **OGs in all 3 species:** {n_all3} ({100*n_all3/n_total:.1f}%)\n\n")
            
            f.write("### Top Conserved OGs (Context-Dependent in Multiple Species)\n\n")
            top_conserved = comparison_df[comparison_df['n_species'] >= 2].head(20)
            
            if not top_conserved.empty:
                f.write("| OG ID | # Species | Mean Strength | A. pulchra | P. evermanni | P. tuahiniensis |\n")
                f.write("|-------|-----------|---------------|------------|--------------|------------------|\n")
                
                for _, row in top_conserved.iterrows():
                    f.write(f"| {row['og_id']} | {row['n_species']} | {row['mean_context_strength']:.3f} | ")
                    f.write(f"{row.get('apul_context_strength', np.nan):.3f} | ")
                    f.write(f"{row.get('peve_context_strength', np.nan):.3f} | ")
                    f.write(f"{row.get('ptua_context_strength', np.nan):.3f} |\n")
                f.write("\n")
        
        f.write("## Interpretation\n\n")
        f.write("OGs (orthologous groups) that show context-dependent regulation in multiple species suggest:\n\n")
        f.write("1. **Conserved regulatory mechanisms** - These genes may have fundamental roles where regulation is evolutionarily maintained\n")
        f.write("2. **Robust biological signals** - Conservation across species reduces the likelihood of false positives\n")
        f.write("3. **Candidates for functional validation** - These genes are high-priority targets for experimental follow-up\n\n")
        
        f.write("## Files Generated\n\n")
        f.write(f"- `cross_species_methylation_mirna_comparison.csv` - Full comparison table\n")
        f.write(f"- `plots/species_comparison_overview.png` - Visualization of cross-species patterns\n")
        f.write(f"- `cross_species_comparison_report.md` - This report\n")
    
    print(f"✅ Saved comparison report: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run context-dependent analysis on all biomin species and compare results"
    )
    parser.add_argument(
        '--skip-analysis',
        action='store_true',
        help='Skip running per-species analysis (use most recent existing outputs)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory for comparison results'
    )
    parser.add_argument(
        '--species-outputs',
        type=str,
        nargs=3,
        metavar=('APUL_DIR', 'PEVE_DIR', 'PTUA_DIR'),
        help='Specify existing output directories for each species'
    )
    
    args = parser.parse_args()
    
    workspace_root = get_workspace_root()
    
    # Create comparison output directory
    if args.output_dir:
        comparison_output = args.output_dir
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        comparison_output = os.path.join(workspace_root, f"output/biomin_comparison_{timestamp}")
    
    os.makedirs(comparison_output, exist_ok=True)
    
    print("="*70)
    print("🧬 CROSS-SPECIES BIOMIN CONTEXT-DEPENDENT ANALYSIS")
    print("="*70)
    print(f"Workspace: {workspace_root}")
    print(f"Comparison output: {comparison_output}")
    print("="*70)
    
    # Track output directories for each species
    species_output_dirs = {}
    
    if args.species_outputs:
        # Use provided output directories
        species_output_dirs = {
            'apul': args.species_outputs[0],
            'peve': args.species_outputs[1],
            'ptua': args.species_outputs[2]
        }
    elif args.skip_analysis:
        # Find most recent outputs (user will need to specify or we use latest)
        print("\n⚠️  --skip-analysis specified. Looking for existing output directories...")
        print("   For best results, use --species-outputs to specify exact directories.\n")
        
        # This is a simplified approach - in practice you'd want better tracking
        output_base = os.path.join(workspace_root, 'output')
        dirs = sorted([d for d in os.listdir(output_base) 
                      if d.startswith('context_dependent_analysis_')], reverse=True)
        
        if len(dirs) >= 3:
            # Assume last 3 are for our species (user should verify)
            for i, species in enumerate(['apul', 'peve', 'ptua']):
                if i < len(dirs):
                    species_output_dirs[species] = os.path.join(output_base, dirs[i])
                    print(f"  Using {dirs[i]} for {species}")
        else:
            print("  ❌ Not enough existing output directories found. Run without --skip-analysis.")
            return 1
    else:
        # Run analysis for each species
        for species_key in SPECIES_CONFIG.keys():
            output_dir = run_species_analysis(species_key, workspace_root)
            if output_dir:
                species_output_dirs[species_key] = output_dir
            else:
                print(f"  ⚠️ Skipping {species_key} due to analysis failure")
    
    if len(species_output_dirs) < 2:
        print("\n❌ Need at least 2 species with successful analysis for comparison")
        return 1
    
    # Load results for each species
    print("\n" + "="*60)
    print("LOADING SPECIES RESULTS")
    print("="*60)
    
    species_results = {}
    for species_key, output_dir in species_output_dirs.items():
        print(f"\n📂 Loading {SPECIES_CONFIG[species_key]['name']}...")
        species_results[species_key] = load_species_results(output_dir, species_key)
    
    # Compare results across species
    comparison_df = compare_methylation_mirna_context(species_results, comparison_output)
    
    # Generate visualizations
    if not comparison_df.empty:
        print("\n" + "="*60)
        print("GENERATING VISUALIZATIONS")
        print("="*60)
        plot_species_comparison(comparison_df, species_results, comparison_output)
    
    # Generate report
    generate_comparison_report(comparison_df, species_results, comparison_output)
    
    print("\n" + "="*70)
    print("🎉 CROSS-SPECIES COMPARISON COMPLETE!")
    print("="*70)
    print(f"Results saved to: {comparison_output}")
    print("="*70)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
