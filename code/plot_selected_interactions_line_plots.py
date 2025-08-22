import os
import ast
import re
import argparse
from datetime import datetime
from typing import List, Tuple

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def load_counts(base_dir: str):
    genes = pd.read_csv(os.path.join(base_dir, 'data', 'cleaned_datasets', 'gene_counts_cleaned.csv'), index_col=0)
    lncrna = pd.read_csv(os.path.join(base_dir, 'data', 'cleaned_datasets', 'lncrna_counts_cleaned.csv'), index_col=0)
    mirna = pd.read_csv(os.path.join(base_dir, 'data', 'cleaned_datasets', 'mirna_counts_cleaned.csv'), index_col=0)
    methyl = pd.read_csv(os.path.join(base_dir, 'data', 'cleaned_datasets', 'wgbs_counts_cleaned.csv'), index_col=0)
    return genes, lncrna, mirna, methyl


def aggregate_timepoints(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse replicates per timepoint to mean and rename TPx -> Tx."""
    tp_labels = ["TP1", "TP2", "TP3", "TP4"]
    agg = {}
    for tp in tp_labels:
        cols = [c for c in df.columns if c.endswith(tp)]
        if not cols:
            raise ValueError(f"No columns for {tp}")
        agg['T'+tp[-1]] = df[cols].mean(axis=1)
    return pd.DataFrame(agg)


def parse_interaction_rows(interactions_path: str, row_indices: List[int]) -> List[Tuple[str, List[str]]]:
    results = []
    with open(interactions_path, 'r') as f:
        _ = next(f)  # header
        for i, line in enumerate(f):
            if i not in row_indices:
                continue
            parts = line.rstrip('\n').split(',', 7)
            if len(parts) < 8:
                continue
            gene = parts[0]
            reg_field = parts[7]
            regs = []
            try:
                parsed = ast.literal_eval(reg_field)
                if isinstance(parsed, str) and parsed.strip().startswith('['):
                    parsed = ast.literal_eval(parsed)
                if isinstance(parsed, list):
                    regs = parsed
            except Exception:
                pass
            results.append((gene, regs))
    return results


def collect_expression(gene: str, regulators: List[str], genes_tp: pd.DataFrame, lncrna_tp: pd.DataFrame, mirna_tp: pd.DataFrame, methyl_tp: pd.DataFrame) -> pd.DataFrame:
    records = []
    def add_row(eid: str, etype: str, source: pd.DataFrame):
        if eid not in source.index:
            return
        for tp in ["T1", "T2", "T3", "T4"]:
            records.append({"entity": eid, "type": etype, "timepoint": tp, "expression": source.loc[eid, tp]})

    # gene
    add_row(gene, 'gene', genes_tp)
    # regulators
    for r in regulators:
        if r.startswith('mirna_'):
            cluster = r.replace('mirna_', '')  # Cluster_XXXX
            add_row(cluster, 'miRNA', mirna_tp)
        elif r.startswith('lncrna_'):
            lncrna_id = r.replace('lncrna_', '')
            add_row(lncrna_id, 'lncRNA', lncrna_tp)
        elif r.startswith('methylation_'):
            cpg_id = r.replace('methylation_', '')
            add_row(cpg_id, 'methylation', methyl_tp)
    return pd.DataFrame(records)


def collect_expression_samples(gene: str, regulators: List[str], genes_df: pd.DataFrame, lncrna_df: pd.DataFrame, mirna_df: pd.DataFrame, methyl_df: pd.DataFrame) -> pd.DataFrame:
    """Return long dataframe with raw per-sample expression for gene and regulators.

    Columns: entity, type, condition, timepoint, expression
    condition = sample prefix without trailing -TPx
    timepoint = T1..T4
    """
    sample_pattern = re.compile(r"^(.*)-TP([1-4])$")

    def melt_entity(eid: str, etype: str, source: pd.DataFrame) -> pd.DataFrame:
        if eid not in source.index:
            return pd.DataFrame()
        row = source.loc[eid]
        records = []
        for col, val in row.items():
            m = sample_pattern.match(col)
            if not m:
                continue
            condition = m.group(1)
            tp = 'T' + m.group(2)
            records.append({
                'entity': eid,
                'type': etype,
                'condition': condition,
                'timepoint': tp,
                'expression': val
            })
        return pd.DataFrame(records)

    dfs = [melt_entity(gene, 'gene', genes_df)]
    for r in regulators:
        if r.startswith('mirna_'):
            dfs.append(melt_entity(r.replace('mirna_', ''), 'miRNA', mirna_df))
        elif r.startswith('lncrna_'):
            dfs.append(melt_entity(r.replace('lncrna_', ''), 'lncRNA', lncrna_df))
        elif r.startswith('methylation_'):
            dfs.append(melt_entity(r.replace('methylation_', ''), 'methylation', methyl_df))
    out = pd.concat([d for d in dfs if not d.empty], ignore_index=True)
    return out


def plot_line_chart(expr_df: pd.DataFrame, gene: str, out_dir: str, row_label: str, color_by_entity: bool = False):
    """Plot per-entity z-scores across timepoints.

    Parameters
    ----------
    expr_df : DataFrame with columns entity, type, timepoint, expression
    gene : str gene identifier
    out_dir : output directory
    row_label : str row number label
    color_by_entity : if True, assign a unique color to every non-gene entity instead of coloring by type.
    """
    # Compute z-scores per entity
    wide = expr_df.pivot_table(index='entity', columns='timepoint', values='expression', aggfunc='mean')
    # Ensure consistent column order
    wide = wide[[c for c in ['T1', 'T2', 'T3', 'T4'] if c in wide.columns]]
    zwide = wide.subtract(wide.mean(axis=1), axis=0)
    stds = wide.std(axis=1).replace(0, 1)
    zwide = zwide.divide(stds, axis=0)
    # Melt back
    zdf = zwide.reset_index().melt(id_vars='entity', var_name='timepoint', value_name='zscore')
    # Merge type info
    types = expr_df[['entity', 'type']].drop_duplicates()
    zdf = zdf.merge(types, on='entity', how='left')
    tp_order = ["T1", "T2", "T3", "T4"]
    zdf['timepoint'] = pd.Categorical(zdf['timepoint'], categories=tp_order, ordered=True)

    plt.figure(figsize=(8.5, 5))
    # Build color mapping
    if color_by_entity:
        import math
        from matplotlib import cm
        entities = [e for e in zdf['entity'].unique() if e != gene]
        n = max(len(entities), 1)
        # Use tab20 first then fallback to hsv if more needed
        base_cmap = cm.get_cmap('tab20')
        colors = [base_cmap(i % base_cmap.N) for i in range(n)]
        if n > base_cmap.N:
            extra_cmap = cm.get_cmap('hsv')
            colors = [extra_cmap(i / n) for i in range(n)]
        entity_colors = {ent: colors[i] for i, ent in enumerate(entities)}
        type_colors = {'gene': '#000000'}  # only used for gene
    else:
        type_colors = {
            'gene': '#000000',
            'miRNA': '#1f77b4',
            'lncRNA': '#2ca02c',
            'methylation': '#d62728'
        }
        entity_colors = {}
    # Gene first (thicker, always black)
    gene_df = zdf[zdf['entity'] == gene]
    if not gene_df.empty:
        plt.plot(gene_df['timepoint'], gene_df['zscore'], marker='o', linewidth=3, color='#000000', label=f"{gene} (gene)")
    # Regulators
    if color_by_entity:
        regulators_df = zdf[zdf['entity'] != gene]
        for ent in regulators_df['entity'].unique():
            edat = regulators_df[regulators_df['entity'] == ent]
            etype = edat['type'].iloc[0]
            plt.plot(edat['timepoint'], edat['zscore'], marker='o', linewidth=1.4, alpha=0.9, color=entity_colors.get(ent, '#555555'), label=f"{ent} ({etype})")
    else:
        for etype in ['miRNA', 'lncRNA', 'methylation']:
            sub = zdf[zdf['type'] == etype]
            for ent in sub['entity'].unique():
                if ent == gene:
                    continue
                edat = sub[sub['entity'] == ent]
                plt.plot(edat['timepoint'], edat['zscore'], marker='o', linewidth=1, alpha=0.7, color=type_colors[etype], label=f"{ent} ({etype})")

    handles, labels = plt.gca().get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    max_labels = 25
    disp_items = list(unique.items())
    truncated = False
    if len(disp_items) > max_labels:
        truncated = True
        disp_items = [i for i in disp_items if '(gene)' in i[0]] + disp_items[1:max_labels]
    plt.legend([h for _, h in disp_items], [l for l, _ in disp_items], fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')
    title_mode = 'entity colors' if color_by_entity else 'type colors'
    plt.title(f"Row {row_label}: {gene} and regulators (z-score, {title_mode})")
    plt.xlabel('Time Point')
    plt.ylabel('Z-score (per entity)')
    if truncated:
        plt.suptitle('Legend truncated', fontsize=8, y=0.02)
    plt.axhline(0, color='grey', linewidth=0.8, linestyle='--')
    plt.tight_layout()
    fname = f"row{row_label}_{gene}_line_plot_zscore{'_entitycolors' if color_by_entity else ''}.png".replace('/', '_')
    out_path = os.path.join(out_dir, fname)
    plt.savefig(out_path, dpi=300)
    plt.close()
    return out_path


def plot_facet_by_type(expr_df: pd.DataFrame, gene: str, out_dir: str, row_label: str):
    """Create faceted line plots (z-scores) per regulator type with gene overlaid in each panel."""
    # Reuse z-score transformation logic
    wide = expr_df.pivot_table(index='entity', columns='timepoint', values='expression', aggfunc='mean')
    wide = wide[[c for c in ['T1', 'T2', 'T3', 'T4'] if c in wide.columns]]
    zwide = wide.subtract(wide.mean(axis=1), axis=0)
    stds = wide.std(axis=1).replace(0, 1)
    zwide = zwide.divide(stds, axis=0)
    zdf = zwide.reset_index().melt(id_vars='entity', var_name='timepoint', value_name='zscore')
    types = expr_df[['entity', 'type']].drop_duplicates()
    zdf = zdf.merge(types, on='entity', how='left')
    tp_order = ["T1", "T2", "T3", "T4"]
    zdf['timepoint'] = pd.Categorical(zdf['timepoint'], categories=tp_order, ordered=True)

    type_colors = {
        'miRNA': '#1f77b4',
        'lncRNA': '#2ca02c',
        'methylation': '#d62728'
    }

    regulator_types = ['miRNA', 'lncRNA', 'methylation']
    n_types_present = sum([(zdf['type'] == t).any() for t in regulator_types])
    if n_types_present == 0:
        return None

    fig, axes = plt.subplots(1, n_types_present, figsize=(5 * n_types_present, 4), sharey=True)
    if n_types_present == 1:
        axes = [axes]

    idx_ax = 0
    for rtype in regulator_types:
        sub = zdf[zdf['type'] == rtype]
        if sub.empty:
            continue
        ax = axes[idx_ax]
        # Plot regulators
        for ent in sub['entity'].unique():
            edat = sub[sub['entity'] == ent]
            ax.plot(edat['timepoint'], edat['zscore'], marker='o', linewidth=1, alpha=0.7, color=type_colors[rtype], label=ent)
        # Overlay gene
        gene_df = zdf[zdf['entity'] == gene]
        if not gene_df.empty:
            ax.plot(gene_df['timepoint'], gene_df['zscore'], marker='o', linewidth=3, color='black', label=f"{gene} (gene)")
        ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
        ax.set_title(f"{rtype}")
        ax.set_xlabel('Time Point')
        if idx_ax == 0:
            ax.set_ylabel('Z-score')
        # Manage legend (limit)
        handles, labels = ax.get_legend_handles_labels()
        if len(labels) > 12:
            # Keep gene + first 11 regulators
            new_items = []
            for h, l in zip(handles, labels):
                if '(gene)' in l or len(new_items) < 12:
                    new_items.append((h, l))
                if len(new_items) >= 13:
                    break
            handles, labels = zip(*new_items)
        ax.legend(handles, labels, fontsize=7, loc='upper left', bbox_to_anchor=(1.0, 1.0))
        idx_ax += 1

    fig.suptitle(f"Row {row_label}: {gene} regulators (z-score by type)", y=1.02)
    fig.tight_layout()
    fname = f"row{row_label}_{gene}_facet_zscore.png".replace('/', '_')
    out_path = os.path.join(out_dir, fname)
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path


def plot_sample_small_multiples(sample_df: pd.DataFrame, out_dir: str, gene: str, row_label: str):
    """Plot per-sample z-score trajectories for each entity as small multiples.

    For each entity: faint lines per condition (alpha), thick mean line.
    """
    if sample_df.empty:
        return None
    # Compute z-scores per entity across all conditions+timepoints
    sample_df = sample_df.copy()
    # pivot and zscore
    z_dfs = []
    for ent, sub in sample_df.groupby('entity'):
        vals = sub['expression'].astype(float)
        mean = vals.mean()
        std = vals.std()
        if std == 0:
            std = 1.0
        sub = sub.copy()
        sub['zscore'] = (sub['expression'] - mean) / std
        z_dfs.append(sub)
    zdf = pd.concat(z_dfs, ignore_index=True)
    entities = sorted(zdf['entity'].unique())
    n = len(entities)
    # dynamic layout
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.2, nrows * 2.2), sharey=True)
    axes = axes.flatten()
    tp_order = ['T1', 'T2', 'T3', 'T4']
    for ax in axes[n:]:
        ax.axis('off')
    for idx, ent in enumerate(entities):
        ax = axes[idx]
        sub = zdf[zdf['entity'] == ent]
        # Lines per condition
        for cond, csub in sub.groupby('condition'):
            csub = csub.sort_values('timepoint')
            ax.plot(csub['timepoint'], csub['zscore'], color='#888888', linewidth=0.8, alpha=0.4)
        # Mean line
        mean_line = sub.groupby('timepoint')['zscore'].mean().reindex(tp_order)
        ax.plot(tp_order, mean_line, color='black', linewidth=2)
        # Highlight gene vs regulators differently by title style
        if ent == gene:
            ax.set_title(ent + ' (gene)', fontsize=9, fontweight='bold')
        else:
            etype = sub['type'].iloc[0]
            ax.set_title(f"{ent} ({etype})", fontsize=8)
        ax.axhline(0, color='grey', linewidth=0.6, linestyle='--')
        ax.set_xticks(tp_order)
        if idx % ncols == 0:
            ax.set_ylabel('z')
    fig.suptitle(f"Row {row_label}: per-sample z-score trajectories", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_path = os.path.join(out_dir, f"row{row_label}_{gene}_samples_small_multiples.png".replace('/', '_'))
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


def plot_condition_facets(sample_df: pd.DataFrame, out_dir: str, gene: str, row_label: str):
    """Facet by condition (sample): each panel shows z-score trajectories (T1-T4) for gene + regulators.

    Z-score is computed per entity over all conditions/timepoints (same scaling across facets).
    """
    if sample_df.empty:
        return None
    df = sample_df.copy()
    # Compute z-score per entity
    z_parts = []
    for ent, sub in df.groupby('entity'):
        vals = sub['expression'].astype(float)
        mean = vals.mean()
        std = vals.std()
        if std == 0:
            std = 1.0
        sub = sub.copy()
        sub['zscore'] = (sub['expression'] - mean) / std
        z_parts.append(sub)
    zdf = pd.concat(z_parts, ignore_index=True)
    # Order timepoints
    tp_order = ['T1', 'T2', 'T3', 'T4']
    zdf['timepoint'] = pd.Categorical(zdf['timepoint'], categories=tp_order, ordered=True)
    conditions = sorted(zdf['condition'].unique())
    n_cond = len(conditions)
    # Layout heuristic: up to 8 columns
    ncols = 8 if n_cond >= 32 else 6 if n_cond >= 24 else 5 if n_cond >= 15 else 4
    nrows = (n_cond + ncols - 1) // ncols
    import math
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.4, nrows * 2.2), sharey=True)
    axes = axes.flatten()
    color_map = {'gene': '#000000', 'miRNA': '#1f77b4', 'lncRNA': '#2ca02c', 'methylation': '#d62728'}
    # Pre-split for speed
    by_entity = {ent: sub.sort_values('timepoint') for ent, sub in zdf.groupby('entity')}
    type_map = {ent: sub['type'].iloc[0] for ent, sub in zdf.groupby('entity')}
    # Determine global y-limits (clip extremes to keep readability)
    all_z = zdf['zscore']
    low, high = all_z.quantile(0.01), all_z.quantile(0.99)
    y_min = float(min(low, -2.5))
    y_max = float(max(high, 2.5))
    for ax in axes[n_cond:]:
        ax.axis('off')
    for i, cond in enumerate(conditions):
        ax = axes[i]
        subset = zdf[zdf['condition'] == cond]
        # Plot regulators first
        for ent, esub in subset.groupby('entity'):
            ent_type = type_map.get(ent, '')
            line_width = 3 if ent == gene else 1
            alpha = 1.0 if ent == gene else 0.6
            color = color_map.get(ent_type, '#555555') if ent != gene else '#000000'
            esub = esub.sort_values('timepoint')
            ax.plot(esub['timepoint'], esub['zscore'], marker='o', linewidth=line_width, alpha=alpha, color=color)
        ax.set_title(cond, fontsize=8, fontweight='bold')
        if i % ncols == 0:
            ax.set_ylabel('z')
        ax.set_ylim(y_min, y_max)
        ax.axhline(0, color='grey', linewidth=0.6, linestyle='--')
        ax.set_xticks(tp_order)
        ax.tick_params(axis='x', labelsize=7)
        ax.tick_params(axis='y', labelsize=7)
    # Adjust layout first to avoid tight_layout cutting off legend
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    # Shared legend inside reserved top area
    legend_elements = [
        Line2D([0], [0], color='#000000', lw=2, label='Gene'),
        Line2D([0], [0], color='#1f77b4', lw=2, label='miRNA'),
        Line2D([0], [0], color='#2ca02c', lw=2, label='lncRNA'),
        Line2D([0], [0], color='#d62728', lw=2, label='Methylation')
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.965), fontsize=9)
    fig.suptitle(f"Row {row_label}: Faceted by condition (gene + regulators z-score)", y=0.995)
    out_path = os.path.join(out_dir, f"row{row_label}_{gene}_condition_facets.png".replace('/', '_'))
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


def main():
    parser = argparse.ArgumentParser(description='Generate interaction expression plots.')
    parser.add_argument('--include-samples', action='store_true', help='Include per-sample small multiple plots (entities).')
    parser.add_argument('--facet-by-condition', action='store_true', help='Facet plots by condition (sample) instead of by regulator type (adds extra figure).')
    parser.add_argument('--single-gene', type=str, help='Run plots for a single gene ID or suffix (e.g., FUN_022503).')
    parser.add_argument('--color-by-entity', action='store_true', help='Use a unique color per entity instead of per type in the combined line plot.')
    parser.add_argument('--gene-with-regulators', type=str, help='Plot a specific gene using its regulators from the interactions CSV.')
    parser.add_argument('--auto-max-improvement', action='store_true', help='Automatically select the gene with the highest improvement_from_regulators.')
    parser.add_argument('--interactions-path', type=str, help='Path to multi_way_interactions.csv (optional).')
    args = parser.parse_args()
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    genes, lncrna, mirna, methyl = load_counts(base_dir)

    # Determine interactions list: list of (gene, [regulators])
    interactions = []
    # Resolve interactions file path
    default_interactions_path = os.path.join(base_dir, 'output', 'subset_context_dependent_analysis_20250818_095433', 'tables', 'multi_way_interactions.csv')
    interactions_path = args.interactions_path or default_interactions_path

    if args.auto_max_improvement:
        if not os.path.exists(interactions_path):
            print(f"Interactions file not found: {interactions_path}")
            return
        df_int = pd.read_csv(interactions_path)
        if 'improvement_from_regulators' not in df_int.columns:
            print('Column improvement_from_regulators missing in interactions file.')
            return
        top_row = df_int.loc[df_int['improvement_from_regulators'].idxmax()]
        gene = top_row['gene']
        regs_raw = top_row.get('regulator_types', '[]')
        try:
            regs = ast.literal_eval(regs_raw)
            if isinstance(regs, str) and regs.startswith('['):
                regs = ast.literal_eval(regs)
            if not isinstance(regs, list):
                regs = []
        except Exception:
            regs = []
        interactions = [(gene, regs)]
        print(f"Auto-selected gene with max improvement: {gene} (improvement={top_row['improvement_from_regulators']}) with {len(regs)} regulators.")
    elif args.gene_with_regulators:
        if not os.path.exists(interactions_path):
            print(f"Interactions file not found: {interactions_path}")
            return
        df_int = pd.read_csv(interactions_path)
        sub = df_int[df_int['gene'] == args.gene_with_regulators]
        if sub.empty:
            print(f"Gene {args.gene_with_regulators} not found in interactions file.")
            return
        row = sub.iloc[0]
        regs_raw = row.get('regulator_types', '[]')
        try:
            regs = ast.literal_eval(regs_raw)
            if isinstance(regs, str) and regs.startswith('['):
                regs = ast.literal_eval(regs)
            if not isinstance(regs, list):
                regs = []
        except Exception:
            regs = []
        interactions = [(row['gene'], regs)]
        print(f"Loaded gene {row['gene']} with {len(regs)} regulators from interactions file (improvement={row.get('improvement_from_regulators', 'NA')}).")
    elif args.single_gene:
        gene_suffix = args.single_gene
        gene_ids = [g for g in genes.index if g == gene_suffix or g.endswith(gene_suffix)]
        if not gene_ids:
            print(f"No gene found matching '{gene_suffix}'. Exiting.")
            return
        interactions = [(g, []) for g in gene_ids]
        print(f"Found {len(gene_ids)} gene(s) for pattern '{gene_suffix}': {gene_ids}")
    else:
        selected_rows = list(range(0, 8))  # default preview subset
        interactions = parse_interaction_rows(interactions_path, selected_rows)
    genes_tp = aggregate_timepoints(genes)
    lncrna_tp = aggregate_timepoints(lncrna)
    mirna_tp = aggregate_timepoints(mirna)
    methyl_tp = aggregate_timepoints(methyl)

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    if args.auto_max_improvement:
        out_tag = 'max_improvement'
    elif args.gene_with_regulators:
        out_tag = 'gene_with_regulators'
    elif args.single_gene:
        out_tag = 'single_gene'
    else:
        out_tag = 'selected_rows'
    out_dir = os.path.join(base_dir, 'output', f'{out_tag}_line_plots_{timestamp}')
    os.makedirs(out_dir, exist_ok=True)

    summary_records = []
    for idx, (gene, regulators) in enumerate(interactions, start=1):
        expr_df = collect_expression(gene, regulators, genes_tp, lncrna_tp, mirna_tp, methyl_tp)
        csv_path = os.path.join(out_dir, f'row{idx}_{gene}_expression.csv'.replace('/', '_'))
        expr_df.to_csv(csv_path, index=False)
        plot_path = plot_line_chart(expr_df, gene, out_dir, str(idx), color_by_entity=args.color_by_entity)
        facet_path = plot_facet_by_type(expr_df, gene, out_dir, str(idx))
        sample_plot_path = ''
        condition_facet_path = ''
        sample_df = None
        if args.include_samples or args.facet_by_condition:
            sample_df = collect_expression_samples(gene, regulators, genes, lncrna, mirna, methyl)
        if args.include_samples and sample_df is not None:
            sample_plot_path = plot_sample_small_multiples(sample_df, out_dir, gene, str(idx)) or ''
        if args.facet_by_condition and sample_df is not None:
            condition_facet_path = plot_condition_facets(sample_df, out_dir, gene, str(idx)) or ''
        summary_records.append({
            'row': idx,
            'gene': gene,
            'n_regulators_found': expr_df['entity'].nunique() - 1,
            'csv': csv_path,
            'plot_zscore': plot_path,
            'color_mode': 'entity' if args.color_by_entity else 'type',
            'facet_plot_zscore': facet_path if facet_path else '',
            'samples_small_multiples': sample_plot_path,
            'condition_facets': condition_facet_path
        })
        print(f"Generated row {idx} plot: {plot_path}")
        if facet_path:
            print(f"Generated row {idx} facet plot: {facet_path}")
        if sample_plot_path:
            print(f"Generated row {idx} sample small multiples: {sample_plot_path}")
        if condition_facet_path:
            print(f"Generated row {idx} condition facets: {condition_facet_path}")

    pd.DataFrame(summary_records).to_csv(os.path.join(out_dir, 'summary.csv'), index=False)
    print(f"All done. Output directory: {out_dir}")


if __name__ == '__main__':
    main()
