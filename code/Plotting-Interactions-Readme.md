# Plotting Interaction Expression Trajectories

This directory contains the script `plot_selected_interactions_line_plots.py` which generates a variety of time‑series expression visualizations for genes and their regulators across four time points (T1–T4). The script ingests pre‑processed count matrices (genes, lncRNA, miRNA, methylation) plus an interactions table (`multi_way_interactions.csv`) and produces several plot types and summary CSVs.

## Data Assumptions

- Cleaned count tables live at `data/cleaned_datasets/` with filenames:
  - `gene_counts_cleaned.csv`
  - `lncrna_counts_cleaned.csv`
  - `mirna_counts_cleaned.csv`
  - `wgbs_counts_cleaned.csv` (methylation)
- Columns follow the pattern `<CONDITION>-TPx` where `x` = 1..4.
- Interaction file (`multi_way_interactions.csv`) includes columns:
  - `gene`
  - `improvement_from_regulators`
  - `regulator_types` (serialized list with prefixes `mirna_`, `lncrna_`, `methylation_`)

## Plot Types Generated

1. Combined z‑score line plot for a gene + its regulators (type or entity colors).
2. Faceted z‑score line plots by regulator type (panel per type).
3. Optional per‑entity small multiples with per‑condition trajectories (`--include-samples`).
4. Optional condition facets (panel per condition) showing all entities (`--facet-by-condition`).
5. Expression and summary CSV exports.

All plots use per‑entity z‑scores (mean 0, SD 1 across its own time points; for sample-level plots, across all condition/timepoint combinations) to emphasize trajectory shape over absolute magnitude.

## Invocation Modes

You can drive the script in several mutually exclusive modes controlling which genes / regulators are plotted:

| Mode | Flags | Description |
|------|-------|-------------|
| Preview subset (default) | (no selection flags) | Plots the first 8 interaction rows from the interactions file. |
| Single gene (no regulators) | `--single-gene <ID or suffix>` | Plots a gene found in the gene counts index whose full ID equals or ends with the provided string; regulators not loaded. Useful for quick inspection. |
| Specific gene with regulators | `--gene-with-regulators <gene_id>` | Extracts that gene and its regulator list from the interactions file. |
| Max improvement auto-pick | `--auto-max-improvement` | Automatically selects the interaction row with the largest `improvement_from_regulators` value. |

Only one of `--single-gene`, `--gene-with-regulators`, or `--auto-max-improvement` should be supplied. Priority order:

1. `--auto-max-improvement`
2. `--gene-with-regulators`
3. `--single-gene`
4. Default subset preview

## Command-Line Arguments

```text
usage: plot_selected_interactions_line_plots.py [-h]
  [--include-samples] [--facet-by-condition]
  [--single-gene SINGLE_GENE]
  [--color-by-entity]
  [--gene-with-regulators GENE_WITH_REGULATORS]
  [--auto-max-improvement]
  [--interactions-path INTERACTIONS_PATH]
```

### Flags & Parameters

- `--include-samples`  Generate per‑entity small multiples (all conditions + mean).
- `--facet-by-condition`  Produce panel-per-condition figure with all entities.
- `--single-gene <STRING>`  Plot a gene by exact ID or suffix (no regulators).
- `--color-by-entity`  Unique color per non‑gene entity in combined plot.
- `--gene-with-regulators <GENE_ID>`  Use that gene's interaction row regulators.
- `--auto-max-improvement`  Auto-pick gene with max improvement.
- `--interactions-path <PATH>`  Custom path to interactions CSV.

### Output Structure

Creates a timestamped directory under `output/` with pattern:

```text
output/<mode>_line_plots_<YYYYMMDD_HHMMSS>/
  summary.csv
  row<n>_<GENE>_expression.csv
  row<n>_<GENE>_line_plot_zscore[_entitycolors].png
  row<n>_<GENE>_facet_zscore.png
  (optional) row<n>_<GENE>_samples_small_multiples.png
  (optional) row<n>_<GENE>_condition_facets.png
```
`summary.csv` columns:

- row
- gene
- n_regulators_found
- csv (path to expression CSV)
- plot_zscore
- color_mode ("type" or "entity")
- facet_plot_zscore
- samples_small_multiples
- condition_facets

### Examples

1. Default (first 8 interactions, type colors):

  ```bash
  python3 code/plot_selected_interactions_line_plots.py
  ```

1. Single gene by suffix (no regulators) + per‑sample detail:

  ```bash
  python3 code/plot_selected_interactions_line_plots.py --single-gene FUN_022503 --include-samples --facet-by-condition
  ```

1. Specific gene with its regulators, entity colors:

  ```bash
  python3 code/plot_selected_interactions_line_plots.py --gene-with-regulators FUN_033824 --color-by-entity
  ```

1. Auto-select max improvement gene with all extras:

  ```bash
  python3 code/plot_selected_interactions_line_plots.py --auto-max-improvement --include-samples --facet-by-condition --color-by-entity
  ```

1. Custom interactions file:

  ```bash
  python3 code/plot_selected_interactions_line_plots.py --interactions-path path/to/multi_way_interactions.csv --auto-max-improvement
  ```

## Color Logic

- Default: fixed palette by type
  - gene: black (#000000)
  - miRNA: blue (#1f77b4)
  - lncRNA: green (#2ca02c)
  - methylation: red (#d62728)
- Entity mode (`--color-by-entity`): unique color per regulator (Tab20 then HSV fallback), gene remains black.

 
## Z-Score Definition

For combined and type-facet plots: z-score per entity across its four timepoint means.
For small multiples & condition facets: z-score per entity across all condition+timepoint values (shared scaling across panels within each figure).

## Edge Cases & Notes

- Regulators absent from their respective count matrix indices are silently skipped.
- Entities with zero variance receive z-scores of 0 (std is replaced with 1 to avoid divide-by-zero).
- Legend truncation occurs if too many entities (>25) in combined plot to prevent clutter.
- Interaction file regulator list may be double-encoded; script attempts nested `ast.literal_eval` parsing.
- If using `--single-gene` only the gene itself is plotted (no regulators) unless you choose a mode that loads regulators.

## Extending

Potential additions (not yet implemented):

- Toggle for raw expression (non z-score).
- Export of z-score matrices.
- Filtering regulators by type via flags.

## Troubleshooting

| Issue | Cause | Remedy |
|-------|-------|--------|
| "Interactions file not found" | Wrong default path or moved output | Provide `--interactions-path` explicitly. |
| Gene not found for `--single-gene` | Suffix doesn’t match any index entry | Confirm exact ID in `gene_counts_cleaned.csv`. |
| Empty plots (only gene) | Regulators missing from count tables | Verify IDs / preprocessing steps. |
| Legend very long | Many regulators | Use facets or rely on small multiples for detail. |

## License

See root `LICENSE` file.
