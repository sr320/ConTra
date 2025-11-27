# ConTra Facet Z-Score Shiny App

This Shiny app provides an interactive web interface for visualizing gene regulator relationships using facet z-score plots, as described in issue #14.

## Features

- **Gene Selection**: Choose from available genes in the dataset via dropdown
- **Interactive Plotting**: Displays facet z-score plots showing regulator types (miRNA, lncRNA, methylation) in separate panels
- **Gene Overlay**: Option to overlay the selected gene in all panels
- **Customizable Display**: Adjustable plot height and display options
- **Data Table**: View underlying expression data
- **Plot Information**: Summary statistics about regulators found

## Files

- `app.R` - Main Shiny application file
- `data_functions.R` - Data loading and processing functions
- `install_packages.R` - Script to install required R packages
- `README.md` - This documentation

## Requirements

### R Packages
- shiny
- ggplot2
- dplyr
- readr
- tidyr

### Data Files
The app expects the ConTra data files to be available in the parent directory structure:
```
../data/cleaned_datasets/
  ├── gene_counts_cleaned.csv
  ├── lncrna_counts_cleaned.csv
  ├── mirna_counts_cleaned.csv
  └── wgbs_counts_cleaned.csv
```

## Installation

1. Install required R packages:
```bash
sudo apt install r-cran-shiny r-cran-ggplot2 r-cran-dplyr r-cran-readr r-cran-tidyr
```

Or run the installation script:
```bash
Rscript install_packages.R
```

2. Ensure data files are in the correct location relative to the app directory.

## Running the App

### Local Development
```r
# In R console, from the shiny_app directory
shiny::runApp()
```

### Command Line
```bash
# From the shiny_app directory
R -e "shiny::runApp()"
```

### Specify Host and Port
```r
# To run on all interfaces (for server deployment)
shiny::runApp(host = "0.0.0.0", port = 3838)
```

## How It Works

1. **Data Loading**: The app loads the four ConTra datasets (genes, lncRNA, miRNA, methylation) on startup
2. **Timepoint Aggregation**: Replicates the Python logic to aggregate timepoints from TP1-TP4 to T1-T4
3. **Regulator Selection**: For demonstration, the app selects sample regulators from each type (miRNA, lncRNA, methylation)
4. **Z-score Calculation**: Computes z-scores per entity across timepoints (matching the Python implementation)
5. **Facet Plotting**: Creates ggplot2 faceted plots with separate panels for each regulator type

## Extending the App

To make this a full-featured implementation:

1. **Load Interaction Data**: Integrate with `multi_way_interactions.csv` to get actual regulators for each gene
2. **Advanced Filtering**: Add options to filter regulators by significance, type, or other criteria  
3. **Download Features**: Allow users to download plots and data
4. **Multiple Plot Types**: Add other plot types from the Python script (condition facets, small multiples)
5. **Deployment**: Deploy to a server using Shiny Server or shinyapps.io

## Plot Description

The facet z-score plot shows:
- **Separate panels** for each regulator type (miRNA, lncRNA, methylation)
- **Z-score trajectories** across timepoints T1-T4
- **Gene overlay** (thick black line) showing the selected gene's trajectory
- **Regulator lines** (colored by type) showing related regulators
- **Zero reference line** (dashed gray) for z-score context

This matches the functionality of the `plot_facet_by_type()` function in the original Python code.

## Troubleshooting

- **Data not loading**: Check that data files exist in `../data/cleaned_datasets/`
- **Package errors**: Ensure all required R packages are installed
- **Performance**: The app loads 100 genes by default; modify `get_gene_choices()` for more/fewer options
- **No regulator data**: The current implementation uses sample regulators; integrate interaction data for gene-specific regulators