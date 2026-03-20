# Gene Expression in Human Developmental Brain

A Shiny web application for visualizing gene expression throughout the development of the human brain. This tool leverages the **BrainSpan Atlas of the Developing Human Brain** dataset and provides a spatial interface using the `ggseg` package to map expression levels onto the Desikan-Killiany (DK) atlas.

## 🚀 Features

- **Gene Expression Mapping**: Search for any HGNC gene name to see its spatial expression pattern across 5 broad developmental stages (1st Trimester to Adult).
- **GSVA-Enrichr Mode**: Upload a custom list of genes (gene signature) to calculate and visualize its activity across brain regions using Single Sample GSEA (ssGSEA).
- **High-Quality Exports**: Download brain maps in SVG or TIFF formats for publications.
- **Interactive UI**: Built with a modern, responsive interface using `shinythemes` (Cerulean).

## 🛠 Prerequisites

The application requires R (≥ 4.3.0) and several system-level dependencies for spatial data processing.

### R Packages

Install the required packages from CRAN, Bioconductor, and GitHub:

```R
# CRAN
install.packages(c("shiny", "dplyr", "tidyr", "ggplot2", "shinythemes", "svglite", "shinycssloaders", "remotes", "ggseg"))

# Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GSVA")
```

## 📂 Data Preparation

Before running the app, you must process the raw BrainSpan CSV files into an optimized format.

1.  Ensure the `genes_matrix_csv/` directory contains:
    - `expression_matrix.csv`
    - `rows_metadata.csv`
    - `columns_metadata.csv`
2.  Run the preparation script:
    ```bash
    Rscript preparaDados.R
    ```
    This script generates `dados_otimizados.rds`, which includes mapped brain regions and pre-calculated age groups.

## 🏃 Running the Application

Once the data is prepared, you can launch the app:

```R
shiny::runApp()
```

## 🐳 Docker Support

A `Dockerfile` is provided for containerized deployment. To build and run:

```bash
# Build the image
docker build -t devbrain-markers .

# Run the container
docker run -p 3838:3838 devbrain-markers
```

## 📊 Data Source

Data is sourced from the [BrainSpan Atlas of the Developing Human Brain](https://www.brainspan.org/), specifically the RNA-Seq RPKM values averaged to genes.

## 🤝 Acknowledgements

Made by the **[MCB Lab at UFRN](https://mcblab.github.io/)**.

If you use this tool in your research, please cite the BrainSpan Atlas and the `ggseg` package.
