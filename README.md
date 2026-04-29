# Gene Expression in Human Developmental Brain

A Shiny web application for visualizing gene expression throughout the development of the human brain. This tool leverages the **BrainSpan Atlas of the Developing Human Brain** dataset and provides a spatial interface using the `ggseg` package to map expression levels onto the Desikan-Killiany (DK) atlas.

## 🚀 Features

- **Gene Expression Mapping**: Search for any HGNC gene name to see its spatial expression pattern across 5 broad developmental stages (1st Trimester to Adult).
- **GSVA-Enrichr Mode**: Upload a custom list of genes (gene signature) to calculate and visualize its activity across brain regions using GSVA.
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

Before running the application, you must process the raw BrainSpan data and calculate the ontology enrichment scores. This is a two-step process based on the provided scripts:

### 1. Prepare the Expression Matrix
Ensure the `genes_matrix_csv/` directory contains:
- `expression_matrix.csv`
- `rows_metadata.csv`
- `columns_metadata.csv`

Run the preparation script from your terminal:
```bash
Rscript preparaDados.R
```
This script generates `dados_otimizados.rds`, which includes mapped brain regions, filtered gene lists, and pre-calculated age groups.

### 2. Generate Ontology Scores (ssGSEA)
The API also requires the `ontologyssGSEA.csv` file to serve the pathway visualization endpoints. This file is generated via the **`ETL.qmd`** document.

1. Open `ETL.qmd` (e.g., in RStudio, VSCode, or via the command line).
2. Ensure you have the `dados_otimizados.rds` file from the previous step.
3. Run the **"Ontology create dataset"** code chunk inside the document. 

This will load the optimized data, fetch gene sets from `msigdbr`, run the GSVA calculations, and output the final `ontologyssGSEA.csv` file to your root directory.

## 🏃 Running the Application

Once the data is prepared, you can launch the app:

```R
shiny::runApp()
```

## 🐳 Docker Support

The full application stack (React Frontend + Plumber API) can be easily served via Docker Compose.

### Using Docker Compose (Recommended)

The easiest way to build and run the complete application is using Docker Compose:

```bash
docker-compose up --build -d
```

This will automatically start both services:
- **React Frontend**: `http://localhost:5173`
- **Plumber API**: `http://localhost:33857` *(Swagger Docs available at `http://localhost:33857/__docs__/`)*

### Live Reloading during Development

Docker Compose is configured to automatically sync your code changes during development without needing full container rebuilds. 

Instead of `docker-compose up`, use the `watch` command:

```bash
docker compose watch
```

With this running:
- Any edits saved to `plumber.R` will instantly sync and restart the API container.
- Any edits saved in the `front/` directory will instantly trigger Vite's Hot Module Replacement (HMR) and update the frontend live in your browser.

### Using Docker Manually (API Only)

If you wish to run *only* the backend API, you can build and run it directly:

```bash
# Build the API image
docker build -t devbrain-markers-api -f Dockerfile.api .

# Run the API container
docker run -p 33857:33857 devbrain-markers-api
```

## 📊 Data Source

Data is sourced from the [BrainSpan Atlas of the Developing Human Brain](https://www.brainspan.org/), specifically the RNA-Seq RPKM values averaged to genes.

## 🤝 Acknowledgements

Made by the **[MCB Lab at UFRN](https://mcblab.github.io/)**.

If you use this tool in your research, please cite the BrainSpan Atlas and the `ggseg` package.
