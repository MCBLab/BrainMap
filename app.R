#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseg)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(msigdbr)
library(scuttle)
library(tibble)
library(shinythemes)
library(svglite)
library(vroom)
# library(RSQLite)

# # Connect to the SQLite database
# con <- dbConnect(RSQLite::SQLite(), "brain_data.sqlite")
#
# # Get gene list from the database
# gene_list <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM rows_metadata")$gene_symbol
# Get gene list from the CSV file

setwd("~/Documentos/BrainSpan/")
rows <- vroom::vroom("genes_matrix_csv/rows_metadata.csv")
columns <- vroom::vroom("genes_matrix_csv/columns_metadata.csv")
gene_list <- unique(rows$gene_symbol)

# Read and assemble the expression data
counts <- vroom::vroom(
  "genes_matrix_csv/expression_matrix.csv",
  col_names = FALSE
) %>%
  dplyr::select(-1)
colnames(counts) <- columns$column_num
counts <- cbind(rows, counts)


# Mapping dataframes
input_values <- c(
  NA,
  "bankssts",
  "caudal middle frontal",
  "fusiform",
  "inferior parietal",
  "inferior temporal",
  "lateral occipital",
  "lateral orbitofrontal",
  "middle temporal",
  "pars opercularis",
  "pars orbitalis",
  "pars triangularis",
  "postcentral",
  "precentral",
  "rostral middle frontal",
  "superior frontal",
  "superior parietal",
  "superior temporal",
  "supramarginal",
  "temporal pole",
  "transverse temporal",
  "insula",
  "caudal anterior cingulate",
  "corpus callosum",
  "cuneus",
  "entorhinal",
  "isthmus cingulate",
  "lingual",
  "medial orbitofrontal",
  "parahippocampal",
  "paracentral",
  "pericalcarine",
  "posterior cingulate",
  "precuneus",
  "rostral anterior cingulate",
  "frontal pole"
)

output_values <- c(
  NA,
  "occipital neocortex",
  "dorsolateral prefrontal cortex",
  "inferolateral temporal cortex (area TEv, area 20)",
  "posteroventral (inferior) parietal cortex",
  "temporal neocortex",
  "primary visual cortex (striate cortex, area V1/17)",
  "orbital frontal cortex",
  "posterior (caudal) superior temporal cortex (area 22c)",
  "ventrolateral prefrontal cortex",
  "orbital frontal cortex",
  "ventrolateral prefrontal cortex",
  "primary somatosensory cortex (area S1, areas 3,1,2)",
  "primary motor cortex (area M1, area 4)",
  "dorsolateral prefrontal cortex",
  "anterior (rostral) cingulate (medial prefrontal) cortex",
  "parietal neocortex",
  "temporal neocortex",
  "posteroventral (inferior) parietal cortex",
  "amygdaloid complex",
  "primary auditory cortex (core)",
  "striatum",
  "anterior (rostral) cingulate (medial prefrontal) cortex",
  "cerebellum",
  "primary visual cortex (striate cortex, area V1/17)",
  "hippocampus (hippocampal formation)",
  "mediodorsal nucleus of thalamus",
  "primary visual cortex (striate cortex, area V1/17)",
  "orbital frontal cortex",
  "hippocampus (hippocampal formation)",
  "primary motor-sensory cortex (samples)",
  "primary visual cortex (striate cortex, area V1/17)",
  "mediodorsal nucleus of thalamus",
  "parietal neocortex",
  "anterior (rostral) cingulate (medial prefrontal) cortex",
  "frontal pole"
)

mapping_df <- data.frame(
  region = input_values,
  structure_name = output_values,
  stringsAsFactors = FALSE
)

ages <- c(
  "8 pcw",
  "9 pcw",
  "12 pcw",
  "13 pcw",
  "16 pcw",
  "17 pcw",
  "19 pcw",
  "21 pcw",
  "24 pcw",
  "25 pcw",
  "26 pcw",
  "35 pcw",
  "37 pcw",
  "4 mos",
  "10 mos",
  "1 yrs",
  "2 yrs",
  "3 yrs",
  "4 yrs",
  "8 yrs",
  "11 yrs",
  "13 yrs",
  "15 yrs",
  "18 yrs",
  "19 yrs",
  "21 yrs",
  "23 yrs",
  "30 yrs",
  "36 yrs",
  "37 yrs",
  "40 yrs"
)

age_mapping <- c(
  rep("1st trimester (n = 5)", 3),
  rep("2nd trimester (n = 10)", 5),
  rep("3rd trimester (n = 5)", 5),
  rep("Infant (n = 8)", 10),
  rep("Adult (n = 14)", 8)
)

age_df <- data.frame(
  age = ages,
  broad_age = age_mapping,
  stringsAsFactors = FALSE
)

#GSVA matrix
filtragem <- counts %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "") %>%
  column_to_rownames(var = "gene_symbol")

X <- filtragem %>%
  dplyr::select(-row_num, -gene_id, -ensembl_gene_id, -entrez_id) %>%
  as.matrix()

X <- log2(X + 1)

# UI
ui <- navbarPage(
  title = "Gene Expression in Human Developmental Brain",
  theme = shinytheme("cerulean"),

  header = tags$head(
    tags$link(
      rel = "stylesheet",
      type = "text/css",
      href = "https://fonts.googleapis.com/css?family=Lato"
    ),
    tags$style(HTML(
      "
            .navbar-brand {
                font-family: 'Lato', sans-serif;
                font-weight: 500;
            }
        "
    ))
  ),

  tabPanel(
    "App",
    sidebarLayout(
      sidebarPanel(
        helpText(
          "Start typing a HGNC gene name to see its expression across developmental stages"
        ),
        selectizeInput(
          "gene",
          "Select Gene:",
          choices = NULL,
          selected = "GFAP"
        ),
        downloadButton("downloadSVG", "Download SVG"),
        downloadButton("downloadTIFF", "Download TIFF")
      ),
      mainPanel(
        plotOutput("brainPlot", height = "1200px"),
        style = "padding-bottom: 70px;"
      )
    )
  ),

  tabPanel(
    "About",
    fluidPage(
      fluidRow(
        column(
          12,
          h3("About This App"),
          p(
            "This application visualizes gene expression data from the Human Developmental Brain RNA-Seq dataset."
          ),
          p(
            "The data is sourced from the",
            a(
              "BrainSpan Atlas of the Developing Human Brain",
              target = "_blank",
              href = "https://www.brainspan.org/"
            ),
            "."
          ),
          h4("Brain Regions Reference"),
          plotOutput("dkPlot", height = "500px")
        )
      )
    )
  ),

  tabPanel(
    "GSVA-Enrichr Mode",
    sidebarLayout(
      sidebarPanel(
        helpText(
          "Cole uma lista de genes para mapear a atividade dessa assinatura no cérebro."
        ),
        textAreaInput(
          "gene_input_list",
          "Lista de Genes:",
          height = "150px",
          placeholder = "TP53\nEGFR\nVEGFA"
        ),
        actionButton(
          "run_enrichment",
          "Mapear Assinatura",
          class = "btn-primary"
        ),
        hr(),
      ),
      mainPanel(
        plotOutput("gsvaPlot", height = "1200px")
      )
    )
  )
)

# Footer
footer = tags$footer(
  style = "position: fixed; bottom: 0; width: 100%; height: 50px; background-color: #f5f5f5; text-align: center; padding: 15px;",
  a(
    href = "https://mcblab.github.io/",
    target = "_blank",
    "Made by MCB Lab at UFRN"
  )
)

# Server
server <- function(input, output, session) {
  updateSelectizeInput(session, "gene", choices = gene_list, server = TRUE)

  brainPlotObject <- reactive({
    req(input$gene)
    gene <- input$gene

    # The `counts` data frame now holds all the data
    gene_data <- counts[counts$gene_symbol == gene, ]

    if (nrow(gene_data) == 0) {
      return(
        ggplot() +
          annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "Data not available for this gene",
            size = 8
          ) +
          theme_void()
      )
    }

    # Reshape the data to a long format
    # The columns to pivot are the ones that are not in rows_metadata
    row_metadata_cols <- colnames(rows)
    gene_data_long <- gene_data %>%
      pivot_longer(
        cols = -all_of(row_metadata_cols),
        names_to = "column_num",
        values_to = "expression_value"
      ) %>%
      mutate(column_num = as.integer(column_num)) %>%
      dplyr::select(column_num, expression_value)

    # Join with metadata
    plot_data <- gene_data_long %>%
      left_join(columns, by = "column_num") %>%
      left_join(mapping_df, by = "structure_name") %>%
      left_join(age_df, by = "age")

    plot_data %>%
      group_by(broad_age) %>%
      ggseg(
        .data = .,
        hemisphere = "left",
        colour = "black",
        mapping = aes(fill = log2(expression_value + 1))
      ) +
      scale_fill_gradientn(
        colours = c("royalblue", "firebrick", "goldenrod"),
        na.value = "white"
      ) +
      labs(
        fill = "Log2(RPKM + 1)",
      ) +
      facet_wrap(
        ~ factor(broad_age, unique(age_df$broad_age)),
        nrow = 5
      ) +
      theme_void() +
      theme(
        strip.text = element_text(size = 18),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA),
        legend.text = element_text(size = 12)
      )
  })

  output$brainPlot <- renderPlot({
    brainPlotObject()
  })

  output$downloadSVG <- downloadHandler(
    filename = function() {
      paste0("brain_plot-", input$gene, ".svg")
    },
    content = function(file) {
      ggsave(
        file,
        plot = brainPlotObject(),
        device = "svg",
        width = 10,
        height = 12
      )
    }
  )

  output$downloadTIFF <- downloadHandler(
    filename = function() {
      paste0("brain_plot-", input$gene, ".tiff")
    },
    content = function(file) {
      ggsave(
        file,
        plot = brainPlotObject(),
        device = "tiff",
        width = 10,
        height = 12,
        dpi = 300
      )
    }
  )

  output$dkPlot <- renderPlot({
    plot(dk) +
      theme_void() +
      theme(
        strip.text = element_text(size = 18),
        legend.position = "bottom",
        title = element_blank(),
        legend.text = element_text(size = 12)
      )
  })

  #GVSA
  # Calculating GVSA
  gsva_result_reactive <- eventReactive(input$run_enrichment, {
    req(input$gene_input_list)

    # Filtragem do texto de entrada, deixando o vetor de genes limpo
    user_genes <- unlist(strsplit(input$gene_input_list, "\n"))
    user_genes <- trimws(user_genes) # Remove espaços extras
    user_genes <- user_genes[user_genes != ""] # Remove linhas vazias
    user_genes <- unique(user_genes)

    # Validar se os genes existem na matriz X
    genes_validos <- intersect(user_genes, rownames(X))

    if (length(genes_validos) < 3) {
      showNotification(
        "Atenção: Menos de 3 genes da lista foram encontrados nos dados. O GSVA pode falhar ou ser impreciso.",
        type = "warning"
      )
    }

    # Lista de Gene Sets
    custom_gene_set <- list("Custom_Signature" = genes_validos)

    # Configurar parametros e rodar GSVA
    params_custom <- gsvaParam(
      exprData = X,
      geneSets = custom_gene_set,
      kcdf = "Gaussian"
    ) # Gaussian -> dados em log2

    gsva_score_matrix <- gsva(params_custom)

    return(gsva_score_matrix)
  })

  # Renderizando o Plot do GSVA
  output$gsvaPlot <- renderPlot({
    req(gsva_result_reactive())

    # Pegar a matriz de scores
    score_matrix <- gsva_result_reactive()

    plot_data_gsva <- data.frame(
      column_num = colnames(score_matrix),
      gsva_score = as.numeric(score_matrix[1, ])
    ) %>%
      mutate(column_num = as.integer(column_num)) %>%
      left_join(columns, by = "column_num") %>%
      left_join(mapping_df, by = "structure_name") %>%
      left_join(age_df, by = "age")

    plot_data_gsva %>%
      group_by(broad_age) %>%
      ggseg(
        .data = .,
        hemisphere = "left",
        colour = "black",
        mapping = aes(fill = gsva_score)
      ) +
      scale_fill_gradient2(
        low = "royalblue",
        mid = "firebrick",
        high = "goldenrod",
        midpoint = 0,
        name = "GSVA Score"
      ) +
      labs(
        title = "Enriquecimento da Assinatura Personalizada nas Áreas Cerebrais",
        subtitle = paste(
          "Baseado em",
          length(intersect(
            rownames(X),
            unique(trimws(unlist(strsplit(input$gene_input_list, "\n"))))
          )),
          "genes encontrados"
        )
      ) +
      facet_wrap(
        ~ factor(broad_age, unique(age_df$broad_age)),
        nrow = 5
      ) +
      theme_void() +
      theme(
        strip.text = element_text(size = 18),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA),
        legend.text = element_text(size = 12)
      )
  })

  # session$onSessionEnded(function() {
  #     dbDisconnect(con)
  # })
}
# Run the application
shinyApp(ui = ui, server = server)
