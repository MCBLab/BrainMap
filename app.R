#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseg)
library(GSVA)
library(shinythemes)
library(svglite)
library(shinycssloaders)

dados <- readRDS("dados_otimizados.rds")

X <- dados$expression_matrix
columns <- dados$col_meta
gene_list <- dados$gene_list

mapping_df <- dados$mapping_info$region_to_structure
age_df <- dados$mapping_info$age_groups

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
        withSpinner(plotOutput("brainPlot", height = "550px"), type = 6),
        style = "padding-bottom: 70px;"
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
          placeholder = "PAX6\nEGFR\nVEGFA"
        ),
        actionButton(
          "run_enrichment",
          "Mapear Assinatura",
          class = "btn-primary"
        ),
        hr(),
      ),
      mainPanel(
        withSpinner(plotOutput("gsvaPlot", height = "550px"), type = 6),
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
  )
)

# Footer
ui <- tagList(
  ui,
  tags$footer(
    style = "position: fixed; bottom: 0; width: 100%; height: 50px; background-color: #f5f5f5; text-align: center; padding: 15px;",
    a(
      href = "https://mcblab.github.io/",
      target = "_blank",
      "Made by MCB Lab at UFRN"
    )
  )
)

# Server
server <- function(input, output, session) {

      updateSelectizeInput(
        session,
        "gene",
        choices = gene_list,
        selected = "GFAP",
        server = TRUE
      )

  
  brainPlotObject <- reactive({
    req(input$gene)
    gene <- input$gene
    
    if (!gene %in% rownames(X)) {
      return(
        ggplot() +
          annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "Data not available",
            size = 8
          ) +
          theme_void()
      )
    }
    
    expression_values <- X[gene, ]
    
    # Criar dataframe leve apenas para plotagem
    plot_data <- data.frame(
      column_num = colnames(X),
      expression_value = as.numeric(expression_values)
    ) %>%
      mutate(column_num = as.integer(column_num)) %>%
      left_join(columns, by = "column_num") %>%
      rename(region = structure_mapped)
    
    plot_data %>%
      filter(!is.na(region)) %>%
      group_by(broad_age) %>%
      ggseg(
        .data = .,
        hemisphere = "left",
        colour = "black",
        mapping = aes(fill = expression_value)
      ) +
      scale_fill_gradient(
        low = "blue", high = "orange",
        na.value = "white"
      ) +
      labs(
        fill = "Log2(RPKM + 1)",
      ) +
      facet_wrap(
        ~ factor(broad_age, unique(age_df$broad_age)),
        nrow = 1,
        ncol = 5
      ) +
      theme_void() +
      theme(
        strip.text = element_text(
          size = 14,
          face = "bold",
          angle = 0,
          hjust = 0.5,
          margin = margin(b = 5, t = 2)
        ),
        strip.placement = "outside",
        strip.background = element_blank(),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA),
        legend.text = element_text(size = 12),
        panel.spacing = unit(0.3, "lines")
      )
  })
  
  output$brainPlot <- renderPlot(
    {
      brainPlotObject()
    },
    height = 450,
    width = 820
  )
  
  output$downloadSVG <- downloadHandler(
    filename = function() {
      paste0("brain_plot-", input$gene, ".svg")
    },
    content = function(file) {
      ggsave(
        file,
        plot = brainPlotObject(),
        device = "svg",
        width = 18,
        height = 5
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
        width = 18,
        height = 5,
        dpi = 300,
        compression = "lzw"
      )
    }
  )
  
  output$dkPlot <- renderPlot({
    plot(dk) +
      theme_void() +
      theme(
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        title = element_blank(),
        legend.text = element_text(size = 10)
      )
  })
  
  #GVSA
  # Calculating GVSA
  gsva_result_reactive <- eventReactive(input$run_enrichment, {
    req(input$gene_input_list)
    
    # Filtragem do texto de entrada, deixando o vetor de genes limpo
    user_genes <- unique(trimws(unlist(strsplit(input$gene_input_list, "\n")))) #trimws -> remove espaços extras
    user_genes <- user_genes[user_genes != ""] #remove linhas vazias
    
    matrix_dados_allen <- as.matrix(dados_app$expression_matrix)
    matrix_dados_allen[is.infinite(matrix_dados_allen)] <- 0
    
    matrix_dados_allen <- matrix_dados_allen[rowSums(matrix_dados_allen != 0) > 0, ]
    
    genes_validos <- intersect(gene_list, rownames(matrix_dados_allen))
    genes_ausentes <- setdiff(gene_list, rownames(matrix_dados_allen))
    
    if (length(genes_validos) == 0) {
      showNotification(
        "Erro: Nenhum dos genes da lista foi encontrado nos dados.",
        type = "error"
      )
      return(NULL)
    }
    
    if (length(genes_validos) < 5) {
      showNotification(
        "Atenção: Menos de 5 genes da lista foram encontrados nos dados. O GSVA pode falhar ou ser impreciso.",
        type = "warning"
      )
    }
  
    custom_gene_set <- list("Custom_Signature" = genes_validos)
    
    params_custom_ssgsea <- ssgseaParam(
      exprData = matrix_dados_allen,
      geneSets = custom_gene_set,
      normalize = TRUE
    )
    
    gsva_score_matrix_ssgsea <- gsva(params_custom_ssgsea)
    
    return(list(scores = gsva_score_matrix_ssgsea, missing = genes_ausentes))
  })
  
  # Renderizando o Plot do GSVA
  output$gsvaPlot <- renderPlot(
    {
      # Peganado a matriz de scores
      resultados <- gsva_result_reactive()
      req(resultados) # Garantir que resultados não é NULL
      
      score_matrix <- resultados$scores
      genes_ausentes <- resultados$missing
      
      tabela_final_ssgsea_allen <- as.data.frame(gsva_score_matrix_ssgsea) %>%
        rownames_to_column(var = "GeneSet") %>%
        pivot_longer(
          cols = -GeneSet,          
          names_to = "Amostra",
          values_to = "Score_ssGSEA"
        )
      
      #Condicional para a ausencia dos genes
      if (length(genes_ausentes) == 0) {
        texto_legenda <- "Todos os genes apresentados estão listados"
      } else {
        lista_genes <- paste(genes_ausentes, collapse = ", ")
        if (nchar(lista_genes) > 80) {
          lista_genes <- paste0(substr(lista_genes, 1, 80), "...")
        }
        texto_legenda <- paste(
          "Genes não listados na base de dados:",
          lista_genes
        )
      }
      
      dados_prontos <- tabela_final_gsva_allen %>%
        mutate(Amostra = as.numeric(Amostra)) %>%
        left_join(dados_app$col_meta, by = c("Amostra" = "column_num")) %>%
        rename(region = structure_mapped) %>%
        filter(
          GeneSet == "Custom_Signature", 
          !is.na(region),        
          !is.na(broad_age),
          !is.na(Score_GSVA)
        ) %>%
        group_by(broad_age, region) %>%
        summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
      
      plot_data_gsva %>%
        filter(!is.na(region)) %>%
        group_by(broad_age) %>%
        dados_prontos %>% 
        ggplot() +
        geom_brain(
          atlas = ggseg::dk(), 
          position = position_brain(c(
            "right lateral", "right medial")),
          mapping = aes(fill = Score_Medio_ssGSEA)
        ) +
        scale_fill_viridis_c(
          option = "viridis", 
          name = "ssGSEA Score"
        ) +
        labs(
          title = "Enriquecimento da Assinatura Personalizada nas Áreas Cerebrais",
          subtitle = texto_legenda
        ) +
        facet_wrap(
          ~ factor(broad_age, unique(age_df$broad_age)),
          nrow = 1,
          ncol = 5
        ) +
        theme_void() +
        theme(
          plot.title = element_text(
            size = 16,
            face = "bold",
            hjust = 0.5,
            vjust = 1,
            margin = margin(t = 15, b = 10),
            color = "black"
          ),
          plot.subtitle = element_text(
            size = 14,
            hjust = 0.5,
            vjust = 1,
            margin = margin(b = 20),
            color = "gray40",
            face = "italic"
          ),
          strip.text = element_text(
            size = 14,
            face = "bold",
            angle = 0,
            hjust = 0.5,
            margin = margin(b = 5, t = 2)
          ),
          strip.placement = "outside",
          strip.background = element_blank(),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = NA),
          legend.text = element_text(size = 12),
          panel.spacing = unit(0.3, "lines")
        )
    },
    height = 450,
    width = 820
  )
}

# Run the application
shinyApp(ui = ui, server = server)