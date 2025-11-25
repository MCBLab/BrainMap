library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseg)
library(shinythemes)
library(svglite)
library(RSQLite)

# Connect to the SQLite database
con <- dbConnect(RSQLite::SQLite(), "brain_data.sqlite")

# Get gene list from the database
gene_list <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM rows_metadata")$gene_symbol

# Mapping dataframes
input_values <- c(NA, "bankssts", "caudal middle frontal", "fusiform", "inferior parietal", 
                  "inferior temporal", "lateral occipital", "lateral orbitofrontal", 
                  "middle temporal", "pars opercularis", "pars orbitalis", "pars triangularis", 
                  "postcentral", "precentral", "rostral middle frontal", "superior frontal", 
                  "superior parietal", "superior temporal", "supramarginal", "temporal pole", 
                  "transverse temporal", "insula", "caudal anterior cingulate", "corpus callosum", 
                  "cuneus", "entorhinal", "isthmus cingulate", "lingual", "medial orbitofrontal", 
                  "parahippocampal", "paracentral", "pericalcarine", "posterior cingulate", 
                  "precuneus", "rostral anterior cingulate", "frontal pole")

output_values <- c(NA,"occipital neocortex", "dorsolateral prefrontal cortex", 
                   "inferolateral temporal cortex (area TEv, area 20)", 
                   "posteroventral (inferior) parietal cortex", "temporal neocortex", 
                   "primary visual cortex (striate cortex, area V1/17)", "orbital frontal cortex", 
                   "posterior (caudal) superior temporal cortex (area 22c)", 
                   "ventrolateral prefrontal cortex", "orbital frontal cortex", 
                   "ventrolateral prefrontal cortex", "primary somatosensory cortex (area S1, areas 3,1,2)", 
                   "primary motor cortex (area M1, area 4)", "dorsolateral prefrontal cortex", 
                   "anterior (rostral) cingulate (medial prefrontal) cortex", "parietal neocortex", 
                   "temporal neocortex", "posteroventral (inferior) parietal cortex", 
                   "amygdaloid complex", "primary auditory cortex (core)", "striatum", 
                   "anterior (rostral) cingulate (medial prefrontal) cortex", "cerebellum", 
                   "primary visual cortex (striate cortex, area V1/17)", "hippocampus (hippocampal formation)", 
                   "mediodorsal nucleus of thalamus", "primary visual cortex (striate cortex, area V1/17)", 
                   "orbital frontal cortex", "hippocampus (hippocampal formation)", 
                   "primary motor-sensory cortex (samples)", "primary visual cortex (striate cortex, area V1/17)", 
                   "mediodorsal nucleus of thalamus", "parietal neocortex", 
                   "anterior (rostral) cingulate (medial prefrontal) cortex", "frontal pole")

mapping_df <- data.frame(region = input_values, structure_name = output_values, stringsAsFactors = FALSE)

ages <- c("8 pcw", "9 pcw", "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", 
          "24 pcw", "25 pcw", "26 pcw", "35 pcw", "37 pcw", "4 mos", "10 mos", "1 yrs", 
          "2 yrs", "3 yrs", "4 yrs", "8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", 
          "19 yrs", "21 yrs", "23 yrs", "30 yrs", "36 yrs", "37 yrs", "40 yrs")

age_mapping <- c(
  rep("1st trimester (n = 5)", 3),
  rep("2nd trimester (n = 10)", 5),
  rep("3rd trimester (n = 5)", 5),
  rep("Infant (n = 8)", 10),
  rep("Adult (n = 14)", 8)
)

age_df <- data.frame(age = ages, broad_age = age_mapping, stringsAsFactors = FALSE)


# UI
ui <- navbarPage(
    title = "Gene Expression in Human Developmental Brain",
    theme = shinytheme("cerulean"),
    
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "https://fonts.googleapis.com/css?family=Lato"),
        tags$style(HTML("
            .navbar-brand {
                font-family: 'Lato', sans-serif;
                font-weight: 500;
            }
        "))
    ),
    
    tabPanel("App",
        sidebarLayout(
            sidebarPanel(
                helpText("Start typing a HGNC gene name to see its expression across developmental stages"),
                selectizeInput("gene", "Select Gene:", choices = NULL, selected = "GFAP"),
                downloadButton("downloadSVG", "Download SVG"),
                downloadButton("downloadTIFF", "Download TIFF")
            ),
            mainPanel(
                plotOutput("brainPlot", height = "1200px")
            )
        )
    ),
    tabPanel("About", 
        mainPanel(
            h3("About This App"),
            p("This application visualizes gene expression data from the Human Developmental Brain RNA-Seq dataset."),
            p("The data is sourced from the",
                a("BrainSpan Atlas of the Developing Human Brain",  target = "_blank", href = "https://www.brainspan.org/"),
                "."),
            h4("Brain Regions Reference"),
            plotOutput("dkPlot", height = "500px")
        )
    ),
    
    # Footer
    tags$footer(
        style = "position: fixed; bottom: 0; width: 100%; height: 50px; background-color: #f5f5f5; text-align: center; padding: 15px;",
        a(href = "https://mcblab.github.io/", target = "_blank", "Made by MCB Lab at UFRN")
    )
)

# Server
server <- function(input, output, session) {
    updateSelectizeInput(session, "gene", choices = gene_list, server = TRUE)

    brainPlotObject <- reactive({
        req(input$gene)
        gene <- input$gene

        query <- paste0(
            "SELECT * FROM expression WHERE gene_symbol = '",
            gene,
            "'"
        )
        gene_data <- dbGetQuery(con, query)
        
        if (nrow(gene_data) == 0) {
            return(
                ggplot() +
                annotate("text", x = 0.5, y = 0.5, label = "Data not available for this gene", size = 8) +
                theme_void()
            )
        }

        # We need to read the columns_metadata table to join with the expression data
        columns_metadata <- dbReadTable(con, "columns_metadata")

        gene_data %>%
            merge(
                .,
                columns_metadata,
                by.x = "column_num",
                by.y = "column_num"
            ) %>%
            merge(mapping_df, ., by = "structure_name") %>%
            merge(age_df, ., by = "age") %>%
            group_by(broad_age) %>%
            ggseg(
                .data = .,
                hemisphere = "left",
                colour = "black",
                mapping = aes(fill = log2(expression_value + 1))
            ) +
            scale_fill_gradientn(
                colours = c("royalblue","firebrick","goldenrod"),
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
            ggsave(file, plot = brainPlotObject(), device = "svg", width = 10, height = 12)
        }
    )

    output$downloadTIFF <- downloadHandler(
        filename = function() {
            paste0("brain_plot-", input$gene, ".tiff")
        },
        content = function(file) {
            ggsave(file, plot = brainPlotObject(), device = "tiff", width = 10, height = 12, dpi = 300)
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

    session$onSessionEnded(function() {
        dbDisconnect(con)
    })
}
# Run the application 
shinyApp(ui = ui, server = server)
