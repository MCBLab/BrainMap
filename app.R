library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseg)
library(shinythemes) # Added shinythemes library

# Load data
columns <- read.csv("genes_matrix_csv/columns_metadata.csv", header = T)
counts <- read.csv("genes_matrix_csv/expression_matrix.csv", header = F)[,-1]
rows <- read.csv("genes_matrix_csv/rows_metadata.csv", header = T)

colnames(counts) = columns$column_num
counts = cbind(rows, counts)

# Get gene list
gene_list <- unique(rows$gene_symbol)

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
ui <- fluidPage(theme = shinytheme("flatly"), # Applied "flatly" theme
    titlePanel("Gene Expression in Human Developmental Brain"),
    sidebarLayout(
        sidebarPanel(
            selectizeInput("gene", "Select Gene:", choices = NULL, selected = "GFAP")
        ),
        mainPanel(
            plotOutput("brainPlot", height = "1200px") # Increased height for vertical facets
        )
    )
)

# Server
server <- function(input, output, session) {
    updateSelectizeInput(session, "gene", choices = gene_list, server = TRUE)

    output$brainPlot <- renderPlot({
        gene <- input$gene
        
        counts %>% filter(gene_symbol == gene) %>% select(-c(gene_id, ensembl_gene_id, gene_symbol, entrez_id)) %>%
          pivot_longer(!row_num) %>% merge(., columns, by.x = "name", by.y = "column_num") %>%
          merge(mapping_df, ., by = "structure_name") %>%
          merge(age_df, ., by = "age") %>%
          group_by(broad_age) %>%
          ggseg(.data=., hemisphere = "left", colour = "black", mapping = aes(fill = log2(value+1))) +
          scale_fill_gradientn(colours = c("royalblue","firebrick","goldenrod"),na.value="white")  +
          labs(fill = "Log2(RPKM + 1)", title = paste0("Human Developmental brain RNAseq - ", gene)) +
          facet_wrap(~factor(broad_age, unique(age_df$broad_age)), nrow = 5) + theme_void() + # Changed ncol to nrow
          theme(strip.text = element_text(size = 18), legend.position = "bottom", title = element_text(size = 24), plot.title = element_text(hjust = .5, vjust = 1))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
