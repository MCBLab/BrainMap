
library(dplyr)
library(tidyr)
library(gtools)
library(ggplot2)
library(ggseg)

# Path to working directory
setwd("/path/to/your/working/directory")
gene = "GFAP"


# Downloaded from:
# https://www.brainspan.org/api/v2/well_known_file_download/267666525
# 
columns <- read.csv("genes_matrix_csv/columns_metadata.csv", header = T)
counts <- read.csv("genes_matrix_csv/expression_matrix.csv", header = F)[,-1]
rows <- read.csv("genes_matrix_csv/rows_metadata.csv", header = T)

colnames(counts) = columns$column_num
counts = cbind(rows, counts)

# Vetor de entrada
input_values <- c(NA, "bankssts", "caudal middle frontal", "fusiform", "inferior parietal", 
                  "inferior temporal", "lateral occipital", "lateral orbitofrontal", 
                  "middle temporal", "pars opercularis", "pars orbitalis", "pars triangularis", 
                  "postcentral", "precentral", "rostral middle frontal", "superior frontal", 
                  "superior parietal", "superior temporal", "supramarginal", "temporal pole", 
                  "transverse temporal", "insula", "caudal anterior cingulate", "corpus callosum", 
                  "cuneus", "entorhinal", "isthmus cingulate", "lingual", "medial orbitofrontal", 
                  "parahippocampal", "paracentral", "pericalcarine", "posterior cingulate", 
                  "precuneus", "rostral anterior cingulate", "frontal pole")

# Vetor de saída corrigido (BrainSpan regions)
output_values <- c(
  NA,
  "posterior (caudal) superior temporal cortex (area 22c)",  # bankssts - parte do sulco temporal superior
  "dorsolateral prefrontal cortex",  # caudal middle frontal
  "inferolateral temporal cortex (area TEv, area 20)",  # fusiform
  "posteroventral (inferior) parietal cortex",  # inferior parietal
  "inferolateral temporal cortex (area TEv, area 20)",  # inferior temporal
  "occipital neocortex",  # lateral occipital - V2/V3, não V1
  "orbital frontal cortex",  # lateral orbitofrontal
  "inferolateral temporal cortex (area TEv, area 20)",  # middle temporal
  "ventrolateral prefrontal cortex",  # pars opercularis
  "ventrolateral prefrontal cortex",  # pars orbitalis - mais VLPFC que OFC
  "ventrolateral prefrontal cortex",  # pars triangularis
  "primary somatosensory cortex (area S1, areas 3,1,2)",  # postcentral
  "primary motor cortex (area M1, area 4)",  # precentral
  "dorsolateral prefrontal cortex",  # rostral middle frontal
  "dorsolateral prefrontal cortex",  # superior frontal
  "posteroventral (inferior) parietal cortex",  # superior parietal
  "posterior (caudal) superior temporal cortex (area 22c)",  # superior temporal
  "posteroventral (inferior) parietal cortex",  # supramarginal
  "inferolateral temporal cortex (area TEv, area 20)",  # temporal pole
  "primary auditory cortex (core)",  # transverse temporal
  "ventrolateral prefrontal cortex",  # insula - insula é única, uso VLPFC como próximo
  "anterior (rostral) cingulate (medial prefrontal) cortex",  # caudal anterior cingulate
  "orbital frontal cortex",  # corpus callosum - melhor usar uma região cortical genérica
  "primary visual cortex (striate cortex, area V1/17)",  # cuneus
  "inferolateral temporal cortex (area TEv, area 20)",  # entorhinal - próximo ao temporal
  "posterior (caudal) superior temporal cortex (area 22c)",  # isthmus cingulate
  "primary visual cortex (striate cortex, area V1/17)",  # lingual
  "orbital frontal cortex",  # medial orbitofrontal
  "inferolateral temporal cortex (area TEv, area 20)",  # parahippocampal - próximo
  "primary motor cortex (area M1, area 4)",  # paracentral - motor/sensory
  "primary visual cortex (striate cortex, area V1/17)",  # pericalcarine
  "posterior (caudal) superior temporal cortex (area 22c)",  # posterior cingulate
  "posteroventral (inferior) parietal cortex",  # precuneus - parietal medial
  "anterior (rostral) cingulate (medial prefrontal) cortex",  # rostral anterior cingulate
  "dorsolateral prefrontal cortex"  # frontal pole
)

# Criando um dataframe para visualizar o mapeamento
mapping_df <- data.frame(region = input_values, structure_name = output_values, stringsAsFactors = FALSE)

# Criando o vetor de idades
ages <- c("8 pcw", "9 pcw", "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", 
          "24 pcw", "25 pcw", "26 pcw", "35 pcw", "37 pcw", "4 mos", "10 mos", "1 yrs", 
          "2 yrs", "3 yrs", "4 yrs", "8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", 
          "19 yrs", "21 yrs", "23 yrs", "30 yrs", "36 yrs", "37 yrs", "40 yrs")


# How many donors per group
merge(columns, age_df, by = "age") %>%
  group_by(broad_age, donor_id) %>%
  tally() %>% ungroup() %>%
  group_by(broad_age) %>% tally()

# Criando um mapeamento das idades para categorias
age_mapping <- c(
  rep("1st trimester (n = 5)", 3),   # 8 pcw - 12 pcw
  rep("2nd trimester (n = 10)", 5),   # 13 pcw - 24 pcw
  rep("3rd trimester (n = 5)", 5),   # 25 pcw - 37 pcw
  rep("Infant (n = 8)", 10),          # 4 mos - 15 yrs
  rep("Adult (n = 14)", 8)           # 18 yrs - 40 yrs
)

# Criando um dataframe com o mapeamento
age_df <- data.frame(age = ages, broad_age = age_mapping, stringsAsFactors = FALSE)


labels_vector <- columns$age %>% unique()

# Criando o labeller
age_labeller <- function(variable, value) {
  return(labels_vector[value])
}

png(paste0(gene,"_dev_brain_RNAseq.png"), width = 2800*1.4, height = 800*1.4, res = 250)

counts %>% filter(gene_symbol == gene) %>% select(-c(gene_id, ensembl_gene_id, gene_symbol, entrez_id)) %>%
  pivot_longer(!row_num) %>% merge(., columns, by.x = "name", by.y = "column_num") %>%
  merge(mapping_df, ., by = "structure_name") %>%
  merge(age_df, ., by = "age") %>%
  group_by(broad_age) %>%
  ggseg(.data=., hemisphere = "left", colour = "black", mapping = aes(fill = log2(value+1))) +
  scale_fill_gradientn(colours = c("royalblue","firebrick","goldenrod"),na.value="white")  +
  labs(fill = "Log2(RPKM + 1)", title = paste0("Human Developmental brain RNAseq - ", gene)) +
  facet_wrap(~factor(broad_age, unique(age_df$broad_age)), ncol = 5) + theme_void() +
  theme(strip.text = element_text(size = 18), legend.position = "bottom", title = element_text(size = 24), plot.title = element_text(hjust = .5, vjust = 1))

dev.off()  

svg(paste0(gene, "_dev_brain_RNAseq.svg"), width = 14, height = 4)

counts %>% as.data.frame() %>% filter(gene_symbol == gene) %>% select(-c(gene_id, ensembl_gene_id, gene_symbol, entrez_id)) %>%
  pivot_longer(!row_num) %>% merge(., columns, by.x = "name", by.y = "column_num") %>%
  merge(mapping_df, ., by = "structure_name") %>%
  merge(age_df, ., by = "age") %>%
  group_by(broad_age) %>%
  ggseg(.data=., hemisphere = "left", colour = "black", mapping = aes(fill = log2(value+1))) +
  scale_fill_gradientn(colours = c("royalblue","firebrick","goldenrod"),na.value="white")  +
  labs(fill = "Log2(RPKM + 1)", title = paste0("Human Developmental brain RNAseq - ", gene)) +
  facet_wrap(~factor(broad_age, unique(age_df$broad_age)), ncol = 5) + theme_void() +
  theme(strip.text = element_text(size = 18), legend.position = "bottom", title = element_text(size = 24), plot.title = element_text(hjust = .5, vjust = 1))

dev.off()  

pdf(paste0(gene,"_dev_brain_RNAseq.pdf"), width = 14, height = 4)

counts %>% as.data.frame() %>% filter(gene_symbol == gene) %>% select(-c(gene_id, ensembl_gene_id, gene_symbol, entrez_id)) %>%
  pivot_longer(!row_num) %>% merge(., columns, by.x = "name", by.y = "column_num") %>%
  merge(mapping_df, ., by = "structure_name") %>%
  merge(age_df, ., by = "age") %>%
  group_by(broad_age) %>%
  ggseg(.data=., hemisphere = "left", colour = "black", mapping = aes(fill = log2(value+1))) +
  scale_fill_gradientn(colours = c("royalblue","firebrick","goldenrod"),na.value="white")  +
  labs(fill = "Log2(RPKM + 1)", title = paste0("Human Developmental brain RNAseq - ", gene)) +
  facet_wrap(~factor(broad_age, unique(age_df$broad_age)), ncol = 5) + theme_void() +
  theme(strip.text = element_text(size = 18), legend.position = "bottom", title = element_text(size = 24), plot.title = element_text(hjust = .5, vjust = 1))

dev.off()  