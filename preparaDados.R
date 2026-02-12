#==============================================
#BrainSpan: Atlas of the Developing Human Brain
#==============================================

#This data set contains RNA-Seq RPKM (reads per kilobase per million; see the whitepaper at
#                                       www.brainspan.org) values averaged to genes.

#expression_matrix.csv -- the rows are genes and the columns samples; the first column is the row number
#rows_metadata.csv -- the genes are listed in the same order as the rows in expression_matrix.csv
#columns_metadata.csv -- the samples are listed in the same order as the columns in expression_matrix

library(vroom)
library(dplyr)
library(tibble)

setwd("~/AppShiny/") # Ajuste para o seu caminho

rows <- vroom("genes_matrix_csv/rows_metadata.csv")
columns <- vroom("genes_matrix_csv/columns_metadata.csv")

counts_raw <- vroom(
  "genes_matrix_csv/expression_matrix.csv",
  col_names = FALSE
) %>%
  dplyr::select(-1)

matriz_expressao <- as.matrix(counts_raw)

colnames(matriz_expressao) <- columns$column_num
rownames(matriz_expressao) <- rows$gene_symbol

# 3. Pré-calcular Log2 (Para não fazer isso a cada clique)
# Lidando com duplicatas de genes ou apenas removendo duplicatas para garantir rownames únicos
genes_unicos <- !duplicated(rownames(matriz_expressao))
matriz_expressao <- matriz_expressao[genes_unicos, ]
matriz_final <- log2(matriz_expressao + 1)

#media dos genes normalizados, caso < 1, é retirado do dataset
media_genes <- rowMeans(matriz_final, na.rm = T)
matriz_final <- matriz_final[media_genes > 1, ]

dados_app <- list(
  expression_matrix = matriz_final,
  col_meta = columns
)

View(dados_app$expression_matrix)

saveRDS(dados_app, "dados_otimizados.rds")


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


#verificação

cat("=== Comparação de estruturas Mapping vs Col_meta ===\n")

estruturas_no_dataset <- unique(dados_app$col_meta$structure_name)
estruturas_no_mapping <- na.omit(mapping_df$structure_name)

existem_no_dataset <- estruturas_no_mapping %in% estruturas_no_dataset

if (!all(existem_no_dataset)) {
  cat("\n❌ Estruturas FALTANTES no dataset:\n")
  faltantes <- estruturas_no_mapping[!existem_no_dataset]
  print(faltantes)
} else {
  cat("✅ TODAS as estruturas do mapping existem no dataset!\n")
}

# 2. Verificar se há estruturas no dataset que NÃO estão no mapping
cat("\n=== ESTRUTURAS EXTRAS NO DATASET ===\n")
extras_no_dataset <- setdiff(estruturas_no_dataset, estruturas_no_mapping)
cat(sprintf("Estruturas no dataset: %d\n", length(estruturas_no_dataset)))
cat(sprintf(
  "Estruturas que NÃO estão no mapping: %d\n",
  length(extras_no_dataset)
))

if (length(extras_no_dataset) > 0) {
  cat("Exemplos de estruturas extras:\n")
  print(head(extras_no_dataset, 10))
}
