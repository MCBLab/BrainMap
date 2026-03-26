library(vroom)
library(dplyr)
library(tibble)


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

genes_unicos <- !duplicated(rownames(matriz_expressao))
matriz_expressao <- matriz_expressao[genes_unicos, ]
matriz_final <- log2(matriz_expressao + 1)

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

# Mapeamento de idades
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

# Renomear a coluna structure_name para evitar conflito
columns_renomeado <- columns %>%
  rename(structure_original = structure_name)

columns_enriquecidas <- columns_renomeado %>%
  left_join(mapping_df, by = c("structure_original" = "structure_name")) %>%
  left_join(age_df, by = "age") %>%
  mutate(
    structure_mapped = region
  ) %>%
  select(
    column_num,
    donor_id,
    donor_name,
    age,
    broad_age,
    structure_original,
    structure_mapped
  )

# Criar uma versão simplificada apenas com o necessário para joins rápidos
mapeamento_rapido <- columns_enriquecidas %>%
  select(column_num, broad_age, structure_mapped)

# Salvar os dados otimizados com todas as informações de mapeamento
dados_app <- list(
  expression_matrix = matriz_final,
  col_meta = columns_enriquecidas, # Agora com os mapeamentos incluídos
  gene_list = rownames(matriz_final),
  mapping_info = list(
    region_to_structure = mapping_df,
    age_groups = age_df,
    quick_lookup = mapeamento_rapido
  )
)

# Salvar o arquivo RDS
saveRDS(dados_app, "dados_otimizados.rds")