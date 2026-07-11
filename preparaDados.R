library(vroom)
library(dplyr)
library(tibble)
library(msigdbr)
library(tidyr) 
library(GSVA)   
library(readr)

#Carregando matrizes de dados do BrainSpan
rows <- vroom("rows_metadata.csv")
columns <- vroom("columns_metadata.csv")
mapeamento_completo <- vroom("mapeamento_regioes.csv")

counts_raw <- vroom(
  "expression_matrix.csv",
  col_names = FALSE
) %>%
  dplyr::select(-1)

matriz_expressao <- as.matrix(counts_raw)
colnames(matriz_expressao) <- columns$column_num
rownames(matriz_expressao) <- rows$gene_symbol

genes_unicos <- !duplicated(rownames(matriz_expressao))
matriz_expressao <- matriz_expressao[genes_unicos, ]

matriz_final <- log2(matriz_expressao + 1)

mapping_df <- mapeamento_completo %>%
  select(region, structure_name) %>%
  as.data.frame()

mapping_macro <- mapeamento_completo %>%
  select(region, macro_region = macro_regions) %>%
  as.data.frame()

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

# Renomeando as colulas
columns_renomeado <- columns %>%
  rename(structure_original = structure_name)

columns_enriquecidas <- columns_renomeado %>%
  left_join(mapping_df, by = c("structure_original" = "structure_name")) %>%
  left_join(age_df, by = "age") %>%
  left_join(mapping_macro, by = "region") %>%
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
    structure_mapped,
    macro_region
  )

mapeamento_rapido <- columns_enriquecidas %>%
  select(column_num, broad_age, structure_mapped, macro_region)

# Salvar os dados otimizados com todas as informações de mapeamento
dados_app <- list(
  expression_matrix = matriz_final,
  col_meta = columns_enriquecidas,
  gene_list = rownames(matriz_final),
  mapping_info = list(
    region_to_structure = mapping_df,
    region_to_macro_structure = mapping_macro,
    age_groups = age_df,
    quick_lookup = mapeamento_rapido
  )
)

#Criando dados para as Ontologias
matrix_dados <- as.matrix(dados_app$expression_matrix)

h_df <- msigdbr(species = "Homo sapiens", category = "H")
biocarta_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
kegg_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY")
reactome_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
geneont_df <- msigdbr(species = "Homo sapiens", collection = "C5")

all_genesets_df <- bind_rows(h_df, biocarta_df, kegg_df, reactome_df, geneont_df)

genesets_list <- split(x = all_genesets_df$gene_symbol, f = all_genesets_df$gs_name)

# ssGSEA Score
ssgsea_param <- ssgseaParam(exprData = matrix_dados, geneSets = genesets_list)
ssgsea_results <- gsva(ssgsea_param)

tabela_final_ssgsea <- as.data.frame(ssgsea_results) %>%
  rownames_to_column(var = "GeneSet") %>%
  pivot_longer(
    cols = -GeneSet,          
    names_to = "Amostra",
    values_to = "Score_ssGSEA"
  )


saveRDS(dados_app, "dados_otimizados.rds")
saveRDS(genesets_list, "genesets_list.rds")
write_csv(tabela_final_ssgsea, "ontologyssGSEA.csv")