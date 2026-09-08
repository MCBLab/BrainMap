# Converte um ontologyssGSEA.csv já existente nos .rds agregados que o
# plumber.R passou a carregar, sem precisar rodar o ssGSEA de novo.
# Uso pontual, só para quem já tem os dados gerados localmente:
#   Rscript agregaOntologias.R

library(vroom)
library(dplyr)

dados_app <- readRDS("dados_otimizados.rds")

# Amostra -> regiao vem de mapeamento_regioes.csv, nao de col_meta$structure_mapped.
# O .rds guarda o mapeamento que valia quando foi gerado; ler o CSV faz esta
# reagregacao valer para a versao atual sem precisar reconstruir o .rds antes.
mapeamento <- read.csv("mapeamento_regioes.csv", stringsAsFactors = FALSE)

meta_agg <- dados_app$col_meta %>%
  distinct(column_num, broad_age, structure_original) %>%
  inner_join(mapeamento %>% select(structure_name, region, macro_region = macro_regions),
             by = c("structure_original" = "structure_name"),
             relationship = "many-to-many") %>%
  filter(!is.na(region), !is.na(broad_age)) %>%
  distinct(column_num, region, broad_age, macro_region)

message("==> Lendo ontologyssGSEA.csv...")
scores_com_meta <- vroom(
  "ontologyssGSEA.csv",
  col_types = list(GeneSet = "c", Amostra = "d", Score_ssGSEA = "d")
) %>%
  inner_join(meta_agg, by = c("Amostra" = "column_num"), relationship = "many-to-many")

message("==> Agregando...")
ontologia_micro <- scores_com_meta %>%
  group_by(GeneSet, broad_age, region) %>%
  summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")

ontologia_macro <- scores_com_meta %>%
  filter(!is.na(macro_region)) %>%
  group_by(GeneSet, broad_age, macro_region) %>%
  summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")

saveRDS(ontologia_micro, "ontologia_micro.rds")
saveRDS(ontologia_macro, "ontologia_macro.rds")

message(sprintf(
  "==> Pronto: ontologia_micro.rds (%d linhas), ontologia_macro.rds (%d linhas)",
  nrow(ontologia_micro), nrow(ontologia_macro)
))
