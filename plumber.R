library(plumber)
library(dplyr)
library(ggplot2)
library(ggseg)

#Use o comando com a porta = 33857 para se vincular com o front!!
#plumb(file='Documentos/BrainMap/plumber.R')$run(port = 33857)

dados_app <- readRDS("~/Documentos/BrainMap/dados_otimizados.rds")
tabela_final_ssgsea <- read.csv("~/Documentos/BrainMap/ontologyssGSEA.csv")

matrix_dados <- dados_app$expression_matrix

meta_limpo <- dados_app$col_meta %>%
  rename(region = structure_mapped) %>%
  filter(!is.na(region), !is.na(broad_age)) %>%
  select(column_num, region, broad_age)

# -- OTIMIZAÇÃO 2: Cruzar a tabela das Ontologias --
# Como a tabela de ontologias é fixa, já fazemos o join e a limpeza de NAs aqui fora
tabela_base_ontologias <- tabela_final_ssgsea %>%
  mutate(Amostra = as.numeric(Amostra)) %>%
  inner_join(meta_limpo, by = c("Amostra" = "column_num"))

setgenes <- unique(dados_app$gene_list)
setontologies <- unique(tabela_final_ssgsea$GeneSet)

#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  plumber::forward()
}

#* @get /list_genes
function(){
  setgenes
}

#* @get /list_ontologies
function() {
  # Retorna um vetor único de GeneSets para o front-end montar a lista
  setontologies
}

#* Retorna o plot do cérebro com a expressão do gene selecionado
#* @param gene Nome do gene selecionado no front-end (ex: SOX10)
#* @serializer png list(width = 800, height = 600, res = 120)
#* @get /plot_brain
function(gene = "SOX10") { 
  
  if (!(gene %in% rownames(matrix_dados))) {
    stop("Erro: Gene não encontrado.") 
  }
  
  # Lógica otimizada: Apenas extrai a linha da matriz e junta com o meta_limpo
  dados_gene_plot <- data.frame(
    column_num = as.numeric(colnames(matrix_dados)),
    Expressao = as.numeric(matrix_dados[gene, ])
  ) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    group_by(broad_age, region) %>%
    summarise(Expressao_Media = mean(Expressao, na.rm = TRUE), .groups = "drop")
  
  p <- dados_gene_plot %>% 
    ggplot() +
    geom_brain(
      atlas = ggseg::dk(),
      position = position_brain(c("right lateral", "right medial")),
      mapping = aes(fill = log2(Expressao_Media))
    ) +
    scale_fill_gradient2(
      low = "royalblue",
      mid = "firebrick",
      high = "goldenrod",
      midpoint = 0,
      name = "Score"
    ) +
    facet_wrap(~ factor(broad_age, levels = unique(dados_gene_plot$broad_age)), ncol = 2) +
    theme_void() +
    theme(
      strip.text = element_text(size = 14, face = "bold", margin = margin(b = 5, t = 2)),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  print(p)
}

#* Plota o mapa cerebral baseado no Score ssGSEA de uma ONTOLOGIA
#* @param geneset A ontologia selecionada
#* @serializer png list(width = 900, height = 700, res = 120)
#* @get /plot_ontology
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS") {
  
  # Lógica Otimizada: Usa a tabela_base_ontologias já cruzada na memória!
  dados_prontos <- tabela_base_ontologias %>%
    filter(GeneSet == geneset) %>%
    group_by(broad_age, region) %>%
    summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  
  if(nrow(dados_prontos) == 0) {
    stop("Ontologia não encontrada ou sem dados para plotagem.")
  }
  
  p <- dados_prontos %>% 
    ggplot() +
    geom_brain(
      atlas = ggseg::dk(),
      position = position_brain(c("right lateral", "right medial")),
      mapping = aes(fill = Score_Medio_ssGSEA)
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      name = "ssGSEA Score"
    ) +
    facet_wrap(~ factor(broad_age, levels = unique(dados_prontos$broad_age)), ncol = 2) +
    theme_void() +
    theme(
      strip.text = element_text(size = 14, face = "bold", margin = margin(b = 5, t = 2)),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  print(p)
}
