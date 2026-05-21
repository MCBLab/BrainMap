library(plumber)
library(dplyr)
library(ggplot2)
library(ggseg)
library(GSVA)
library(cowplot)

#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  plumber::forward()
}

# 1. CARREGAMENTO DOS DADOS
dados_app <- readRDS("dados_otimizados.rds")
tabela_final_ssgsea <- read.csv("ontologyssGSEA.csv")

matrix_dados <- dados_app$expression_matrix

meta_limpo <- dados_app$col_meta %>%
  rename(region = structure_mapped) %>%
  filter(!is.na(region), !is.na(broad_age)) %>%
  select(column_num, region, broad_age)

tabela_base_ontologias <- tabela_final_ssgsea %>%
  mutate(Amostra = as.numeric(Amostra)) %>%
  inner_join(meta_limpo, by = c("Amostra" = "column_num"), relationship = "many-to-many")

# 2. LISTAS OTIMIZADAS (À prova de falhas)
setgene <- rownames(matrix_dados) 
setontologies <- unique(tabela_final_ssgsea$GeneSet)

#* @get /list_genes
function(){
  setgene
}

#* @get /list_ontologies
function() {
  setontologies
}

# 3. ENDPOINTS DE PLOTAGEM (SINTAXE OFICIAL MODERNA DO GGSEG)

#* Retorna o plot do cérebro com a expressão do gene selecionado
#* @param gene Nome do gene selecionado no front-end
#* @serializer svg list(width = 8, height = 4) 
#* @get /plot_brain
function(gene = "SOX10") { 
  
  if (!(gene %in% rownames(matrix_dados))) stop("Erro: Gene não encontrado.") 
  
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  matrix_dados <- as.matrix(dados_app$expression_matrix)
  matrix_dados <- log2(matrix_dados + 1)
  
  dados_gene_plot <- data.frame(
    Amostra = as.numeric(colnames(matrix_dados)),
    Expressao = as.numeric(matrix_dados[gene, ])) %>%
    left_join(dados_app$col_meta, by = c("Amostra" = "column_num")) %>%
    rename(region = structure_mapped) %>% 
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age)) %>%
    group_by(broad_age, region) %>%
    summarise(Expressao_Media = mean(Expressao, na.rm = TRUE), .groups = "drop")
  
  p_lateral <- ggplot() +
    geom_brain(
      data = dados_gene_plot,
      atlas = dk,
      hemi = "right",
      position = position_brain(c("right lateral")),
      mapping = aes(fill = log2(Expressao_Media))
    ) +
    scale_fill_gradient2(
      low = "royalblue", mid = "firebrick", high = "goldenrod", midpoint = 0, name = "Log2 Expr"
    ) +
    # Fatiar por ambas as variáveis. Como temos 5 idades, forçamos ncol = 5.
    # O R organizará as 10 imagens (5 idades x 2 lados) em 2 linhas de 5 colunas.
    facet_wrap(~ broad_age, ncol = 5) + 
    theme_void() +
    theme(
      # Ajustamos o texto para mostrar as labels claramente
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  p_medial <- ggplot() +
    geom_brain(
      data = dados_gene_plot,
      atlas = dk,
      hemi = "right",
      position = position_brain(c("right medial")),
      mapping = aes(fill = log2(Expressao_Media))
    ) +
    scale_fill_gradient2(
      low = "royalblue", mid = "firebrick", high = "goldenrod", midpoint = 0, name = "Log2 Expr"
    ) +
    # Fatiar por ambas as variáveis. Como temos 5 idades, forçamos ncol = 5.
    # O R organizará as 10 imagens (5 idades x 2 lados) em 2 linhas de 5 colunas.
    facet_wrap(~ broad_age, ncol = 5) + 
    theme_void() +
    theme(
      # Ajustamos o texto para mostrar as labels claramente
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(t = 0, r = 5, b = 10, l = 5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  plots_brain <- plot_grid(p_lateral, p_medial, nrow = 2, align = "v", rel_heights = c(1,1))

  print(plots_brain)
}

#* Plota o mapa cerebral baseado no Score ssGSEA de uma ONTOLOGIA
#* @param geneset A ontologia selecionada
#* @serializer svg list(width = 16, height = 6)
#* @get /plot_ontology
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS") {
  
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  
  dados_prontos <- tabela_base_ontologias %>%
    filter(GeneSet == geneset) %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades)) %>%
    group_by(broad_age, region) %>%
    summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  
  if(nrow(dados_prontos) == 0) stop("Ontologia não encontrada ou sem dados para plotagem.")
  
  p <- ggplot() +
    geom_brain(
      data = dados_prontos,
      atlas = dk,
      hemi = "right",
      mapping = aes(fill = Score_Medio_ssGSEA)
    ) +
    scale_fill_viridis_c(option = "viridis", name = "ssGSEA Score") +
    facet_wrap(~ broad_age, ncol = 5) + 
    theme_void() +
    theme(
      strip.text.x = element_text(size = 14, face = "bold", margin = margin(b = 10, t = 2)),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  print(p)
}

#* Plota o mapa cerebral baseado no ssGSEA de uma lista CUSTOMIZADA de genes
#* @param gene_string Uma string de genes 
#* @serializer svg list(width = 16, height = 6)
#* @post /plot_genelist
function(gene_string = "") {
  
  genes_brutos <- unlist(strsplit(gene_string, split = "[[:space:],]+"))
  genes_limpos <- trimws(genes_brutos)
  genes_limpos <- genes_limpos[genes_limpos != ""]
  
  genes_validos <- intersect(genes_limpos, rownames(matrix_dados))
  genes_invalidos <- setdiff(genes_limpos, rownames(matrix_dados))
  
  if (length(genes_validos) == 0) stop("Erro: Nenhum dos genes fornecidos foi encontrado na base de dados.")
  
  texto_aviso <- NULL
  if (length(genes_invalidos) > 0) {
    exemplo <- paste(head(genes_invalidos, 5), collapse = ", ")
    sufixo <- ifelse(length(genes_invalidos) > 5, "...", "")
    texto_aviso <- paste("⚠️ Aviso:", length(genes_invalidos), "gene(s) ignorado(s) (ex:", exemplo, sufixo, ")")
  }
  
  ssgsea_parametros <- GSVA::ssgseaParam(
    exprData = as.matrix(matrix_dados), 
    geneSets = list(UserSet = genes_validos)
  )
  ssgsea_resultado <- GSVA::gsva(ssgsea_parametros, verbose = FALSE)
  
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  
  dados_prontos <- data.frame(
    column_num = as.numeric(colnames(ssgsea_resultado)),
    Score_ssGSEA = as.numeric(ssgsea_resultado["UserSet", ])
  ) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades)) %>%
    group_by(broad_age, region) %>%
    summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot() +
    geom_brain(
      data = dados_prontos,
      atlas = dk,
      hemi = "right",
      mapping = aes(fill = Score_Medio_ssGSEA)
    ) +
    scale_fill_viridis_c(option = "magma", name = "Custom ssGSEA") +
    facet_wrap(~ broad_age, ncol = 5) +
    labs(caption = texto_aviso) +
    theme_void() +
    theme(
      strip.text.x = element_text(size = 14, face = "bold", margin = margin(b = 10, t = 2)),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", colour = NA),
      plot.caption = element_text(color = "firebrick", face = "bold", size = 11, hjust = 0.5, margin = margin(t = 15)) 
    )
  
  print(p)
}

#* Baixar dados do mapa de Gene Único
#* @param gene Nome do gene
#* @serializer csv
#* @get /data_brain
function(gene = "SOX10") {
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  dados <- data.frame(column_num = as.numeric(colnames(matrix_dados)), Expressao = as.numeric(matrix_dados[gene, ])) %>%
    inner_join(meta_limpo, by = "column_num") %>% filter(!is.na(broad_age)) %>% mutate(broad_age = factor(broad_age, levels = ordem_idades)) %>%
    group_by(broad_age, region) %>% summarise(Expressao_Media = mean(Expressao, na.rm = TRUE), .groups = "drop")
  return(dados)
}

#* Baixar dados da Ontologia
#* @param geneset A ontologia
#* @serializer csv
#* @get /data_ontology
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS") {
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  dados <- tabela_base_ontologias %>% filter(GeneSet == geneset) %>% filter(!is.na(broad_age)) %>% mutate(broad_age = factor(broad_age, levels = ordem_idades)) %>%
    group_by(broad_age, region) %>% summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  return(dados)
}

#* Baixar dados da Lista Customizada
#* @param gene_string
#* @serializer csv
#* @post /data_genelist
function(gene_string = "") {
  genes_brutos <- unlist(strsplit(gene_string, split = "[[:space:],]+"))
  genes_limpos <- trimws(genes_brutos)
  genes_validos <- intersect(genes_limpos[genes_limpos != ""], rownames(matrix_dados))
  if (length(genes_validos) == 0) return(data.frame(Erro="Nenhum gene valido"))
  
  ssgsea_parametros <- GSVA::ssgseaParam(exprData = as.matrix(matrix_dados), geneSets = list(UserSet = genes_validos))
  ssgsea_resultado <- GSVA::gsva(ssgsea_parametros, verbose = FALSE)
  
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  dados <- data.frame(column_num = as.numeric(colnames(ssgsea_resultado)), Score_ssGSEA = as.numeric(ssgsea_resultado["UserSet", ])) %>%
    inner_join(meta_limpo, by = "column_num") %>% filter(!is.na(broad_age)) %>% mutate(broad_age = factor(broad_age, levels = ordem_idades)) %>%
    group_by(broad_age, region) %>% summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  return(dados)
}
