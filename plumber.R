library(plumber)
library(dplyr)
library(ggplot2)
library(ggseg)
library(GSVA)
library(cowplot)
library(paletteer)
library(PNWColors)

#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  plumber::forward()
}

dados_app <- readRDS("dados_otimizados.rds")
tabela_final_ssgsea <- read.csv("ontologyssGSEA.csv")
genesets_list <- readRDS("genesets_list.rds")

matrix_dados <- dados_app$expression_matrix

meta_limpo <- dados_app$col_meta %>%
  rename(region = structure_mapped) %>%
  filter(!is.na(region), !is.na(broad_age)) %>%
  select(column_num, region, broad_age, macro_region)

tabela_base_ontologias <- tabela_final_ssgsea %>%
  mutate(Amostra = as.numeric(Amostra)) %>%
  inner_join(meta_limpo, by = c("Amostra" = "column_num"), relationship = "many-to-many")

setgene <- rownames(matrix_dados)
setontologies <- unique(tabela_final_ssgsea$GeneSet)

#* @get /list_genes
function() {
  setgene
}

#* @get /list_ontologies
function() {
  setontologies
}

# Helper function to create the lateral and medial brain plots grid
build_brain_grid <- function(data, fill_var, fill_scale, caption_text = NULL) {
  p_lateral <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::dk(),
      position = position_brain(c("right lateral")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    theme_void() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  p_medial <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::dk(),
      position = position_brain(c("right medial")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    labs(caption = caption_text) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

    p_sagittal <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::aseg(),
      position = position_brain(c("sagittal")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    labs(caption = caption_text) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

    p_coronal <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::aseg(),
      position = position_brain(c("coronal_1")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    labs(caption = caption_text) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.caption = element_text(color = "firebrick", face = "bold", size = 11, hjust = 0.5, margin = margin(t = 15))
    )

  plot_grid(p_lateral, p_medial, p_sagittal, p_coronal, nrow = 4, align = "v", rel_heights = c(1, 0.90, 1.10, 1.40))
}

#* Retorna o plot do cérebro com a expressão do gene selecionado
#* @param gene Nome do gene selecionado no front-end
#* @serializer svg list(width = 11, height = 7)
#* @get /plot_brain
function(gene = "SOX10", escala = "micro") {
  if (!(gene %in% rownames(matrix_dados))) stop("Erro: Gene não encontrado.")

  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  matrix_dados <- as.matrix(dados_app$expression_matrix)

  dados_base <- data.frame(
    Amostra = as.numeric(colnames(matrix_dados)),
    Expressao = as.numeric(matrix_dados[gene, ])
  ) %>%
    left_join(dados_app$col_meta, by = c("Amostra" = "column_num")) %>%
    rename(region = structure_mapped) %>%
    filter(!is.na(region), !is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))

  if(escala == "macro"){
    dados_gene_plot <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(Expressao_Media = mean(Expressao, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(broad_age, region, Expressao_Media)
  } else {
    dados_gene_plot <- dados_base %>%
    group_by(broad_age, region) %>%
    summarise(Expressao_Media = mean(Expressao, na.rm = TRUE), .groups = "drop")
  }

  escala_cores <- scale_fill_paletteer_c("scico::lajjola", name = "Log2 Expr", na.value = "darkgray")


  plots_brain <- build_brain_grid(dados_gene_plot, "Expressao_Media", escala_cores)
  print(plots_brain)
}

#* Plota o mapa cerebral baseado no Score ssGSEA de uma ONTOLOGIA
#* @param geneset A ontologia selecionada
#* @serializer svg list(width = 11, height = 7)
#* @get /plot_ontology
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", escala = "micro") {
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")

  dados_base <- tabela_base_ontologias %>%
    filter(GeneSet == geneset) %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))

  if (nrow(dados_base) == 0) stop("Ontologia não encontrada ou sem dados para plotagem.")

  if(escala == "macro") {
    dados_prontos <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(broad_age, region, Score_Medio_ssGSEA)
  } else {
    dados_prontos <- dados_base %>%
    group_by(broad_age, region) %>%
    summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  }

  escala <- scale_fill_viridis_c(option = "magma", name = "ssGSEA Score", na.value = "darkgray")
  
  p <- build_brain_grid(dados_prontos, "Score_Medio_ssGSEA", escala)
  print(p)
}

#* Plota o mapa cerebral baseado no ssGSEA de uma lista CUSTOMIZADA de genes
#* @param gene_string Uma string de genes
#* @serializer svg list(width = 11, height = 7)
#* @post /plot_genelist
function(gene_string = "", escala = "micro") {
  genes_brutos <- unlist(strsplit(gene_string, split = "[[:space:],]+"))
  genes_limpos <- trimws(genes_brutos)
  genes_limpos <- genes_limpos[genes_limpos != ""]

  genes_validos <- intersect(genes_limpos, rownames(matrix_dados))
  genes_invalidos <- setdiff(genes_limpos, rownames(matrix_dados))

  if (length(genes_validos) == 0) stop("Erro: Nenhum dos genes fornecidos foi encontrado na base de dados.")


  ssgsea_parametros <- GSVA::ssgseaParam(
    exprData = as.matrix(matrix_dados),
    geneSets = list(UserSet = genes_validos)
  )
  ssgsea_resultado <- GSVA::gsva(ssgsea_parametros, verbose = FALSE)

  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")

  dados_base <- data.frame(
    column_num = as.numeric(colnames(ssgsea_resultado)),
    Score_ssGSEA = as.numeric(ssgsea_resultado["UserSet", ])
  ) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades)) 

  if(escala == "macro"){
    dados_prontos <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(broad_age, region, Score_Medio_ssGSEA)
  } else {
    dados_prontos <- dados_base %>%
    group_by(broad_age, region) %>%
    summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  }

  escala <- scale_fill_viridis_c(option = "magma", name = "Custom ssGSEA", na.value = "darkgray")
  
  p <- build_brain_grid(dados_prontos, "Score_Medio_ssGSEA", escala)
  print(p)
}

#* Baixar dados do mapa de Gene Único
#* @param gene Nome do gene
#* @serializer csv
#* @get /data_brain
function(gene = "SOX10", escala = "micro") {
  
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  
  dados_base <- data.frame(
    column_num = as.numeric(colnames(matrix_dados)), 
    Expressao = as.numeric(matrix_dados[gene, ])
  ) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))
    
  if (escala == "macro") {
    dados_finais <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      summarise(Expressao_Media = mean(Expressao, na.rm = TRUE), .groups = "drop")
  } else {
    dados_finais <- dados_base %>%
      group_by(broad_age, region) %>%
      summarise(Expressao_Media = mean(Expressao, na.rm = TRUE), .groups = "drop")
  }

  return(dados_finais)
}

#* Baixar dados da Ontologia
#* @param geneset A ontologia
#* @serializer csv
#* @get /data_ontology
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", escala = "micro") {
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
  dados_base <- tabela_base_ontologias %>%
    filter(GeneSet == geneset) %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades)) %>%
    group_by(broad_age, region) %>%
    summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")

  if (nrow(dados_base) == 0) stop("Ontologia não encontrada ou sem dados para plotagem.")

  if(escala == "macro") {
    dados_prontos <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(broad_age, region, Score_Medio_ssGSEA)
  } else {
    dados_prontos <- dados_base %>%
    group_by(broad_age, region) %>%
    summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  }

  return(dados_prontos)
}

#* Baixar dados da Lista Customizada
#* @param gene_string
#* @serializer csv
#* @post /data_genelist
function(gene_string = "", escala = "micro") {
  genes_brutos <- unlist(strsplit(gene_string, split = "[[:space:],]+"))
  genes_limpos <- trimws(genes_brutos)
  genes_validos <- intersect(genes_limpos[genes_limpos != ""], rownames(matrix_dados))
  if (length(genes_validos) == 0) {
    return(data.frame(Erro = "Nenhum gene valido"))
  }

  ssgsea_parametros <- GSVA::ssgseaParam(exprData = as.matrix(matrix_dados), geneSets = list(UserSet = genes_validos))
  ssgsea_resultado <- GSVA::gsva(ssgsea_parametros, verbose = FALSE)

  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 8)", "Adult (n = 14)")
dados_base <- data.frame(
    column_num = as.numeric(colnames(ssgsea_resultado)),
    Score_ssGSEA = as.numeric(ssgsea_resultado["UserSet", ])
  ) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades)) 

  if(escala == "macro"){
    dados_prontos <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(broad_age, region, Score_Medio_ssGSEA)
  } else {
    dados_prontos <- dados_base %>%
    group_by(broad_age, region) %>%
    summarise(Score_Medio_ssGSEA = mean(Score_ssGSEA, na.rm = TRUE), .groups = "drop")
  }
  return(dados_prontos)
}

#* Verifica quais genes da lista existem no banco de dados
#* @param gene_string
#* @post /validate_genelist
function(gene_string = "") {
  genes_brutos <- unlist(strsplit(gene_string, split = "[[:space:],]+"))
  genes_limpos <- trimws(genes_brutos)
  genes_limpos <- genes_limpos[genes_limpos != ""]

  genes_validos <- intersect(genes_limpos, rownames(matrix_dados))
  genes_invalidos <- setdiff(genes_limpos, rownames(matrix_dados))

  return(list(
    validos = genes_validos,
    invalidos = genes_invalidos
  ))
}

#* Retorna os genes da ontologia selecionada
#* @param geneset Nome da ontologia
#* @get /genes_da_via
function(geneset = "") {
  if (geneset == "" || is.null(genesets_list[[geneset]])) {
    return(list(genes = character(0)))
  }
  
  genes_da_via <- genesets_list[[geneset]]
  genes_presentes <- intersect(genes_da_via, rownames(matrix_dados))
  
  return(list(genes = as.character(genes_presentes)))
}

#* Retorna o mapa de referência de uma aba específica
#* @param view Qual corte mostrar (lateral, medial, sagittal, coronal)
#* @serializer svg list(width = 8, height = 6)
#* @get /plot_reference
function(view = "lateral") {
  
  mapa_referencia <- dados_app$mapping_info$region_to_structure %>%
    filter(!is.na(structure_name))
  
  todas_areas <- sort(unique(mapa_referencia$structure_name))
  
  mapa_referencia$structure_name <- factor(mapa_referencia$structure_name, levels = todas_areas)
  
  escala_inteligente <- scale_fill_discrete(drop = TRUE, name = "")

  if (view == "lateral") {
    df_dk <- as.data.frame(ggseg::dk())
    areas <- df_dk$region[df_dk$hemi == "right" & df_dk$view == "lateral"]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::dk(), position = position_brain("right lateral"), mapping = aes(fill = structure_name), color = "black", size = 0.3) + ggtitle("Córtex Lateral")
    
  } else if (view == "medial") {
    df_dk <- as.data.frame(ggseg::dk())
    areas <- df_dk$region[df_dk$hemi == "right" & df_dk$view == "medial"]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::dk(), position = position_brain("right medial"), mapping = aes(fill = structure_name), color = "black", size = 0.3) + ggtitle("Córtex Medial")
    
  } else if (view == "sagittal") {
    df_aseg <- as.data.frame(ggseg::aseg())
    areas <- df_aseg$region[grepl("sagittal", df_aseg$view, ignore.case = TRUE)]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::aseg(), position = position_brain("sagittal"), mapping = aes(fill = structure_name), color = "black", size = 0.3) + ggtitle("Subcortical (Sagital)")
    
  } else if (view == "coronal") {
    df_aseg <- as.data.frame(ggseg::aseg())
    areas <- df_aseg$region[grepl("coronal_1", df_aseg$view, ignore.case = TRUE)]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::aseg(), position = position_brain("coronal_1"), mapping = aes(fill = structure_name), color = "black", size = 0.3) + ggtitle("Subcortical (Coronal)")
    
  } else {
    stop("Corte não encontrado.")
  }

  p_final <- p + 
    escala_inteligente + 
    theme_void() + 
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    guides(fill = guide_legend(ncol = 2, keywidth = 1, keyheight = 1))

  print(p_final)
}