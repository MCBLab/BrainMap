library(plumber)
library(dplyr)
library(ggplot2)
library(ggseg)
library(cowplot)

#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  plumber::forward()
}

dados_app <- readRDS("dados_otimizados.rds")
# Médias já agregadas por preparaDados.R (o CSV bruto de ~9.7M linhas não é
# mais lido no boot, o que mantém o startup dentro do limite do Cloud Run)
ontologia_micro <- readRDS("ontologia_micro.rds")
ontologia_macro <- readRDS("ontologia_macro.rds")
genesets_list <- readRDS("genesets_list.rds")

matrix_dados <- dados_app$expression_matrix

meta_limpo <- dados_app$col_meta %>%
  rename(region = structure_mapped) %>%
  filter(!is.na(region), !is.na(broad_age)) %>%
  select(column_num, region, broad_age, macro_region)

# Combinações (idade, região) efetivamente amostradas: usadas para espalhar o
# score da macro-região de volta nas regiões que o ggseg desenha
regioes_por_macro <- meta_limpo %>%
  filter(!is.na(macro_region)) %>%
  distinct(broad_age, region, macro_region)

# 2. LISTAS OTIMIZADAS (À prova de falhas)
setgene <- rownames(matrix_dados)
setontologies <- sort(unique(ontologia_micro$GeneSet))

# --- Significancia por estrutura (asteriscos no mapa) ------------------------
# Mesma estatistica da publicacao (publication/figures/enrichment_stats.R): por
# (estrutura, janela), o desvio da lista em relacao a media da janela e comparado
# com o de 4.000 listas aleatorias de mesmo tamanho e mesmo nivel de expressao,
# BH sobre as celulas testadas, e a celula ainda precisa continuar significativa
# quando qualquer doador e retirado. Nao ha referencia pre-computada: tudo sai
# da matriz de expressao e dos metadados das amostras.
source("publication/figures/enrichment_stats.R")
ref_enriquecimento <- cells_from_meta(dados_app$col_meta)

# Tabelas de postos (ssGSEA por algebra, ~300 MB) e z-scores (assinatura
# ponderada, ~200 MB): construidas na primeira vez em que alguem pedir, nao no
# boot, que o Cloud Run cronometra. Com as tabelas de postos o ssGSEA de uma
# lista sai em < 1 s, no lugar dos ~45 s do GSVA::gsva().
# O nulo sorteia genes do universo do MSigDB (o mesmo tipo de gene de uma lista
# curada ou de usuario), nao das ~30 mil linhas nao anotadas da matriz.
universo_genes <- intersect(unique(unlist(genesets_list, use.names = FALSE)), setgene)
cache_estat <- new.env()
tabelas_rank <- function() {
  if (is.null(cache_estat$rt)) cache_estat$rt <- rank_tables(matrix_dados, universe = universo_genes)
  cache_estat$rt
}
referencia_ponderada <- function() {
  if (is.null(cache_estat$wref)) cache_estat$wref <- weighted_reference(matrix_dados, universe = universo_genes)
  cache_estat$wref
}

pediu_significancia <- function(x) tolower(as.character(x)) %in% c("true", "1", "yes", "sim")

# Devolve (broad_age, region, Score_Medio_ssGSEA) para uma ontologia. Na escala
# macro cada região recebe a média da macro-região a que pertence.
scores_ontologia <- function(geneset, escala) {
  if (escala == "macro") {
    ontologia_macro %>%
      filter(GeneSet == geneset) %>%
      inner_join(regioes_por_macro,
                 by = c("broad_age", "macro_region"),
                 relationship = "many-to-many") %>%
      select(broad_age, region, Score_Medio_ssGSEA)
  } else {
    ontologia_micro %>%
      filter(GeneSet == geneset) %>%
      select(broad_age, region, Score_Medio_ssGSEA)
  }
}

# --- Aba Gene List: assinatura simples ou ponderada por log2FC --------------

# Redeclarada dentro de cada endpoint; aqui em cima porque agrega_por_regiao()
# resolve o nome pelo ambiente onde foi definida, nao pelo de quem a chama.
ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)",
                  "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")

# Le a caixa de texto da aba Gene List. Cada linha e "GENE" ou "GENE <log2FC>",
# com espaco, tab, virgula ou ponto-e-virgula entre os dois. Uma linha unica
# com varios simbolos separados por virgula continua valendo (formato antigo).
# tinha_peso marca quem realmente trouxe um numero: e o que decide, mais na
# frente, entre o ssGSEA de sempre e o score ponderado.
parse_assinatura <- function(gene_string) {
  vazio <- data.frame(gene = character(0), peso = numeric(0), tinha_peso = logical(0))
  if (is.null(gene_string) || !nzchar(trimws(gene_string))) return(vazio)

  linhas <- trimws(unlist(strsplit(gene_string, "\r?\n")))
  linhas <- linhas[nzchar(linhas)]
  if (length(linhas) == 0) return(vazio)

  divide <- function(linha) {
    tokens <- unlist(strsplit(linha, "[ \t,;]+"))
    tokens[nzchar(tokens)]
  }

  # Cabecalho de export do DESeq2/edgeR/limma ("gene,log2FoldChange")
  primeiro <- divide(linhas[1])
  if (length(primeiro) >= 2 &&
      grepl("^(gene|genes|symbol|id)$", primeiro[1], ignore.case = TRUE) &&
      !is.finite(suppressWarnings(as.numeric(primeiro[2])))) {
    linhas <- linhas[-1]
  }

  partes <- lapply(linhas, function(linha) {
    tokens <- divide(linha)
    if (length(tokens) == 0) return(NULL)

    if (length(tokens) == 2) {
      valor <- suppressWarnings(as.numeric(tokens[2]))
      if (is.finite(valor)) {
        return(data.frame(gene = tokens[1], peso = valor, tinha_peso = TRUE))
      }
    }
    # Um simbolo so, ou varios sem numero nenhum: lista simples.
    data.frame(gene = tokens, peso = 1, tinha_peso = FALSE)
  })

  do.call(rbind, c(list(vazio), partes))
}

# Resolve a caixa digitada para o nome exato da matriz. Um toupper cego nao
# serve: 1.168 dos 47.808 simbolos tem minuscula (a familia C1orf112), e
# maiusculizar transformaria 1.167 genes hoje validos em "nao encontrado".
# Um unico par colide ignorando a caixa (5S_rRNA / 5S_RRNA); fica com o
# primeiro, e quem digitar a grafia exata cai nela pelo teste de igualdade.
indice_genes <- setgene[!duplicated(toupper(setgene))]
names(indice_genes) <- toupper(indice_genes)

# Parse + confronto com a matriz. Concentra tudo o que os tres endpoints da
# aba Gene List precisam saber sobre a lista submetida.
assinatura_valida <- function(gene_string) {
  bruta <- parse_assinatura(gene_string)

  # Grafia exata primeiro; so o que nao bater cai no casamento sem caixa.
  exato <- bruta$gene %in% setgene
  canonico <- as.character(ifelse(exato, bruta$gene,
                                  unname(indice_genes[toupper(bruta$gene)])))
  no_banco <- !is.na(canonico)

  validos <- bruta[no_banco, , drop = FALSE]
  validos$gene <- canonico[no_banco]
  validos <- validos[!duplicated(validos$gene), , drop = FALSE]

  ponderado <- any(validos$tinha_peso)

  list(
    validos    = validos,
    # O texto original, para a pessoa achar o simbolo na propria lista.
    invalidos  = unique(bruta$gene[!no_banco]),
    ponderado  = ponderado,
    # Sem log2FC nenhum a assinatura nao tem direcao: nao ha up nem down.
    n_up       = if (ponderado) sum(validos$peso > 0) else 0L,
    n_down     = if (ponderado) sum(validos$peso < 0) else 0L,
    n_sem_peso = if (ponderado) sum(!validos$tinha_peso) else 0L
  )
}

# score_j = sum_i(w_i * z_ij) / sum_i(|w_i|), com z = expressao padronizada
# entre as amostras. Genes sem variancia nao informam nada e saem dos dois
# lados da fracao: mante-los no denominador achataria o mapa inteiro em direcao
# ao zero.
score_ponderado <- function(genes, pesos) {
  sub <- matrix_dados[genes, , drop = FALSE]
  desv <- apply(sub, 1, sd)
  ok <- is.finite(desv) & desv > 0
  if (!any(ok)) stop("Erro: nenhum dos genes informados varia entre as amostras.")

  sub <- sub[ok, , drop = FALSE]
  pesos <- pesos[ok]
  if (sum(abs(pesos)) == 0) stop("Erro: todos os log2FC informados sao zero.")

  z <- (sub - rowMeans(sub)) / desv[ok]
  as.numeric(crossprod(z, pesos) / sum(abs(pesos)))
}

# Recebe um score por amostra e devolve a media por (janela, regiao). Na escala
# macro cada regiao desenhada recebe a media da macro-regiao a que pertence.
agrega_por_regiao <- function(score_por_amostra, escala, nome_col,
                              amostras = colnames(matrix_dados)) {
  dados_base <- data.frame(
    column_num = as.numeric(amostras),
    Score = score_por_amostra
  ) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))

  if (escala == "macro") {
    dados_prontos <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(Score_Medio = mean(Score, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(broad_age, region, Score_Medio)
  } else {
    dados_prontos <- dados_base %>%
      group_by(broad_age, region) %>%
      summarise(Score_Medio = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  names(dados_prontos)[names(dados_prontos) == "Score_Medio"] <- nome_col
  dados_prontos
}

# Estatistica por (estrutura, janela) para a assinatura submetida. Devolve
# list(celulas, motivo): `celulas` e a saida de enrich() / enrich_weighted(),
# `motivo` explica, para a legenda, por que nao houve teste.
significancia_assinatura <- function(assinatura) {
  genes <- assinatura$validos$gene
  if (length(genes) < 5) {
    return(list(celulas = NULL,
                motivo = sprintf("needs at least 5 genes found in the data (got %d)",
                                 length(genes))))
  }
  celulas <- if (assinatura$ponderado) {
    enrich_weighted(ref_enriquecimento, referencia_ponderada(), genes, assinatura$validos$peso)
  } else {
    enrich(ref_enriquecimento, tabelas_rank(), genes)
  }
  list(celulas = celulas, motivo = NULL)
}

# Parcelas que recebem asterisco: toda parcela alimentada por uma estrutura
# enriquecida (paint_outward, a mesma regra da Fig. 6 da publicacao). A
# inferencia e por estrutura do BrainSpan; a parcela so herda o veredito.
marcas_significancia <- function(estat) {
  if (is.null(estat$celulas)) return(NULL)
  paint_outward(estat$celulas, dados_app$mapping_info$region_to_structure) %>%
    filter(enriched) %>%
    distinct(broad_age, region) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))
}

legenda_significancia <- function(estat) {
  if (is.null(estat)) return(NULL)
  if (!is.null(estat$motivo)) return(paste0("Significance not tested: ", estat$motivo))
  sprintf(paste0("* q < 0.05 vs 4,000 expression-matched random gene sets (BH over %d cells),\n",
                 "robust to leaving any single donor out; %d enriched structure x window cells"),
          sum(estat$celulas$tested), sum(estat$celulas$enriched))
}

# Score ponderado quando a lista trouxe log2FC, ssGSEA de sempre quando nao.
# Com significancia = TRUE a lista devolve tambem `estat` (ver acima).
scores_assinatura <- function(assinatura, escala, significancia = FALSE) {
  genes <- assinatura$validos$gene
  if (length(genes) == 0) {
    stop("Erro: Nenhum dos genes fornecidos foi encontrado na base de dados.")
  }

  if (assinatura$ponderado) {
    coluna <- "Score_Medio_Ponderado"
    dados <- agrega_por_regiao(score_ponderado(genes, assinatura$validos$peso), escala, coluna)
  } else {
    coluna <- "Score_Medio_ssGSEA"
    # O mesmo ssGSEA do GSVA::gsva(), calculado por algebra sobre as tabelas de
    # postos (identico a 1e-6, ver validate_stats.R). GSVA normaliza dividindo
    # pela amplitude do proprio conjunto; refeito aqui, identico, para a exibicao.
    rt <- tabelas_rank()
    por_amostra <- ssgsea_from_ranks(rt, membership(match(genes, rt$genes), rt$n))[1, ]
    dados <- agrega_por_regiao(por_amostra / diff(range(por_amostra)), escala, coluna,
                               names(por_amostra))
  }

  estat <- if (significancia) significancia_assinatura(assinatura) else NULL
  list(dados = dados, coluna = coluna, estat = estat)
}
#* @get /list_genes
function() {
  setgene
}

#* @get /list_ontologies
function() {
  setontologies
}

# Helper function to create the lateral and medial brain plots grid
# Onde cai o asterisco: um ponto interior do maior poligono de cada parcela,
# calculado sobre a geometria ja reposicionada que o ggseg desenha naquela vista
# -- por isso ha um ramo por versao do ggseg. Hemisferios desenhados um sobre o
# outro (o cerebelo na vista sagital) viram um ponto so.
centroides_parcelas <- function(atlas, pos) {
  # Um ponto garantidamente DENTRO do poligono (sf::st_point_on_surface): o
  # centroide cai fora em parcelas finas e curvas (superior parietal, occipital
  # lateral na vista medial).
  medir <- function(pg, region, hemi, view) {
    bb <- sf::st_bbox(pg)
    pt <- tryCatch(sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(pg))),
                   error = function(e) sf::st_coordinates(sf::st_centroid(pg)))
    tibble(region = region, hemi = hemi, view = view, area = as.numeric(sf::st_area(pg)),
           x = pt[1, 1], y = pt[1, 2],
           xmin = bb[["xmin"]], xmax = bb[["xmax"]], ymin = bb[["ymin"]], ymax = bb[["ymax"]])
  }
  if (exists("prepare_polygon_atlas", envir = asNamespace("ggseg"), inherits = FALSE)) {
    # ggseg >= 2.2 desenha poligonos planos (x, y)
    plano <- ggseg:::prepare_polygon_atlas(atlas, position = position_brain(pos)) %>%
      filter(!is.na(region)) %>%
      mutate(chave = paste(region, hemi, .feature_id, subgroup, sep = "|"))
    aneis <- bind_rows(lapply(split(plano, plano$chave), function(d) {
      xy <- cbind(d$x, d$y)
      if (any(xy[1, ] != xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
      medir(sf::st_sfc(sf::st_polygon(list(xy))), d$region[1], d$hemi[1], d$view[1])
    }))
  } else {
    # ggseg 2.1.x (fixado no Dockerfile) desenha geometria sf reposicionada
    rep <- ggseg::reposition_brain(as.data.frame(atlas), position = pos)
    rep <- rep[!is.na(rep$region), ]
    geo <- sf::st_geometry(rep)
    aneis <- bind_rows(lapply(seq_len(nrow(rep)), function(i) {
      polys <- suppressWarnings(sf::st_cast(geo[i], "POLYGON"))
      bind_rows(lapply(seq_along(polys), function(k) medir(polys[k], rep$region[i], rep$hemi[i], rep$view[i])))
    }))
  }
  # a parcel that is only a sliver in this view (precentral seen from the medial
  # side: 0.05% of the panel) is starred where it is actually drawn, not here
  area_vista <- (max(aneis$xmax) - min(aneis$xmin)) * (max(aneis$ymax) - min(aneis$ymin))
  aneis %>%
    group_by(region, hemi, view) %>%
    slice_max(area, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(area >= 0.0025 * area_vista) %>%
    group_by(region, view) %>%
    # os dois hemisferios ficam so quando as caixas nao se sobrepoem
    filter(n() == 1 | max(xmin) > min(xmax) | max(ymin) > min(ymax) | area == max(area)) %>%
    ungroup() %>%
    select(region, hemi, view, x, y)
}


# Camada de asteriscos: `marcas` tem (broad_age, region), uma linha por parcela
# enriquecida por janela. Glifo branco sobre um halo escuro, para ler em cima de
# qualquer tom da paleta. NULL quando nao ha o que marcar (somar NULL a um
# ggplot nao faz nada).
camada_asteriscos <- function(marcas, atlas, pos, tamanho = 9) {
  if (is.null(marcas) || nrow(marcas) == 0) return(NULL)
  pontos <- marcas %>%
    inner_join(centroides_parcelas(atlas, pos), by = "region", relationship = "many-to-many")
  if (nrow(pontos) == 0) return(NULL)
  list(
    geom_text(data = pontos, aes(x = x, y = y), inherit.aes = FALSE, label = "*",
              size = tamanho * 1.2, colour = "grey20", fontface = "bold", vjust = 0.3),
    geom_text(data = pontos, aes(x = x, y = y), inherit.aes = FALSE, label = "*",
              size = tamanho, colour = "white", fontface = "bold", vjust = 0.3)
  )
}

# Os dados vao explicitamente para o geom_brain, agrupados pela faceta: o ggseg
# >= 2.2 deixou de ler os dados do ggplot(), e na 2.1.x o resultado e o mesmo.
build_brain_grid <- function(data, fill_var, fill_scale, caption_text = NULL, marcas = NULL) {
  data <- group_by(data, broad_age)
  p_lateral <- ggplot() +
    geom_brain(
      data = data,
      atlas = ggseg::dk(),
      position = position_brain(c("right lateral")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    camada_asteriscos(marcas, ggseg::dk(), "right lateral") +
    facet_wrap(~broad_age, ncol = 5) +
    theme_void() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  p_medial <- ggplot() +
    geom_brain(
      data = data,
      atlas = ggseg::dk(),
      position = position_brain(c("right medial")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    camada_asteriscos(marcas, ggseg::dk(), "right medial") +
    facet_wrap(~broad_age, ncol = 5) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

    p_sagittal <- ggplot() +
    geom_brain(
      data = data,
      atlas = ggseg::aseg(),
      position = position_brain(c("sagittal")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    camada_asteriscos(marcas, ggseg::aseg(), "sagittal") +
    facet_wrap(~broad_age, ncol = 5) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

    p_coronal <- ggplot() +
    geom_brain(
      data = data,
      atlas = ggseg::aseg(),
      position = position_brain(c("coronal_1")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    # parcelas pequenas e proximas: glifo menor para nao se sobreporem
    camada_asteriscos(marcas, ggseg::aseg(), "coronal_1", tamanho = 6) +
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

  plot_grid(p_lateral, p_medial, p_sagittal, p_coronal, nrow = 4, align = "v", rel_heights = c(1, 0.89, 1.10, 1.4))  #1, 0.90, 1.10, 1.40
}

# 3. ENDPOINTS DE PLOTAGEM (SINTAXE OFICIAL MODERNA DO GGSEG)

# Cada mapa e servido em dois formatos: SVG (vetorial, para publicacao) e PNG
# (o que a tela exibe - o SVG do ggseg carrega milhares de poligonos e pesa no
# navegador). O calculo fica nas funcoes grafico_*, e cada endpoint so escolhe
# o device do serializer, para que os dois formatos nunca divirjam.

grafico_gene <- function(gene = "SOX10", escala = "micro") {
  if (!(gene %in% rownames(matrix_dados))) stop("Erro: Gene não encontrado.")

  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")
  matrix_dados <- as.matrix(dados_app$expression_matrix)

  dados_base <- data.frame(
    Amostra = as.numeric(colnames(matrix_dados)),
    Expressao = as.numeric(matrix_dados[gene, ])
  ) %>%
    left_join(dados_app$col_meta, by = c("Amostra" = "column_num")) %>%
    rename(region = structure_mapped) %>%
    filter(!is.na(region), !is.na(broad_age)) %>%
    mutate(Expressao = ifelse(Expressao > 6, 6, Expressao)) %>%
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

  escala <- scale_fill_gradientn( colors = c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845"), 
    values = scales::rescale(c(0, 1.5, 3, 4.5, 6)),limits = c(0, 6), 
    name = "Log2 Expr", na.value = "darkgray")  

  build_brain_grid(dados_gene_plot, "Expressao_Media", escala)
}

#* Retorna o plot do cérebro com a expressão do gene selecionado (SVG)
#* @param gene Nome do gene selecionado no front-end
#* @serializer svg list(width = 11, height = 7)
#* @get /plot_brain
function(gene = "SOX10", escala = "micro") {
  print(grafico_gene(gene, escala))
}

#* Mesmo mapa do /plot_brain, em PNG
#* @param gene Nome do gene selecionado no front-end
#* @serializer png list(width = 11, height = 7, units = "in", res = 150)
#* @get /plot_brain_png
function(gene = "SOX10", escala = "micro") {
  print(grafico_gene(gene, escala))
}

grafico_ontologia <- function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", escala = "micro") {
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")

  dados_prontos <- scores_ontologia(geneset, escala)

  if (nrow(dados_prontos) == 0) stop("Ontologia não encontrada ou sem dados para plotagem.")

  dados_prontos <- dados_prontos %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))

  escala_cor <- scale_fill_gradientn( colors = c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845"),
    #values = scales::rescale(c(0, 2, 4, 6, 8)),limits = c(0, 8),
    name = "ssGSEA Score", na.value = "darkgray")

  build_brain_grid(dados_prontos, "Score_Medio_ssGSEA", escala_cor)
}

#* Plota o mapa cerebral baseado no Score ssGSEA de uma ONTOLOGIA (SVG)
#* @param geneset A ontologia selecionada
#* @serializer svg list(width = 11, height = 7)
#* @get /plot_ontology
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", escala = "micro") {
  print(grafico_ontologia(geneset, escala))
}

#* Mesmo mapa do /plot_ontology, em PNG
#* @param geneset A ontologia selecionada
#* @serializer png list(width = 11, height = 7, units = "in", res = 150)
#* @get /plot_ontology_png
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", escala = "micro") {
  print(grafico_ontologia(geneset, escala))
}

grafico_genelist <- function(gene_string = "", escala = "micro", significancia = FALSE) {
  assinatura <- assinatura_valida(gene_string)
  resultado <- scores_assinatura(assinatura, escala, significancia)

  if (assinatura$ponderado) {
    # O score soma zero entre as amostras por construcao, entao limites
    # simetricos deixam o tom claro da paleta exatamente no zero: azul e down,
    # vinho e up. Sem isso o mapa nao distingue as duas direcoes.
    limite <- max(abs(resultado$dados[[resultado$coluna]]), na.rm = TRUE)
    if (!is.finite(limite) || limite == 0) limite <- 1

    escala_cor <- list(
      scale_fill_gradientn(
        colors = c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845"),
        limits = c(-limite, limite),
        breaks = c(-limite, 0, limite),
        labels = scales::label_number(accuracy = 0.01),
        name = "Weighted Score", na.value = "darkgray"),
      # Com as cinco marcas padrao o sinal de menos faz os rotulos colidirem
      # na barra estreita; tres marcas numa barra mais larga cabem.
      guides(fill = guide_colorbar(barwidth = unit(5, "cm")))
    )

    # Um mapa ponderado e um nao ponderado sao indistinguiveis depois de
    # baixados; a legenda carimba o metodo e a composicao da assinatura.
    legenda <- sprintf("Weighted signature - %d genes (%d up / %d down)",
                       nrow(assinatura$validos), assinatura$n_up, assinatura$n_down)
    n_fora <- length(assinatura$invalidos)
    if (n_fora > 0) {
      legenda <- paste0(legenda, sprintf(" - %d symbol%s not found",
                                         n_fora, if (n_fora == 1) "" else "s"))
    }
  } else {
    escala_cor <- scale_fill_gradientn( colors = c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845"), 
      #values = scales::rescale(c(0, 2, 4, 6, 8)),limits = c(0, 8), 
      name = "Custom ssGSEA", na.value = "darkgray")
    legenda <- NULL
  }

  # A nota de significancia entra na mesma legenda carimbada no mapa: uma imagem
  # baixada tem de dizer sozinha o que o asterisco quer dizer.
  legenda <- c(legenda, legenda_significancia(resultado$estat))
  if (length(legenda) == 0) legenda <- NULL else legenda <- paste(legenda, collapse = "\n")

  build_brain_grid(resultado$dados, resultado$coluna, escala_cor, legenda,
                   marcas = marcas_significancia(resultado$estat))
}

#* Plota o mapa cerebral baseado no ssGSEA de uma lista CUSTOMIZADA de genes (SVG)
#* @param gene_string Uma string de genes
#* @param significancia "true" marca com * as estruturas enriquecidas (q < 0.05)
#* @serializer svg list(width = 11, height = 7)
#* @post /plot_genelist
function(gene_string = "", escala = "micro", significancia = "false") {
  print(grafico_genelist(gene_string, escala, pediu_significancia(significancia)))
}

#* Mesmo mapa do /plot_genelist, em PNG
#* @param gene_string Uma string de genes
#* @param significancia "true" marca com * as estruturas enriquecidas (q < 0.05)
#* @serializer png list(width = 11, height = 7, units = "in", res = 150)
#* @post /plot_genelist_png
function(gene_string = "", escala = "micro", significancia = "false") {
  print(grafico_genelist(gene_string, escala, pediu_significancia(significancia)))
}

#* Baixar dados do mapa de Gene Único
#* @param gene Nome do gene
#* @serializer csv
#* @get /data_brain
function(gene = "SOX10", escala = "micro") {
  
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")
  
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
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")

  dados_prontos <- scores_ontologia(geneset, escala)

  if (nrow(dados_prontos) == 0) stop("Ontologia não encontrada ou sem dados para plotagem.")

  dados_prontos <- dados_prontos %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades)) %>%
    arrange(broad_age, region)

  return(dados_prontos)
}

#* Baixar dados da Lista Customizada
#* @param gene_string
#* @param significancia "true" acrescenta tested / enriched / q_competitive / p_loo por regiao
#* @serializer csv
#* @post /data_genelist
function(gene_string = "", escala = "micro", significancia = "false") {
  assinatura <- assinatura_valida(gene_string)

  if (nrow(assinatura$validos) == 0) {
    return(data.frame(Erro = "Nenhum gene valido"))
  }

  resultado <- scores_assinatura(assinatura, escala, pediu_significancia(significancia))
  dados <- resultado$dados
  if (!is.null(resultado$estat$celulas)) {
    # Mesma regra do mapa: a regiao herda o veredito da(s) estrutura(s) que a
    # alimentam; o q reportado e o menor entre elas.
    por_regiao <- paint_outward(resultado$estat$celulas,
                                dados_app$mapping_info$region_to_structure) %>%
      transmute(broad_age = factor(broad_age, levels = ordem_idades), region,
                tested, enriched, q_competitive, p_loo)
    dados <- left_join(dados, por_regiao, by = c("broad_age", "region"))
  }
  dados
}

#* Verifica quais genes da lista existem no banco de dados
#* @param gene_string
#* @post /validate_genelist
function(gene_string = "") {
  assinatura <- assinatura_valida(gene_string)

  # unbox: sem isso o serializer manda [false] / [3] no lugar de escalares,
  # e o front teria de desembrulhar campo a campo.
  return(list(
    validos    = assinatura$validos$gene,
    invalidos  = assinatura$invalidos,
    ponderado  = jsonlite::unbox(assinatura$ponderado),
    n_up       = jsonlite::unbox(assinatura$n_up),
    n_down     = jsonlite::unbox(assinatura$n_down),
    n_sem_peso = jsonlite::unbox(assinatura$n_sem_peso)
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
#* @param escala Qual nível de detalhe (micro, macro)
#* @serializer svg list(width = 9, height = 7)
#* @get /plot_reference
function(view = "lateral", escala = "micro") {
  
  if (escala == "macro") {
    mapa_referencia <- dados_app$mapping_info$region_to_macro_structure
    coluna <- "macro_region"
  } else {
    mapa_referencia <- dados_app$mapping_info$region_to_structure
    coluna <- "structure_name"
  }

  mapa_referencia <- mapa_referencia %>%
    filter(!is.na(.data[[coluna]]))
  
  todas_areas <- sort(unique(mapa_referencia[[coluna]]))
  
  mapa_referencia[[coluna]] <- factor(mapa_referencia[[coluna]], levels = todas_areas)
  
  escala_sem_NA <- scale_fill_discrete(drop = TRUE, na.translate = FALSE, name = "")

  if (view == "lateral") {
    df_dk <- as.data.frame(ggseg::dk())
    areas <- df_dk$region[df_dk$hemi == "right" & df_dk$view == "lateral"]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::dk(), position = position_brain("right lateral"), mapping = aes(fill = .data[[coluna]]), color = "black", size = 0.3) + ggtitle(paste("Córtex Lateral -", ifelse(escala == "macro", "Macro", "Micro")))
    
  } else if (view == "medial") {
    df_dk <- as.data.frame(ggseg::dk())
    areas <- df_dk$region[df_dk$hemi == "right" & df_dk$view == "medial"]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::dk(), position = position_brain("right medial"), mapping = aes(fill = .data[[coluna]]), color = "black", size = 0.3) + ggtitle(paste("Córtex Medial -", ifelse(escala == "macro", "Macro", "Micro")))
    
  } else if (view == "sagittal") {
    df_aseg <- as.data.frame(ggseg::aseg())
    areas <- df_aseg$region[grepl("sagittal", df_aseg$view, ignore.case = TRUE)]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::aseg(), position = position_brain("sagittal"), mapping = aes(fill = .data[[coluna]]), color = "black", size = 0.3) + ggtitle(paste("Subcortical (Sagital) -", ifelse(escala == "macro", "Macro", "Micro")))
    
  } else if (view == "coronal") {
    df_aseg <- as.data.frame(ggseg::aseg())
    areas <- df_aseg$region[grepl("coronal_1", df_aseg$view, ignore.case = TRUE)]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::aseg(), position = position_brain("coronal_1"), mapping = aes(fill = .data[[coluna]]), color = "black", size = 0.3) + ggtitle(paste("Subcortical (Coronal) -", ifelse(escala == "macro", "Macro", "Micro")))
    
  } else {
    stop("Corte não encontrado.")
  }

  # 4. Finalização do Plot
  p_final <- p + 
    escala_sem_NA + 
    theme_void() + 
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    guides(fill = guide_legend(ncol = 2, keywidth = 1, keyheight = 1))

  print(p_final)
}
