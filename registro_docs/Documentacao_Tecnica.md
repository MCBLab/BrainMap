# Documentação Técnica - VEGAbrain

## Especificações do Sistema

O sistema *VEGAbrain* é estruturado em uma arquitetura cliente-servidor (Frontend e API Backend) rodando em contêineres Docker independentes para escalabilidade e isolamento.

- **Backend (API):** Desenvolvido na linguagem estatística R (v4.3.0+), rodando o pacote `plumber` para expor rotas HTTP (RESTful). A API carrega em memória um grande volume de matrizes otimizadas provenientes do BrainSpan Atlas para resposta às requisições de plotagem e download numérico.
- **Frontend:** Aplicação do tipo Single Page Application (SPA), construída com React e pacote de bundler Vite, utilizando requisições assíncronas para interface com a API back-end.
- **Mecanismos de Processamento:** A análise estatística de enriquecimento celular utiliza a biblioteca R Bioconductor `GSVA`.

## Fluxograma de Funcionamento da Aplicação

1. **Acesso do Usuário:** O usuário interage com o front-end React na porta local / web padrão.
2. **Entrada de Dados:** O usuário escolhe um gene único ou providencia uma lista customizada (assinatura gênica).
3. **Requisição HTTP:** O frontend envia a lista encapsulada em um POST/GET request para a API Plumber rodando na porta alocada (ex. 33857).
4. **Processamento Back-end:**
   - Para genes únicos: A API localiza as quantificações na matriz BrainSpan, filtra e constrói um *data frame*.
   - Para listas de genes (GSVA): O pacote R `GSVA` constrói pontuações ssGSEA de enriquecimento amostra-a-amostra.
5. **Renderização Visual:** A biblioteca geométrica gráfica `ggseg` mapeia os escores nas subdivisões anatômicas de imagens hemisféricas laterais e mediais cerebrais para as 5 faixas etárias.
6. **Resposta:** O plot codificado como SVG e os respectivos dados na forma de `.csv` são devolvidos ao Front-end, que o apresenta na tela do usuário para análise e exportação.

---

## Listagem do Código Fonte (Parcial)

Abaixo é apresentada uma listagem do principal orquestrador de requisições, o script da API em `plumber.R`:

```r
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

dados_app <- readRDS("dados_otimizados.rds")
tabela_final_ssgsea <- read.csv("ontologyssGSEA.csv")

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

# Helper function to create the lateral and medial brain plots grid
build_brain_grid <- function(data, fill_var, fill_scale, caption_text = NULL) {
  p_lateral <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::dk(),
      position = position_brain(c("right lateral")),
      mapping = aes(fill = .data[[fill_var]])
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
      mapping = aes(fill = .data[[fill_var]])
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    labs(caption = caption_text) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(t = 0, r = 5, b = 10, l = 5),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.caption = element_text(color = "firebrick", face = "bold", size = 11, hjust = 0.5, margin = margin(t = 15))
    )

  plot_grid(p_lateral, p_medial, nrow = 2, align = "v", rel_heights = c(1, 1))
}

#* Retorna o plot do cérebro com a expressão do gene selecionado
#* @param gene Nome do gene selecionado no front-end
#* @serializer svg list(width = 8, height = 4)
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

  escala_cores <- scale_fill_gradient(low = "blue", high = "orange", name = "Log2 Expr", na.value = "darkgray")
  
  plots_brain <- build_brain_grid(dados_gene_plot, "Expressao_Media", escala_cores)
  print(plots_brain)
}
```
