library(vroom)
library(dplyr)
library(tibble)

setwd("~/BrainMap/") 

rows <- vroom("genes_matrix_csv/rows_metadata.csv")
columns <- vroom("genes_matrix_csv/columns_metadata.csv")

counts_raw <- vroom(
  "genes_matrix_csv/expression_matrix.csv",
  col_names = FALSE
) %>%
  dplyr::select(-1)

#Matrix numerica
matriz_expressao <- as.matrix(counts_raw)

colnames(matriz_expressao) <- columns$column_num
rownames(matriz_expressao) <- rows$gene_symbol

# 3. Pré-calcular Log2 (Para não fazer isso a cada clique)
# Lidando com duplicatas de genes (somando ou tirando média) se houver,
# ou apenas removendo duplicatas para garantir rownames únicos
genes_unicos <- !duplicated(rownames(matriz_expressao))
matriz_expressao <- matriz_expressao[genes_unicos, ]
matriz_final <- log2(matriz_expressao + 1)

dados_app <- list(
  expression_matrix = matriz_final,
  col_meta = columns,
  gene_list = rownames(matriz_final)
)

saveRDS(dados_app, "dados_otimizados.rds")
