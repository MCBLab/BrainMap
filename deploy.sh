#!/usr/bin/env bash
# Build local + push para o Artifact Registry + deploy no Cloud Run.
# O build roda localmente de proposito: preparaDados.R faz o ssGSEA completo e
# levaria horas no Cloud Build (que ainda estoura o timeout padrao de 10 min).
set -euo pipefail

PROJECT_ID="${PROJECT_ID:?defina PROJECT_ID=seu-projeto}"
REGION="${REGION:-southamerica-east1}"
REPO="${REPO:-brainmap}"
TAG="${TAG:-v1}"

IMG="${REGION}-docker.pkg.dev/${PROJECT_ID}/${REPO}/api:${TAG}"

docker build --platform linux/amd64 -f Dockerfile.api -t "$IMG" .
docker push "$IMG"

# --concurrency=4: o plumber e single-threaded, entao o default (80) so faria
#   fila dentro de um unico processo R ate estourar timeout.
# --timeout=900: /plot_genelist roda GSVA::gsva() por requisicao.
# --min-instances=0: uso esperado e de 10-200 requisicoes/mes, esporadicas.
#   Manter instancia quente custaria ~US$115-160/mes parada; escalando a zero
#   o volume cabe na free tier. O front ja tem retry ("API ligando...") que
#   absorve o cold start no carregamento das listas.
# --cpu-boost: CPU extra durante o startup, encurta o cold start.
# --memory=4Gi --cpu=2: nao vale reduzir. Nesse volume o custo e o mesmo
#   (free tier) e mais CPU deixa o cold start e o gsva() mais rapidos.
gcloud run deploy brainmap-api \
  --image="$IMG" \
  --region="$REGION" \
  --allow-unauthenticated \
  --memory=4Gi --cpu=2 \
  --concurrency=4 \
  --timeout=900 \
  --min-instances=0 --max-instances=5 \
  --cpu-throttling \
  --cpu-boost
