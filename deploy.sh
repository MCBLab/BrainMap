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
# --min-instances=1 --no-cpu-throttling: mantem uma instancia quente.
# --memory=4Gi: o boot mede ~450 MB de dados; o resto e folga para os pacotes
#   e para o gsva() sob demanda. Confira a metrica de memoria no console e
#   ajuste depois do primeiro deploy.
gcloud run deploy brainmap-api \
  --image="$IMG" \
  --region="$REGION" \
  --allow-unauthenticated \
  --memory=4Gi --cpu=2 \
  --concurrency=4 \
  --timeout=900 \
  --min-instances=1 --max-instances=5 \
  --no-cpu-throttling
