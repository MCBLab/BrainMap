#!/usr/bin/env bash
# Build publication/main.pdf with the Springer Nature template.
#
#   bash build.sh          # build
#   bash build.sh clean    # remove LaTeX intermediates
#
# LaTeX comes from TinyTeX (installed via R: tinytex::install_tinytex()).
# tinytex::latexmk() auto-installs any LaTeX package the template pulls in,
# which is why the build goes through R rather than calling latexmk directly.
set -euo pipefail
cd "$(dirname "$0")"

if [ "${1:-}" = "clean" ]; then
  rm -f main.aux main.bbl main.blg main.log main.out main.spl main.fls main.fdb_latexmk
  echo "cleaned"
  exit 0
fi

export PATH="$HOME/.TinyTeX/bin/x86_64-linux:$PATH"

command -v pdflatex >/dev/null 2>&1 || {
  echo "error: no LaTeX found. Install it with:" >&2
  echo "  Rscript -e 'tinytex::install_tinytex()'" >&2
  exit 1
}

# clean = FALSE keeps main.log / main.blg so the build can be checked for
# undefined citations and references (see the verification notes in
# IMPLEMENTATION.md). "bash build.sh clean" removes them.
Rscript -e 'tinytex::latexmk("main.tex", engine = "pdflatex", bib_engine = "bibtex", clean = FALSE)'

echo
echo "built: $(pwd)/main.pdf"
pdfinfo main.pdf 2>/dev/null | grep -E '^(Pages|Page size)' || true
