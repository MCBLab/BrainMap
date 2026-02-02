# Use a imagem base oficial do R com Shiny
FROM rocker/shiny-verse:4.3.2

# Instale dependências do sistema necessárias para ggseg e outros pacotes
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libgsl-dev \
    libgslcblas0 \
    libudunits2-dev \
    libproj-dev \
    libgdal-dev \
    libgeos-dev \
    libsqlite3-dev \
    libspatialindex-dev \
    python3-dev \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Instale pacotes Python necessários para ggseg (caso use reticulate)
RUN pip3 install numpy pandas

# Instale pacotes R específicos necessários para o seu script
RUN R -e "install.packages(c('shiny', 'dplyr', 'tidyr', 'ggplot2', 'GSVA', 'shinythemes', 'svglite', 'shinycssloaders', 'remotes', 'devtools'), dependencies = TRUE, repos='https://cloud.r-project.org/')"

# Instale dependências do Bioconductor primeiro
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c('GSVA', 'Biobase', 'GSEABase'), update = FALSE, ask = FALSE)"

# Instale dependências específicas do ggseg antes
RUN R -e "install.packages(c('sf', 'rgeos', 'rgdal', 'sp', 'maps', 'mapdata', 'maptools', 'rnaturalearth', 'rnaturalearthdata', 'ggrepel', 'viridis'), dependencies = TRUE, repos='https://cloud.r-project.org/')"

# Instale ggseg e ggsegExtra do GitHub
RUN R -e "remotes::install_github('LCBC-UiO/ggseg@main', dependencies = TRUE)"
RUN R -e "remotes::install_github('LCBC-UiO/ggsegExtra@main', dependencies = TRUE)"

# Instale pacotes adicionais que ggseg pode precisar
RUN R -e "install.packages(c('brainconn', 'cowplot', 'ggridges', 'patchwork', 'scales', 'RColorBrewer', 'viridisLite'), dependencies = TRUE, repos='https://cloud.r-project.org/')"

# Verifique a instalação do ggseg
RUN R -e "if (!require('ggseg')) { cat('ggseg not installed!\\n'); } else { cat('ggseg installed successfully!\\n'); library(ggseg); print(packageVersion('ggseg')) }"

# Crie diretório para a aplicação
RUN mkdir -p /srv/shiny-server/app

# Copie o arquivo do aplicativo
COPY app.R /srv/shiny-server/app/
# Copie o arquivo de dados
# COPY dados_otimizados.rds /srv/shiny-server/app/

# Configure permissões
RUN chown -R shiny:shiny /srv/shiny-server/app

# Exponha a porta
EXPOSE 3838

# Comando para rodar o aplicativo
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app', host = '0.0.0.0', port = 3838)"]
