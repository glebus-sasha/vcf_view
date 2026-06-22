FROM rocker/r-ver:4.5.1

# System dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libuv1-dev \
    dos2unix \
    && rm -rf /var/lib/apt/lists/* 

# CRAN packages
RUN R -e "install.packages(c('dplyr','tidyr','DT','htmlwidgets','htmltools'), repos='https://cloud.r-project.org')"

# Bioconductor package
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('vcfR', update = FALSE, ask = FALSE)"

# Scripts
RUN mkdir -p /usr/local/my-scripts
COPY scripts/ /usr/local/my-scripts/

# Permissions + CRLF fix
RUN chmod +x /usr/local/my-scripts/*.R && \
    dos2unix /usr/local/my-scripts/*.R

# Symlinks
RUN for f in /usr/local/my-scripts/*.R; do \
    ln -s "$f" "/usr/local/bin/$(basename ${f%.R})"; \
done

WORKDIR /usr/local/my-scripts

ENTRYPOINT ["/bin/bash"]