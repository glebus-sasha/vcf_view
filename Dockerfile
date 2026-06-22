FROM rocker/r-ver:4.3.2

# --- System dependencies ---
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libglpk-dev \
    dos2unix \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# --- CRAN packages ---
RUN R -e "install.packages(c(
  'dplyr',
  'tidyr',
  'stringr',
  'readr',
  'ggplot2',
  'plotly',
  'htmlwidgets',
  'htmltools',
  'optparse',
  'RobustRankAggreg'
), repos='https://cloud.r-project.org')"

# --- Bioconductor packages ---
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c(
  'EnhancedVolcano',
  'WebGestaltR',
  'vcfR'
), update=FALSE, ask=FALSE)"

# --- scripts directory ---
RUN mkdir -p /usr/local/my-scripts
COPY scripts/ /usr/local/my-scripts/

# --- permissions + CRLF fix ---
RUN chmod +x /usr/local/my-scripts/*.R && dos2unix /usr/local/my-scripts/*.R

# --- symlinks (run without .R) ---
RUN for f in /usr/local/my-scripts/*.R; do \
    ln -s "$f" "/usr/local/bin/$(basename ${f%.R})"; \
done

# --- working dir ---
WORKDIR /usr/local/my-scripts

# --- default shell ---
ENTRYPOINT ["/bin/bash"]