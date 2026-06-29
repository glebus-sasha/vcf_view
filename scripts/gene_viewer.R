#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(htmltools)
})

# =========================================================
# DEBUG
# =========================================================
args <- c(
  "raw/R_PTA_22_2.vcf",
  "raw/R_PTA_22.vcf",
  "raw/carcinoma_vs_normal_gene_names_added.tsv",
  "raw/symbol_list_wes_3.txt",
  "gene_view.html"
)

# =========================================================
# CLI
# =========================================================
# args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: <special.vcf> <other.vcf> <deg.tsv> <gene_list.txt> <out.html>")
}

vcf_special <- args[1]
vcf_other   <- args[2]
deg_path    <- args[3]
gene_list_p <- args[4]
out_html    <- args[5]

# =========================================================
# INPUTS
# =========================================================
gene_list <- readLines(gene_list_p)
gene_list <- gene_list[gene_list != ""]

deg <- read.delim(deg_path, stringsAsFactors = FALSE)
deg$gene_id_clean <- sub("\\..*$", "", deg$gene_id)

# =========================================================
# VCF PARSER
# =========================================================
parse_vcf <- function(path, label) {
  
  vcf <- read.vcfR(path, verbose = FALSE)
  df <- as.data.frame(vcf@fix)
  
  csq_header <- grep("##INFO=<ID=CSQ", vcf@meta, value = TRUE)
  
  if (length(csq_header) > 0) {
    
    csq_header <- strsplit(csq_header, "Format: ")[[1]][2]
    csq_cols <- strsplit(csq_header, "\\|")[[1]]
    
    df <- df %>%
      mutate(
        CSQ_split = strsplit(INFO, ";") %>%
          lapply(function(x) {
            csq <- grep("^CSQ=", x, value = TRUE)
            if (length(csq) == 0) return(NA)
            sub("^CSQ=", "", csq)
          })
      ) %>%
      unnest(CSQ_split) %>%
      separate(CSQ_split,
               into = csq_cols,
               sep = "\\|",
               fill = "right",
               extra = "drop")
  }
  
  df$gene_id_clean <- sub("\\..*$", "", df$Gene)
  df$wes_source <- label
  
  df
}

# =========================================================
# LOAD
# =========================================================
s <- parse_vcf(vcf_special, "special")
o <- parse_vcf(vcf_other, "other")

vcf_all <- bind_rows(s, o)

# =========================================================
# IMPACT AGGREGATION PER GENE
# =========================================================
impact <- vcf_all %>%
  group_by(SYMBOL, gene_id_clean, wes_source) %>%
  summarise(
    H = sum(IMPACT == "HIGH", na.rm = TRUE),
    M = sum(IMPACT == "MODERATE", na.rm = TRUE),
    L = sum(IMPACT == "LOW", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    WES_IMPACT = paste0("H=", H, " | M=", M, " | L=", L),
    
    # scoring for coloring
    IMPACT_SCORE = H * 3 + M * 2 + L * 1
  )

# =========================================================
# DEG
# =========================================================
deg_small <- deg %>%
  select(gene_id_clean, log2FoldChange, padj)

# =========================================================
# FINAL MERGE (ONE TABLE)
# =========================================================
final <- impact %>%
  left_join(deg_small, by = "gene_id_clean") %>%
  mutate(
    group = ifelse(SYMBOL %in% gene_list, "selected", "other")
  ) %>%
  select(
    SYMBOL,
    log2FoldChange,
    padj,
    WES_IMPACT,
    IMPACT_SCORE,
    wes_source,
    group
  )

# =========================================================
# COLORING HELPERS
# =========================================================
log2_col <- function(x) {
  ifelse(is.na(x), "",
         ifelse(x > 1, "background-color:#d7191c;color:white;",
                ifelse(x < -1, "background-color:#2c7bb6;color:white;",
                       "")))
}

padj_col <- function(x) {
  ifelse(is.na(x), "",
         ifelse(x < 0.05, "background-color:#fdae61;",
                ""))
}

impact_col <- function(x) {
  ifelse(is.na(x), "",
         ifelse(x >= 5, "background-color:#7f0000;color:white;",
                ifelse(x >= 3, "background-color:#d73027;color:white;",
                       ifelse(x >= 1, "background-color:#fee08b;",
                              ""))))
}

# =========================================================
# TABLES
# =========================================================
make_table <- function(df) {
  
  datatable(df, filter = "top", options = list(pageLength = 25, scrollX = TRUE)) %>%
    formatStyle("log2FoldChange", backgroundColor = styleInterval(c(-1, 1),
                                                                  c("#2c7bb6", "white", "#d7191c"))) %>%
    formatStyle("padj", backgroundColor = styleInterval(0.05,
                                                        c("#fdae61", "white"))) %>%
    formatStyle("IMPACT_SCORE", backgroundColor = styleInterval(c(1, 3, 5),
                                                                c("#ffffcc", "#fed976", "#d95f0e", "#7f0000")),
                color = styleInterval(c(1, 3, 5),
                                      c("black", "black", "black", "white")))
}

tbl_selected <- make_table(final %>% filter(group == "selected"))
tbl_other    <- make_table(final %>% filter(group == "other"))

# =========================================================
# UI
# =========================================================
ui <- tagList(
  
  h2("WES Gene Integrated Viewer"),
  
  tags$ul(
    class = "nav nav-tabs",
    
    tags$li(class = "active",
            tags$a(href = "#tab1", `data-toggle` = "tab", "Selected genes")),
    
    tags$li(tags$a(href = "#tab2", `data-toggle` = "tab", "Other genes"))
  ),
  
  tags$div(class = "tab-content",
           
           tags$div(
             class = "tab-pane active",
             id = "tab1",
             tags$h3("Selected genes"),
             tbl_selected
           ),
           
           tags$div(
             class = "tab-pane",
             id = "tab2",
             tags$h3("Other genes"),
             tbl_other
           )
  )
)

# =========================================================
# SAVE
# =========================================================
htmltools::save_html(ui, file = out_html)

cat("Saved:", out_html, "\n")