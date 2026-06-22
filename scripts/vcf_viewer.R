#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(htmlwidgets)
  library(htmltools)
})

# =========================================================
# DEBUG MODE
# =========================================================
# args <- c(
#   "raw/R_PTA_22_2.vcf",
#   "raw/R_PTA_22.vcf",
#   "raw/symbol_list_wes_3.txt",
#   "vcf_output.html"
# )

# =========================================================
# CLI MODE (uncomment in production)
# =========================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage:\n",
    "Rscript vcf_viewer.R <vcf1> <vcf2> [gene_list] <out_html>\n"
  )
}

vcf_path1 <- args[1]
vcf_path2 <- args[2]

gene_list_path <- if (length(args) == 4) args[3] else NA
out_html <- args[length(args)]

# =========================================================
# gene list
# =========================================================
gene_list <- NULL

if (!is.na(gene_list_path)) {
  gene_list <- readLines(gene_list_path)
  gene_list <- gene_list[gene_list != ""]
}

gene_list_str <- if (!is.null(gene_list)) {
  paste(gene_list, collapse = ", ")
} else {
  "No gene filter provided"
}

# =========================================================
# helper: bcftools command
# =========================================================
get_bcftools_cmd <- function(path) {
  vcf <- read.vcfR(path, verbose = FALSE)
  meta <- vcf@meta
  
  cmd <- grep("##bcftools_viewCommand", meta, value = TRUE)
  
  if (length(cmd) == 0) return("No bcftools command found")
  
  sub("^##bcftools_viewCommand=", "", cmd)
}

# =========================================================
# VCF parser
# =========================================================
make_vcf_table <- function(path, gene_filter = NULL) {
  
  vcf <- read.vcfR(path, verbose = FALSE)
  df <- as.data.frame(vcf@fix)
  
  # INFO
  df$DP_INFO <- suppressWarnings(as.numeric(vcfR::extract.info(vcf, "DP")))
  
  # FORMAT
  GT  <- vcfR::extract.gt(vcf, "GT")[, 1]
  GQ  <- vcfR::extract.gt(vcf, "GQ")[, 1]
  DP  <- vcfR::extract.gt(vcf, "DP")[, 1]
  AD  <- vcfR::extract.gt(vcf, "AD")[, 1]
  VAF <- vcfR::extract.gt(vcf, "VAF")[, 1]
  PL  <- vcfR::extract.gt(vcf, "PL")[, 1]
  
  df$GT <- GT
  df$GQ <- suppressWarnings(as.numeric(GQ))
  df$DP <- suppressWarnings(as.numeric(DP))
  df$VAF <- suppressWarnings(as.numeric(VAF))
  df$PL <- PL
  
  # AD split
  ad_split <- strsplit(AD, ",")
  
  df$AD_REF <- suppressWarnings(as.numeric(sapply(ad_split, `[`, 1)))
  df$AD_ALT <- suppressWarnings(as.numeric(sapply(ad_split, function(x)
    if(length(x) >= 2) x[2] else NA
  )))
  
  # HGVS
  df$HGVS <- paste0(df$CHROM, ":g.", df$POS, df$REF, ">", df$ALT)
  
  # CSQ
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
      separate(CSQ_split, into = csq_cols, sep = "\\|", fill = "right", extra = "drop")
  }
  
  # gene grouping
  if (!is.null(gene_filter) && "SYMBOL" %in% colnames(df)) {
    df$group <- ifelse(df$SYMBOL %in% gene_filter,
                       "selected genes",
                       "other genes")
  } else {
    df$group <- "all"
  }
  
  keep_cols <- c(
    "HGVS","DP","AD_REF","AD_ALT","VAF","QUAL",
    "Consequence","IMPACT","SYMBOL","Gene",
    "gnomADe_AF","CLIN_SIG","PUBMED","group"
  )
  
  df_small <- df %>% select(any_of(keep_cols))
  
  datatable(
    df_small,
    filter = "top",
    rownames = FALSE,
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      autoWidth = TRUE
    )
  )
}

# =========================================================
# load commands (IMPORTANT FIX)
# =========================================================
cmd1 <- get_bcftools_cmd(vcf_path1)
cmd2 <- get_bcftools_cmd(vcf_path2)

# =========================================================
# tables
# =========================================================
tbl1 <- make_vcf_table(vcf_path1, gene_list)
tbl2 <- make_vcf_table(vcf_path2, NULL)

# =========================================================
# UI
# =========================================================
ui <- tagList(
  
  h2("VCF Viewer"),
  
  tags$div(
    style = "padding:10px; background:#f5f5f5;",
    tags$b("Gene list (selected genes): "),
    gene_list_str
  ),
  
  tags$ul(
    class = "nav nav-tabs",
    tags$li(class = "active",
            tags$a(href = "#tab1", `data-toggle` = "tab", "selected genes")),
    tags$li(
      tags$a(href = "#tab2", `data-toggle` = "tab", "other genes"))
  ),
  
  tags$div(class = "tab-content",
           
           tags$div(class = "tab-pane active", id = "tab1",
                    tags$h4("bcftools command (VCF1)"),
                    tags$pre(cmd1),
                    tbl1
           ),
           
           tags$div(class = "tab-pane", id = "tab2",
                    tags$h4("bcftools command (VCF2)"),
                    tags$pre(cmd2),
                    tbl2
           )
  )
)

# =========================================================
# save
# =========================================================
htmltools::save_html(
  ui,
  file = out_html,
  selfcontained = TRUE
)

cat("Saved:", out_html, "\n")