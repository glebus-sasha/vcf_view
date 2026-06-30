#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(htmltools)
})

# =========================================================
# DEBUG (with thresholds)
# =========================================================
# args <- c(
#   "raw/R_PTA_22_2.vcf",
#   "raw/R_PTA_22.vcf",
#   "raw/carcinoma_vs_normal_gene_names_added.tsv",
#   "raw/symbol_list_wes_3.txt",
#   "gene_view.html",
#   "1.0",      # log2FC threshold
#   "0.05"      # padj threshold
# )

# =========================================================
# CLI
# =========================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop("Usage: <special.vcf> <other.vcf> <deg.tsv> <gene_list.txt> <out.html> <log2fc_threshold> <padj_threshold>")
}

vcf_special <- args[1]
vcf_other   <- args[2]
deg_path    <- args[3]
gene_list_p <- args[4]
out_html    <- args[5]
log2fc_cut  <- as.numeric(args[6])
padj_cut    <- as.numeric(args[7])

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
  
  if (!file.exists(path) || file.size(path) == 0) {
    warning("File ", path, " is empty or missing.")
    return(data.frame())
  }
  
  vcf <- tryCatch(
    read.vcfR(path, verbose = FALSE),
    error = function(e) {
      warning("Error reading ", path, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(vcf)) return(data.frame())
  
  df <- as.data.frame(vcf@fix)  # contains CHROM, POS, REF, ALT, etc.
  
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
  
  if (!"Gene" %in% colnames(df)) {
    warning("No Gene column in VCF; cannot extract gene IDs.")
    return(data.frame())
  }
  
  df$gene_id_clean <- sub("\\..*$", "", df$Gene)
  df$wes_source <- label
  
  # Ensure Consequence column exists
  if (!"Consequence" %in% colnames(df)) {
    df$Consequence <- NA_character_
  }
  
  # Keep CHROM, POS, REF, ALT already present
  df
}

# =========================================================
# LOAD
# =========================================================
s <- parse_vcf(vcf_special, "special")
o <- parse_vcf(vcf_other, "other")

vcf_all <- bind_rows(s, o)

# =========================================================
# AGGREGATION: IMPACT + VARIANTS
# =========================================================
if (nrow(vcf_all) > 0 && all(c("SYMBOL", "IMPACT", "gene_id_clean", "wes_source", 
                               "Consequence", "CHROM", "POS", "REF", "ALT") %in% colnames(vcf_all))) {
  
  # First, create a variant identifier
  vcf_all <- vcf_all %>%
    mutate(VariantID = paste(CHROM, POS, REF, ALT, sep = ":"))
  
  # Group by gene, source, and variant to collect consequences
  variant_consequences <- vcf_all %>%
    group_by(SYMBOL, gene_id_clean, wes_source, VariantID, CHROM, POS, REF, ALT) %>%
    summarise(
      Conseq = paste(unique(Consequence[!is.na(Consequence)]), collapse = ", "),
      H = sum(IMPACT == "HIGH", na.rm = TRUE),
      M = sum(IMPACT == "MODERATE", na.rm = TRUE),
      L = sum(IMPACT == "LOW", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Build variant string: "chr:g.POSREF>ALT | Consequence"
      VariantStr = paste0(
        CHROM, ":g.", POS, REF, ">", ALT,
        ifelse(Conseq != "", paste0(" | ", Conseq), "")
      )
    )
  
  # Now aggregate by gene and source, summing H/M/L and collecting variant strings
  impact <- variant_consequences %>%
    group_by(SYMBOL, gene_id_clean, wes_source) %>%
    summarise(
      H = sum(H, na.rm = TRUE),
      M = sum(M, na.rm = TRUE),
      L = sum(L, na.rm = TRUE),
      Variants = paste(VariantStr, collapse = "<br>"),
      .groups = "drop"
    ) %>%
    mutate(
      # Build WES_IMPACT with only non-zero categories
      WES_IMPACT = paste(
        ifelse(H > 0, paste0("H=", H), ""),
        ifelse(M > 0, paste0("M=", M), ""),
        ifelse(L > 0, paste0("L=", L), ""),
        sep = " | "
      ),
      # Clean up extra separators
      WES_IMPACT = gsub("^ \\| | \\| $| \\| \\| ", "", WES_IMPACT),
      WES_IMPACT = ifelse(WES_IMPACT == "", "", WES_IMPACT)
    )
} else {
  impact <- data.frame()
}

# =========================================================
# DEG
# =========================================================
deg_small <- deg %>%
  select(gene_id_clean, log2FoldChange, padj) %>%
  mutate(
    log2FoldChange = round(log2FoldChange, 3),
    padj = round(padj, 3)
  )

# =========================================================
# FINAL MERGE (ONE TABLE)
# =========================================================
if (nrow(impact) > 0) {
  final <- impact %>%
    left_join(deg_small, by = "gene_id_clean") %>%
    select(
      SYMBOL,
      log2FoldChange,
      padj,
      wes_source,
      WES_IMPACT,
      Variants
    )
} else {
  final <- data.frame()
}

# =========================================================
# COLORING HELPERS (soft pastel)
# =========================================================
# Impact coloring: prioritize HIGH > MODERATE > LOW
impact_col <- function(x) {
  if (is.na(x) || x == "") return("")
  # parse counts: "H=2 | M=1" or "H=2" etc.
  h <- as.numeric(gsub(".*H=([0-9]+).*", "\\1", x))
  m <- as.numeric(gsub(".*M=([0-9]+).*", "\\1", x))
  l <- as.numeric(gsub(".*L=([0-9]+).*", "\\1", x))
  if (is.na(h)) h <- 0
  if (is.na(m)) m <- 0
  if (is.na(l)) l <- 0
  
  if (h > 0) return("background-color:#fccde5;")   # soft pink
  if (m > 0) return("background-color:#d9d9d9;")   # soft grey
  if (l > 0) return("background-color:#e0f3f8;")   # soft cyan
  return("")
}

# =========================================================
# BUILD TABLE (with protection and styling)
# =========================================================
if (nrow(final) == 0) {
  # Empty table – show message
  ui <- tagList(
    tags$head(
      tags$style(HTML("
        body { font-family: 'Segoe UI', Arial, sans-serif; background-color: #f5f7fa; padding: 20px; }
        .container { max-width: 1200px; margin: 0 auto; background-color: white; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); padding: 25px; }
        h2 { color: #2c3e50; border-bottom: 2px solid #ecf0f1; padding-bottom: 10px; }
      "))
    ),
    tags$div(class = "container",
             h2("WES Gene Integrated Viewer"),
             div(
               style = "margin-top: 2em; padding: 1em; background-color: #f9f9f9; border: 1px solid #ddd; border-radius: 4px;",
               p("No data available. Please check your input files.")
             )
    )
  )
} else {
  
  # Prepare gene list info
  gene_list_str <- paste(gene_list, collapse = ", ")
  
  dat <- datatable(
    final,
    filter = "top",
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      columnDefs = list(
        list(targets = c("log2FoldChange", "padj"), class = "dt-right"),
        list(targets = which(colnames(final) == "Variants"), render = JS("function(data, type, row) { return data; }"))
      )
    ),
    rownames = FALSE,
    class = "display compact stripe hover row-border",
    escape = FALSE   # allow HTML in Variants column
  ) %>%
    formatStyle(
      "log2FoldChange",
      backgroundColor = styleInterval(
        c(-log2fc_cut, log2fc_cut),
        c("#b3cde3", "white", "#fbb4ae")
      )
    ) %>%
    formatStyle(
      "padj",
      backgroundColor = styleInterval(
        padj_cut,
        c("#fed9a6", "white")
      )
    ) %>%
    formatStyle(
      "WES_IMPACT",
      backgroundColor = styleEqual(
        unique(final$WES_IMPACT),
        sapply(unique(final$WES_IMPACT), impact_col)
      )
    )
  
  ui <- tagList(
    tags$head(
      tags$style(HTML("
        body { font-family: 'Segoe UI', Arial, sans-serif; background-color: #f5f7fa; padding: 20px; }
        .container { max-width: 1200px; margin: 0 auto; background-color: white; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); padding: 25px; }
        h2 { color: #2c3e50; border-bottom: 2px solid #ecf0f1; padding-bottom: 10px; }
        .gene-list-box { background-color: #f8f9fa; border-left: 4px solid #3498db; padding: 12px 18px; margin-bottom: 20px; border-radius: 4px; }
        .gene-list-box strong { color: #2c3e50; }
        .dataTables_wrapper .dataTables_filter input { border-radius: 4px; border: 1px solid #ccc; padding: 4px 8px; }
        .dataTables_wrapper .dataTables_length select { border-radius: 4px; border: 1px solid #ccc; padding: 4px; }
        table.dataTable thead th { background-color: #f1f3f5; color: #2c3e50; font-weight: 600; }
        table.dataTable tbody tr:hover { background-color: #f0f7ff; }
        /* Prevent line breaks in other columns */
        .dataTable td { white-space: nowrap; }
        .dataTable td:last-child { white-space: normal; }  /* allow wrap in Variants */
      "))
    ),
    tags$div(class = "container",
             h2("WES Gene Integrated Viewer"),
             
             tags$div(class = "gene-list-box",
                      tags$strong("Genes of interest:"),
                      tags$span(gene_list_str, style = "margin-left: 8px;")
             ),
             
             dat
    )
  )
}

# =========================================================
# SAVE
# =========================================================
htmltools::save_html(ui, file = out_html)

cat("Saved:", out_html, "\n")