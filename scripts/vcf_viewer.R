#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(htmlwidgets)
  library(htmltools)
  library(GenomicRanges)
  library(rtracklayer)
})

# =========================================================
# CLI MODE (uncomment in production)
# =========================================================
args <- commandArgs(trailingOnly = TRUE)

# =========================================================
# DEBUG MODE (comment out in production)
# =========================================================
# args <- c(
#   "raw/R_PTA_1_special.vcf",
#   "raw/R_PTA_1_other.vcf",
#   "raw/symbol_list_wes_3.txt",
#   "raw/CoreExomePanel.hg38.p12.target.v3.bed",
#   "-f 'PASS' -i '(AD[*:1] > 4) || (DP > 80)'",
#   "-f 'PASS' -i '(AD[*:1] > 4) && ((AD[*:1]/AD[*:0] > 0.05) || (DP > 80))'",
#   "(SYMBOL in symbol_list_wes_3.txt) and (IMPACT is HIGH or IMPACT is MODERATE)",
#   "(not SYMBOL in symbol_list_wes_3.txt) and (MAX_AF < 0.001 or not MAX_AF) and (IMPACT is HIGH)",
#   "vcf_output.html"
# )


if (length(args) != 9) {
  stop(
    paste(
      "Usage:",
      "Rscript vcf_viewer.R",
      "<special.vcf>",
      "<other.vcf>",
      "<gene_list.txt>",
      "<panel.bed>",
      "<bcftools_filter_special>",
      "<bcftools_filter_other>",
      "<vep_filter_special>",
      "<vep_filter_other>",
      "<output.html>"
    )
  )
}

vcf_path1 <- args[1]                  # special.vcf
vcf_path2 <- args[2]                  # other.vcf
gene_list_path <- args[3]             # symbol_list_wes_3.txt
bed_path <- args[4]                   # CoreExomePanel.hg38.p12.target.v3.bed

bcftools_filter_special <- args[5]
bcftools_filter_other   <- args[6]

vep_filter_special <- args[7]
vep_filter_other   <- args[8]

out_html <- args[9]

# =========================================================
# Helper: load gene list
# =========================================================
load_gene_list <- function(path) {
  if (!is.na(path) && file.exists(path)) {
    genes <- readLines(path)
    genes <- genes[genes != ""]
    return(genes)
  }
  return(NULL)
}

gene_list <- load_gene_list(gene_list_path)
gene_list_str <- if (!is.null(gene_list) && length(gene_list) > 0) {
  paste(gene_list, collapse = ", ")
} else {
  "No gene filter provided"
}

# =========================================================
# Helper: load BED and calculate panel size (in Mb)
# =========================================================
load_panel_size <- function(bed_path) {
  if (!file.exists(bed_path)) {
    warning("BED file not found: ", bed_path)
    return(NA_real_)
  }
  
  bed <- tryCatch({
    rtracklayer::import(bed_path)
  }, error = function(e) {
    warning("Could not parse BED file: ", e$message)
    return(NULL)
  })
  
  if (is.null(bed) || length(bed) == 0) {
    return(NA_real_)
  }
  
  # Calculate total covered bases (merge overlapping intervals first)
  gr <- GenomicRanges::reduce(bed)
  total_bp <- sum(as.numeric(width(gr)))
  total_mb <- total_bp / 1e6
  
  return(total_mb)
}

panel_size_mb <- load_panel_size(bed_path)

# =========================================================
# Helper: parse VCF and extract variant info for TMB
# =========================================================
extract_variants_for_tmb <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  
  vcf <- tryCatch({
    read.vcfR(path, verbose = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(vcf) || length(vcf@fix) == 0 || nrow(vcf@fix) == 0) {
    return(data.frame())
  }
  
  df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  
  # Extract basic info
  df$DP_INFO <- suppressWarnings(as.numeric(
    tryCatch(vcfR::extract.info(vcf, "DP"), error = function(e) rep(NA, nrow(df)))
  ))
  
  # Extract genotype AD for filtering
  ad <- tryCatch(vcfR::extract.gt(vcf, "AD")[, 1], error = function(e) rep(NA, nrow(df)))
  
  if (!all(is.na(ad))) {
    ad_split <- strsplit(ad, ",")
    df$AD_ALT <- suppressWarnings(as.numeric(sapply(ad_split, function(x) if (length(x) >= 2) x[2] else NA)))
    df$AD_REF <- suppressWarnings(as.numeric(sapply(ad_split, `[`, 1)))
  } else {
    df$AD_ALT <- NA
    df$AD_REF <- NA
  }
  
  # Get DP from FORMAT
  dp <- tryCatch(vcfR::extract.gt(vcf, "DP")[, 1], error = function(e) rep(NA, nrow(df)))
  df$DP <- suppressWarnings(as.numeric(dp))
  
  # Extract CSQ for IMPACT filtering
  csq_header <- grep("##INFO=<ID=CSQ", vcf@meta, value = TRUE)
  
  if (length(csq_header) > 0) {
    csq_header <- strsplit(csq_header, "Format: ")[[1]][2]
    csq_cols <- strsplit(csq_header, "\\|")[[1]]
    
    df$CSQ_raw <- sapply(df$INFO, function(x) {
      csq <- grep("^CSQ=", strsplit(x, ";")[[1]], value = TRUE)
      if (length(csq) == 0) return(NA_character_)
      sub("^CSQ=", "", csq[1])
    })
    
    csq_split <- strsplit(df$CSQ_raw, "\\|")
    csq_df <- do.call(rbind, lapply(csq_split, function(x) {
      if (length(x) < length(csq_cols)) {
        x <- c(x, rep(NA, length(csq_cols) - length(x)))
      }
      x[1:length(csq_cols)]
    }))
    colnames(csq_df) <- csq_cols
    csq_df <- as.data.frame(csq_df, stringsAsFactors = FALSE)
    
    df <- cbind(df, csq_df)
    df$CSQ_raw <- NULL
  }
  
  # Apply bcftools-like filters (simplified for TMB calculation)
  # Filter: PASS and (AD_ALT > 4 or DP > 80)
  df_filtered <- df %>%
    filter(
      FILTER == "PASS" | FILTER == ".",
      (AD_ALT > 4 | DP > 80)
    )
  
  # For TMB we want non-synonymous mutations
  # Keep only HIGH and MODERATE impact variants
  if ("IMPACT" %in% colnames(df_filtered)) {
    df_filtered <- df_filtered %>%
      filter(IMPACT %in% c("HIGH", "MODERATE"))
  }
  
  return(df_filtered)
}

# =========================================================
# Calculate TMB
# =========================================================
calculate_tmb <- function(variants_df, panel_mb) {
  if (is.na(panel_mb) || panel_mb <= 0) {
    return(NA_real_)
  }
  
  if (is.null(variants_df) || nrow(variants_df) == 0) {
    return(0)
  }
  
  # Count unique variants (by position)
  n_mutations <- nrow(variants_df)
  
  # TMB = mutations per Mb
  tmb <- n_mutations / panel_mb
  
  return(tmb)
}

# =========================================================
# Parse VCF for display table (with empty handling)
# =========================================================
make_vcf_table <- function(path, gene_filter = NULL) {
  
  # ---- file existence ----
  if (!file.exists(path)) {
    msg <- data.frame(Message = paste("VCF file not found:", path))
    return(datatable(msg, rownames = FALSE, options = list(dom = 't', ordering = FALSE)))
  }
  
  # ---- read VCF ----
  vcf <- tryCatch({
    read.vcfR(path, verbose = FALSE)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(vcf) || length(vcf@fix) == 0 || nrow(vcf@fix) == 0) {
    msg <- data.frame(Message = "VCF file is empty or could not be parsed.")
    return(datatable(msg, rownames = FALSE, options = list(dom = 't', ordering = FALSE)))
  }
  
  df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    msg <- data.frame(Message = "No variant records in VCF.")
    return(datatable(msg, rownames = FALSE, options = list(dom = 't', ordering = FALSE)))
  }
  
  # ---- extract INFO fields ----
  dp_info <- tryCatch(vcfR::extract.info(vcf, "DP"), error = function(e) rep(NA, nrow(df)))
  df$DP_INFO <- suppressWarnings(as.numeric(dp_info))
  
  # ---- extract FORMAT fields (first sample only) ----
  gt <- tryCatch(vcfR::extract.gt(vcf, "GT")[, 1], error = function(e) rep(NA, nrow(df)))
  gq <- tryCatch(vcfR::extract.gt(vcf, "GQ")[, 1], error = function(e) rep(NA, nrow(df)))
  dp <- tryCatch(vcfR::extract.gt(vcf, "DP")[, 1], error = function(e) rep(NA, nrow(df)))
  ad <- tryCatch(vcfR::extract.gt(vcf, "AD")[, 1], error = function(e) rep(NA, nrow(df)))
  vaf <- tryCatch(vcfR::extract.gt(vcf, "VAF")[, 1], error = function(e) rep(NA, nrow(df)))
  pl <- tryCatch(vcfR::extract.gt(vcf, "PL")[, 1], error = function(e) rep(NA, nrow(df)))
  
  df$GT <- gt
  df$GQ <- suppressWarnings(as.numeric(gq))
  df$DP <- suppressWarnings(as.numeric(dp))
  df$VAF <- suppressWarnings(as.numeric(vaf))
  df$PL <- pl
  
  # ---- AD split ----
  if (!all(is.na(ad))) {
    ad_split <- strsplit(ad, ",")
    df$AD_REF <- suppressWarnings(as.numeric(sapply(ad_split, `[`, 1)))
    df$AD_ALT <- suppressWarnings(as.numeric(sapply(ad_split, function(x) if (length(x) >= 2) x[2] else NA)))
  } else {
    df$AD_REF <- NA
    df$AD_ALT <- NA
  }
  
  # ---- HGVS ----
  df$HGVS <- paste0(df$CHROM, ":g.", df$POS, df$REF, ">", df$ALT)
  
  # ---- CSQ parsing (safe, keeps all rows) ----
  csq_header <- grep("##INFO=<ID=CSQ", vcf@meta, value = TRUE)
  
  if (length(csq_header) > 0) {
    csq_header <- strsplit(csq_header, "Format: ")[[1]][2]
    csq_cols <- strsplit(csq_header, "\\|")[[1]]
    
    # Extract first CSQ entry per variant (if multiple, take the first)
    df$CSQ_raw <- sapply(df$INFO, function(x) {
      csq <- grep("^CSQ=", strsplit(x, ";")[[1]], value = TRUE)
      if (length(csq) == 0) return(NA_character_)
      sub("^CSQ=", "", csq[1])
    })
    
    # Split into columns
    csq_split <- strsplit(df$CSQ_raw, "\\|")
    csq_df <- do.call(rbind, lapply(csq_split, function(x) {
      if (length(x) < length(csq_cols)) {
        x <- c(x, rep(NA, length(csq_cols) - length(x)))
      }
      x[1:length(csq_cols)]
    }))
    colnames(csq_df) <- csq_cols
    csq_df <- as.data.frame(csq_df, stringsAsFactors = FALSE)
    
    # Bind to original df (row order preserved)
    df <- cbind(df, csq_df)
    df$CSQ_raw <- NULL  # remove temporary column
  }
  
  # ---- group assignment ----
  if (!is.null(gene_filter) && "SYMBOL" %in% colnames(df)) {
    df$group <- ifelse(df$SYMBOL %in% gene_filter, "selected genes", "other genes")
  } else {
    df$group <- "all"
  }
  
  # ---- select relevant columns ----
  keep_cols <- c(
    "HGVS", "DP", "AD_REF", "AD_ALT", "VAF", "QUAL",
    "Consequence", "IMPACT", "SYMBOL", "Gene",
    "gnomADe_AF", "CLIN_SIG", "PUBMED", "group"
  )
  
  df_small <- df %>% select(any_of(keep_cols))
  
  if (nrow(df_small) == 0) {
    msg <- data.frame(Message = "No variants remain after filtering / parsing.")
    return(datatable(msg, rownames = FALSE, options = list(dom = 't', ordering = FALSE)))
  }
  
  # ---- return interactive table ----
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
# Extract variants for TMB calculation
# =========================================================
cat("Extracting variants for TMB calculation...\n")

variants_special <- extract_variants_for_tmb(vcf_path1)
variants_other <- extract_variants_for_tmb(vcf_path2)

# Combine all variants for total TMB
all_variants <- rbind(variants_special, variants_other)

# Calculate TMB
tmb_total <- calculate_tmb(all_variants, panel_size_mb)
tmb_special <- calculate_tmb(variants_special, panel_size_mb)
tmb_other <- calculate_tmb(variants_other, panel_size_mb)

# Count mutations
n_mutations_total <- nrow(all_variants)
n_mutations_special <- nrow(variants_special)
n_mutations_other <- nrow(variants_other)

# =========================================================
# Generate display tables
# =========================================================
cat("Generating tables...\n")
tbl1 <- make_vcf_table(vcf_path1, gene_list)
tbl2 <- make_vcf_table(vcf_path2, NULL)

# =========================================================
# UI
# =========================================================
ui <- tagList(
  
  h2("VCF Viewer with TMB Analysis"),
  
  # =========================================================
  # TMB Summary Panel
  # =========================================================
  tags$div(
    style = "padding:15px; background:#e8f4f8; border:2px solid #2196F3; border-radius:8px; margin-bottom:20px;",
    
    tags$h3("🧬 Tumor Mutational Burden (TMB) Summary"),
    
    tags$div(
      style = "display:flex; flex-wrap:wrap; gap:20px; margin:10px 0;",
      
      tags$div(
        style = "flex:1; min-width:150px; padding:10px; background:white; border-radius:5px;",
        tags$b("Panel Size:"),
        if (!is.na(panel_size_mb) && panel_size_mb > 0) {
          paste0(round(panel_size_mb, 2), " Mb")
        } else {
          "Not available"
        }
      ),
      
      tags$div(
        style = "flex:1; min-width:150px; padding:10px; background:white; border-radius:5px;",
        tags$b("Total TMB:"),
        if (!is.na(tmb_total)) {
          paste0(round(tmb_total, 2), " mutations/Mb")
        } else {
          "Not available"
        }
      ),
      
      tags$div(
        style = "flex:1; min-width:150px; padding:10px; background:white; border-radius:5px;",
        tags$b("Total Mutations:"),
        n_mutations_total
      )
    ),
    
    tags$hr(),
    
    tags$div(
      style = "display:flex; flex-wrap:wrap; gap:20px; margin:10px 0;",
      
      tags$div(
        style = "flex:1; min-width:150px; padding:10px; background:#e8f5e9; border-radius:5px;",
        tags$b("Selected genes TMB:"),
        if (!is.na(tmb_special)) {
          paste0(round(tmb_special, 2), " mutations/Mb")
        } else {
          "Not available"
        },
        tags$br(),
        tags$span(style = "font-size:0.9em; color:#666;",
                  paste0("(", n_mutations_special, " mutations)"))
      ),
      
      tags$div(
        style = "flex:1; min-width:150px; padding:10px; background:#fff3e0; border-radius:5px;",
        tags$b("Other genes TMB:"),
        if (!is.na(tmb_other)) {
          paste0(round(tmb_other, 2), " mutations/Mb")
        } else {
          "Not available"
        },
        tags$br(),
        tags$span(style = "font-size:0.9em; color:#666;",
                  paste0("(", n_mutations_other, " mutations)"))
      )
    ),
    
    tags$div(
      style = "margin-top:10px; font-size:0.9em; color:#555;",
      tags$i("TMB calculated using protein-coding variants (HIGH and MODERATE impact) with PASS filter and depth/allele fraction thresholds")
    )
  ),
  
  # =========================================================
  # Gene List and Filters
  # =========================================================
  tags$div(
    style = "padding:10px; background:#f5f5f5; margin-bottom:15px;",
    tags$b("Gene list (selected genes): "),
    gene_list_str,
    tags$br(),
    tags$b("BED file: "),
    basename(bed_path)
  ),
  
  tags$ul(
    class = "nav nav-tabs",
    tags$li(
      class = "active",
      tags$a(href = "#tab1", `data-toggle` = "tab", "Selected genes")
    ),
    tags$li(
      tags$a(href = "#tab2", `data-toggle` = "tab", "Other genes")
    )
  ),
  
  tags$div(class = "tab-content",
           
           # =========================
           # TAB 1 - SELECTED GENES
           # =========================
           tags$div(
             class = "tab-pane active",
             id = "tab1",
             
             tags$div(
               style = "
          border: 2px solid #4CAF50;
          border-radius: 8px;
          padding: 12px;
          margin-top: 10px;
          background: #f6fff6;
        ",
               
               tags$h3("🧬 Selected genes block"),
               
               tags$h4("Gene list"),
               tags$pre(gene_list_str),
               
               tags$h4("bcftools filter (selected)"),
               tags$pre(bcftools_filter_special),
               
               tags$h4("VEP filter (selected)"),
               tags$pre(vep_filter_special),
               
               tags$hr(),
               
               tags$h4(paste0("Variants table (selected genes) - ", n_mutations_special, " mutations")),
               tbl1
             )
           ),
           
           # =========================
           # TAB 2 - OTHER GENES
           # =========================
           tags$div(
             class = "tab-pane",
             id = "tab2",
             
             tags$div(
               style = "
          border: 2px solid #FF9800;
          border-radius: 8px;
          padding: 12px;
          margin-top: 10px;
          background: #fffaf3;
        ",
               
               tags$h3("🧬 Other genes block"),
               
               tags$h4("Gene list (implicit)"),
               tags$div(
                 style = "color:#666;",
                 "Everything NOT in selected gene list"
               ),
               
               tags$h4("bcftools filter (other)"),
               tags$pre(bcftools_filter_other),
               
               tags$h4("VEP filter (other)"),
               tags$pre(vep_filter_other),
               
               tags$hr(),
               
               tags$h4(paste0("Variants table (other genes) - ", n_mutations_other, " mutations")),
               tbl2
             )
           )
  )
)

# =========================================================
# Save HTML
# =========================================================
cat("Saving HTML report...\n")
htmltools::save_html(
  ui,
  file = out_html
)

cat("✅ Saved:", out_html, "\n")
cat("📊 TMB Summary:\n")
cat("   Total TMB:", round(tmb_total, 2), "mutations/Mb\n")
cat("   Selected genes TMB:", round(tmb_special, 2), "mutations/Mb\n")
cat("   Other genes TMB:", round(tmb_other, 2), "mutations/Mb\n")
cat("   Total mutations:", n_mutations_total, "\n")