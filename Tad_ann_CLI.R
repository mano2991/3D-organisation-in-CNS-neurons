#!/usr/bin/env Rscript

# ==============================================================================
# Tad_ann_CLI.R - Manual Args Version
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(tidyr)
})

# --- Manual Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && length(args) >= idx + 1) {
    return(args[idx + 1])
  }
  return(default)
}

mode <- get_arg("--mode", "both")
n1   <- get_arg("--n1", get_arg("--name1"))
n2   <- get_arg("--n2", get_arg("--name2"))
t1   <- get_arg("--t1", get_arg("--tad1"))
t2   <- get_arg("--t2", get_arg("--tad2"))
l1   <- get_arg("--l1", get_arg("--loop1"))
l2   <- get_arg("--l2", get_arg("--loop2"))
gene_file <- get_arg("--gene", "gene_split.csv")
out_dir   <- get_arg("--out")

if (is.null(n1) || is.null(n2)) {
  cat("Usage: Tad_ann_CLI.R --n1 <name1> --n2 <name2> [options]\n")
  cat("Options: --mode, --l1, --l2, --t1, --t2, --gene, --out\n")
  q(status = 1)
}

if (is.null(out_dir)) out_dir <- paste0(n1, "_vs_", n2, "_comparison")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper Functions ---

annotate_with_genes <- function(bed_df, gene_gr) {
  ann_anchor <- function(chr, start, end) {
    gr <- GRanges(chr, IRanges(start, end))
    ol <- findOverlaps(gr, gene_gr)
    res <- rep(NA_character_, length(gr))
    if (length(ol) > 0) {
      subject_hits <- subjectHits(ol)
      query_hits <- queryHits(ol)
      gene_names <- gene_gr$gene_name[subject_hits]
      genes_pasted <- tapply(gene_names, query_hits, function(x) paste(unique(x), collapse=","))
      res[as.numeric(names(genes_pasted))] <- as.character(genes_pasted)
    }
    return(res)
  }
  if (ncol(bed_df) >= 6) {
    bed_df$gene_a1 <- ann_anchor(bed_df[[1]], bed_df[[2]], bed_df[[3]])
    bed_df$gene_a2 <- ann_anchor(bed_df[[4]], bed_df[[5]], bed_df[[6]])
  } else {
    bed_df$gene_id <- ann_anchor(bed_df[[1]], bed_df[[2]], bed_df[[3]])
  }
  return(bed_df)
}

calculate_overlap_val <- function(s1, e1, s2, e2) {
  inter_s <- pmax(s1, s2)
  inter_e <- pmin(e1, e2)
  inter_w <- pmax(0, inter_e - inter_s)
  avg_w   <- ((e1 - s1) + (e2 - s2)) / 2
  return(round((inter_w / avg_w) * 100, 2))
}

# --- Main Execution ---
cat("Loading Gene Annotations from:", gene_file, "\n")
gene_data <- read.csv(gene_file)
gene_gr <- GRanges(seqnames = gene_data$Chromosome, 
                   ranges = IRanges(start = gene_data$Start, end = gene_data$End), 
                   gene_name = gene_data$Gene_name)

if (mode %in% c("loop", "both") && !is.null(l1) && !is.null(l2)) {
  cat("\n=== Running Loop Analysis ===\n")
  LOOP_DIR <- file.path(out_dir, "Loop")
  dir.create(LOOP_DIR, showWarnings = FALSE)
  
  l1_df <- read.table(l1)
  l2_df <- read.table(l2)
  
  colnames(l1_df)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
  colnames(l2_df)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
  
  l1_df$id <- seq_len(nrow(l1_df))
  l2_df$id <- seq_len(nrow(l2_df))
  
  l1_ann <- annotate_with_genes(l1_df, gene_gr)
  l2_ann <- annotate_with_genes(l2_df, gene_gr)
  
  gr1 <- GRanges(l1_ann$chr1, IRanges(pmin(l1_ann$start1, l1_ann$start2), pmax(l1_ann$end1, l1_ann$end2)))
  gr2 <- GRanges(l2_ann$chr1, IRanges(pmin(l2_ann$start1, l2_ann$start2), pmax(l2_ann$end1, l2_ann$end2)))
  
  ol <- findOverlaps(gr1, gr2)
  
  if (length(ol) > 0) {
    ov_df <- data.frame(q_idx = queryHits(ol), s_idx = subjectHits(ol))
    ov_df$overlap <- mapply(function(q, s) {
      calculate_overlap_val(start(gr1)[q], end(gr1)[q], start(gr2)[s], end(gr2)[s])
    }, ov_df$q_idx, ov_df$s_idx)
    
    # --- Greedy 1-to-1 Matching ---
    potential_matches <- ov_df %>% 
      dplyr::filter(overlap >= 40) %>%
      dplyr::arrange(desc(overlap))
    
    used_q <- logical(nrow(l1_ann))
    used_s <- logical(nrow(l2_ann))
    matched_q <- integer(0); matched_s <- integer(0); matched_ov <- numeric(0)
    
    if (nrow(potential_matches) > 0) {
      for (i in seq_len(nrow(potential_matches))) {
        q_curr <- potential_matches$q_idx[i]
        s_curr <- potential_matches$s_idx[i]
        if (!used_q[q_curr] && !used_s[s_curr]) {
          matched_q <- c(matched_q, q_curr); matched_s <- c(matched_s, s_curr)
          matched_ov <- c(matched_ov, potential_matches$overlap[i])
          used_q[q_curr] <- TRUE; used_s[s_curr] <- TRUE
        }
      }
    }
    best_ov <- data.frame(q_idx = matched_q, s_idx = matched_s, overlap = matched_ov)
    
    l1_res <- l1_ann %>% dplyr::left_join(best_ov, by = c("id" = "q_idx")) %>%
      dplyr::mutate(status = dplyr::case_when(
        is.na(s_idx) ~ paste0("Specific_", n1),
        overlap >= 80 ~ "Conserved",
        TRUE ~ "Shifted"
      ))
    
    l2_match <- l2_ann
    colnames(l2_match) <- paste0(n2, "_", colnames(l2_match))
    
    final_df <- l1_res %>% 
      dplyr::left_join(l2_match, by = setNames(paste0(n2, "_id"), "s_idx"), keep = TRUE) %>%
      dplyr::select(-s_idx)
    
    s1_cols <- setdiff(colnames(l1_ann), "id")
    final_df <- final_df %>% dplyr::rename_with(~paste0(n1, "_", .), dplyr::all_of(s1_cols))
    
    matched_s2 <- unique(na.omit(best_ov$s_idx))
    l2_spec <- l2_ann %>% dplyr::filter(!(id %in% matched_s2))
    
    if (nrow(l2_spec) > 0) {
      l2_spec_prep <- l2_spec
      colnames(l2_spec_prep) <- paste0(n2, "_", colnames(l2_spec_prep))
      l2_spec_prep$status <- paste0("Specific_", n2)
      l2_spec_prep$overlap <- 0
      final_df <- dplyr::bind_rows(final_df, l2_spec_prep)
    }
    
    final_df <- final_df %>% 
      dplyr::select(-dplyr::any_of(c(paste0(n1, "_id"), paste0(n2, "_id")))) %>%
      dplyr::select(status, overlap, dplyr::everything())

  } else {
    final_df <- l1_ann %>% mutate(status = paste0("Specific_", n1), overlap = 0)
    colnames(final_df) <- paste0(n1, "_", colnames(final_df))
    l2_spec <- l2_ann
    colnames(l2_spec) <- paste0(n2, "_", colnames(l2_spec))
    l2_spec$status <- paste0("Specific_", n2)
    l2_spec$overlap <- 0
    final_df <- bind_rows(final_df, l2_spec)
  }
  
  write.csv(final_df, file.path(LOOP_DIR, "Loop_comparison_final.csv"), row.names = FALSE)
  
  p <- ggplot(final_df, aes(x = status, fill = status)) + geom_bar() + 
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5) +
    theme_classic() + labs(title = paste("Loop Dynamics:", n1, "vs", n2))
  ggsave(file.path(LOOP_DIR, "Loop_Summary_Plot.png"), p, width = 8, height = 6)
}

cat("\nAnalysis Complete. Results in:", out_dir, "\n")
