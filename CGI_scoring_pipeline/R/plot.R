######
# PLOTTING CODE
######

#' Plot replicate comparisons.
#' 
#' Plots replicate comparisons for all replicates in a list of screens and outputs
#' plots to a given folder.
#' 
#' @param df Reads or lfc dataframe.
#' @param screens List of screens created with \code{add_screens}.
#' @param output_folder Folder to output plots to. 
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @param display_numbers Whether or not to include PCC values in heatmap (default TRUE).
#' @param show_rownames Whether or not to show row names on the plot (default TRUE).
#' @param show_colnames Whether or not to show column names on the plot (default TRUE).
#' @export 
plot_lfc_qc <- function(df, screens, output_folder, plot_type = "png",
                        display_numbers = TRUE, show_rownames = TRUE,
                        show_colnames = TRUE) {
  
  # Checks for input errors
  check_screen_params(df, screens)
  
  # Gets essential gene QC metrics
  auc <- essential_lfc_qc(df, screens)
  auc_file <- file.path(output_folder, "essential_PR_QC.tsv")
  if (!is.null(auc)) {
    utils::write.table(auc, auc_file, quote = FALSE, sep = "\t",
                       row.names = FALSE, col.names = TRUE) 
  }
  
  # Checks plot type and converts to lowercase
  plot_type <- tolower(plot_type)
  if (plot_type != "png" & plot_type != "pdf") {
    stop("plot_type must be either png or pdf")
  }
  
  
  # Gets sample groups
  all_cols <- c()
  col_groups <- c()
  i <- 1
  for (screen_name in names(screens)) {
    screen <- screens[[screen_name]]
    for (col in screen[["replicates"]]) {
      all_cols <- c(all_cols, col)
      col_groups[i] <- screen_name
      i <- i +1
    }
  }
  
  # Plots heatmap of log-scaled read counts
  df <- df[,all_cols]
  heatmap_file <- file.path(output_folder, paste0("lfc_heatmap.", plot_type))
  plot_heatmap(df, col_groups, heatmap_file, display_numbers,
               show_rownames, show_colnames)
  
  # Builds dataframe of replicate PCCs
  cor_df <- NULL
  
  # Compares replicates across all screens
  for (screen_name in names(screens)) {
    screen <- screens[[screen_name]]
    rep_cols <- screen[["replicates"]]
    if (length(rep_cols) > 1) {
      pairs <- utils::combn(rep_cols, 2)
      for (i in 1:ncol(pairs)) {
        col1 <- pairs[1,i]
        col2 <- pairs[2,i]
        x_label <- paste0(col1, " log fold change")
        y_label <- paste0(col2, " log fold change")
        temp <- plot_samples(df, col1, col2, x_label, y_label, print_cor = TRUE)
        p <- temp[[1]]
        pcc <- temp[[2]]
        scc <- temp[[3]]
        file_name <- paste0(col1, "_vs_", col2, "_replicate_comparison.", plot_type)
        suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300))
        
        # Stores PCC in dataframe
        if (is.null(cor_df)) {
          cor_df <- data.frame(screen = screen_name, rep1 = col1, rep2 = col2, pcc = pcc, scc = scc,
                               stringsAsFactors = FALSE)
        } else {
          cor_df <- rbind(cor_df, c(screen_name, col1, col2, pcc, scc))
        }
      }   
    }
  }
  
  # Writes PCCs to file
  cor_file <- file.path(output_folder, "replicate_cor.tsv")
  utils::write.table(cor_df, cor_file, quote = FALSE, sep = "\t",
                     row.names = FALSE, col.names = TRUE)
}

#' Plot read counts for a screen.
#' 
#' Plots a histogram of read counts for each replicate of all screens. Also
#' plots total reads for all screens.
#' 
#' @param df Reads dataframe.
#' @param screens List of screens created with \code{add_screens}.
#' @param output_folder Folder to output plots to. 
#' @param log_scale If true, log-normalizes data.
#' @param pseudocount Pseudocounts to add to log-normalized data if specified (default 1).
#' @param display_numbers Whether or not to include PCC values in heatmap (default TRUE).
#' @param show_rownames Whether or not to show row names on the plot (default TRUE).
#' @param show_colnames Whether or not to show column names on the plot (default TRUE).
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
plot_reads_qc <- function(df, screens, output_folder, plot_type = "png",
                          log_scale = TRUE, pseudocount = 1,
                          display_numbers = TRUE, show_rownames = TRUE,
                          show_colnames = TRUE) {
  
  # Checks for input errors
  check_screen_params(df, screens)
  
  # Checks plot type and converts to lowercase
  plot_type <- tolower(plot_type)
  if (plot_type != "png" & plot_type != "pdf") {
    stop("plot_type must be either png or pdf")
  }
  
  # Plots read count histograms for all replicates of all screens and stores total reads 
  reads_df <- NULL
  all_cols <- c()
  col_groups <- c()
  all_coverage <- c()
  i <- 1
  for (screen_name in names(screens)) {
    screen <- screens[[screen_name]]
    for (col in screen[["replicates"]]) {
      total_reads <- sum(df[,col], na.rm = TRUE)
      if (is.null(reads_df)) {
        reads_df <- data.frame(rep = col, reads = total_reads,
                               stringsAsFactors = FALSE)
      } else {
        reads_df <- rbind(reads_df, c(col, total_reads))
      }
      all_coverage <- c(all_coverage, screen[["target_coverage"]])
      all_cols <- c(all_cols, col)
      p <- plot_reads(df, col, log_scale, pseudocount)
      file_name <- paste0(col, "_raw_reads_histogram.", plot_type)
      suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300))
      col_groups[i] <- screen_name
      i <- i + 1
    }
  }
  
  # Gets unique coverage breakpoints
  all_coverage <- unique(all_coverage)
  all_coverage <- all_coverage*nrow(df)
  
  # Plots total reads
  reads_df$reads <- as.numeric(reads_df$reads)
  p <- ggplot2::ggplot(reads_df, ggplot2::aes_string(x = "rep", y = "reads")) +
    ggplot2::geom_bar(stat = "identity", color = "Black", fill = "gray30")
  for (coverage in all_coverage) {
    p <- p + ggplot2::geom_hline(yintercept = coverage, linetype = 2, size = 1, alpha = 0.9, color = "Gray")
  }
  p <- p +
    ggplot2::xlab("Replicate") +
    ggplot2::ylab("Total reads") +
    ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    ggthemes::theme_tufte(base_size = 20) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  file_name <- paste0("total_reads.", plot_type)
  suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), plot = p, width = 10, height = 7, dpi = 300))
  
  # Log-scales read counts if specified
  df <- df[,all_cols]
  if (log_scale) {
    for (col in colnames(df)) {
      df[,col] <- log2(df[,col] + 1)
    }
  }
  
  # Plots heatmap of log-scaled read counts
  heatmap_file <- file.path(output_folder, paste0("reads_heatmap.", plot_type))
  plot_heatmap(df, col_groups, heatmap_file, display_numbers,
               show_rownames, show_colnames)
}

#' Plot read counts.
#'
#' Plots a histogram of read counts for a given column.
#'
#' @param df Reads dataframe.
#' @param col Name of column to plot.
#' @param log_scale If true, log-normalizes data.
#' @param pseudocount Pseudocounts to add to log-normalized data if specified (default 1).
#' @return A ggplot object.
plot_reads <- function(df, col, log_scale = TRUE, pseudocount = 1) {
  x_label <- paste(col, "log-normalized read counts")
  y_label <- "Number of read counts"
  if (log_scale) {
    df[,col] <- log2(df[,col] + 1)
    y_label <- "Number of log-normalized read counts"
  }
  p <- ggplot2::ggplot(df, ggplot2::aes_string(col)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::xlab(x_label) +
    ggplot2::ylab(y_label) +
    ggthemes::theme_tufte(base_size = 20)
  return(p)
}

#' Plots sample comparisons.
#'
#' Pretty-plots comparisons between two samples in a scatterplot.
#
#' @param df Reads or lfc dataframe.
#' @param xcol Name of column containing values to plot on the x-axis.
#' @param ycol Name of column containing values to plot on the y-axis.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param color_col Name of column to color points by (optional).
#' @param color_lab Name of color legend (optional, defaults to color_col).
#' @param print_cor If true, prints Pearson correlation between columns 
#'   (default FALSE).
#' @return A list of three elements. The first is a ggplot object, the
#'   second is the correlation between xcol and ycol, and the third is
#'   the Spearman correlation between xcol and ycol.
plot_samples <- function(df, xcol, ycol, xlab, ylab, 
                         color_col = NULL, color_lab = NULL,
                         print_cor = FALSE) {
  
  # Computes correlations and optionally prints Pearson correlation
  pcc <- NA
  scc <- NA
  if (sum(stats::complete.cases(df[,colnames(df) %in% c(xcol, ycol)])) > 10) {
    pcc <- stats::cor(df[[xcol]], df[[ycol]], use = "complete.obs")
    scc <- stats::cor(df[[xcol]], df[[ycol]], method = "spearman", use = "complete.obs")
    if (print_cor) {
      cat(paste("Pearson correlation between", xcol, "and", ycol, ":", pcc, "\n"))
    }
  } else {
    warning(paste("too few complete element pairs to take correlations for", xcol, "and", ycol))
  }
  
  # Makes plot
  p <- NULL
  if (is.null(color_col)) {
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = xcol, y = ycol)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2, size = 1) +
      ggplot2::geom_point(size = 1.5, alpha = 0.7) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggthemes::theme_tufte(base_size = 20) 
  } else {
    if (is.null(color_lab)) {
      color_lab <- color_col
    }
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = xcol, y = ycol, color = color_col)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2, size = 1) +
      ggplot2::geom_point(size = 1.5, alpha = 0.7) +
      ggplot2::scale_color_gradientn(colors = c("blue", "gray"), name = color_lab) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggthemes::theme_tufte(base_size = 20)
  }
  return(list(p, pcc, scc))
}

#' Plots Pearson correlation heatmap.
#' 
#' Plots heatmap of Pearson correlations for a dataframe with labels
#' for groups of samples and a blue-yellow color scheme.
#' 
#' @param df Reads or lfc dataframe
#' @param col_groups List of grouping labels for each column in df. 
#' @param filename Output filename for plot. 
#' @param display_numbers Whether or not to include PCC values in heatmap (default TRUE).
#' @param show_rownames Whether or not to show row names on the plot (default TRUE).
#' @param show_colnames Whether or not to show column names on the plot (default TRUE).
#' @return Writes plot to file with no return value.
plot_heatmap <- function(df, col_groups, filename, display_numbers = TRUE, 
                         show_rownames = TRUE, show_colnames = TRUE) {
  
  # Returns warning if not enough complete observations
  if (sum(stats::complete.cases(df)) <= 10) {
    warning(paste("too few complete observations to construct heatmap"))
    return()
  }
  
  # Gets colors for different screens
  screen_colors <- NA
  n_colors <- length(unique(col_groups))
  if (n_colors < 10) {
    screen_colors <- list(group = RColorBrewer::brewer.pal(length(unique(col_groups)), "Set1"))
  } else {
    pal <- RColorBrewer::brewer.pal(9, "Set1")
    pal <- grDevices::colorRampPalette(pal)(n_colors)
    screen_colors <- list(group = pal)
  }
  names(screen_colors$group) <- unique(col_groups)
  
  # Gets PCCs for heatmap
  cor_mat <- data.matrix(stats::cor(df, use = "complete.obs"))
  
  # Gets annotation for heatmap
  col_groups <- data.frame("Screen" = col_groups)
  rownames(col_groups) <- colnames(df)
  colnames(cor_mat) <- colnames(df)
  
  # Gets color for heatmap values
  breaks <- seq(-1, 1, by = (1/150))
  pal <- grDevices::colorRampPalette(c("#7fbf7b", "#f7f7f7", "#af8dc3"))(n = length(breaks))
  
  # Plots heatmap of raw reads to file
  pheatmap::pheatmap(cor_mat,
                     border_color = NA,
                     annotation_col = col_groups,
                     annotation_colors = screen_colors,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     display_numbers = display_numbers,
                     color = pal, 
                     breaks = breaks,
                     filename = filename)
}


#' Plots chromosomal QC. 
#' 
#' Plots qGI effects by chromosomal location.
#' 
#' @param df LFC dataframe passed into \code{score_drugs_vs_controls} with chromosomal location 
#'   specified in a column named "chr", start coordinates specified in "start_loc" and stop 
#'   coordinates specified in "stop_loc." 
#' @param guide_df Dataframe of guide-level qGI scores. 
#' @param output_folder Folder to output plots to. 
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
plot_chromosome <- function(df, guide_df, output_folder, plot_type = "png") {
  
  # Gets guide density for each relevant chromosome
  chrom <- names(sort(table(df$chr)))
  chr_features <- data.frame(matrix(nrow = length(chrom), ncol = 4))
  rownames(chr_features) <- chrom
  colnames(chr_features) <- c("length", "n_genes", "length_per", "roll_window")
  for (i in 1:length(chrom)) {
    chr_features[i,"length"] <- max(df$stop_loc[df$chr == chrom[i]]) - min(df$start_loc[df$chr == chrom[i]])
    chr_features[i,"n_genes"] <- length(which(df$chr == chrom[i]))
  }
  chr_features$length_per <- chr_features$n_genes / chr_features$length
  chr_features$roll_window <- ceiling(log2(chr_features$length_per / stats::median(chr_features$length_per) + 2) * 150)
  chr_features$labels <- chrom
  
  # Plots chromosomal features information
  file_name <- paste0("chromosome_genes.", plot_type)
  p <- ggplot2::ggplot(chr_features, ggplot2::aes_string(x = "length", y = "n_genes")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Chromosome length") +
    ggplot2::ylab("Number of genes") +
    ggplot2::geom_text(ggplot2::aes_string(label = "labels"), hjust = -0.1, vjust = -0.1) +
    ggthemes::theme_tufte(base_size = 20)
  suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), 
                                   plot = p, width = 10, height = 7, dpi = 300))
  file_name <- paste0("chromosome_adjusted_genes.", plot_type)
  p <- ggplot2::ggplot(chr_features, ggplot2::aes_string(x = "length", y = "length_per")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Chromosome length") +
    ggplot2::ylab("Number of guides adjusted by chromosomal length") +
    ggplot2::geom_text(ggplot2::aes_string(label = "labels"), hjust = -0.1, vjust = -0.1) +
    ggthemes::theme_tufte(base_size = 20)
  suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), 
                                   plot = p, width = 10, height = 7, dpi = 300))
  file_name <- paste0("chromosome_rolling_window.", plot_type)
  p <- ggplot2::ggplot(chr_features, ggplot2::aes_string(x = "length", y = "roll_window")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Chromosome length") +
    ggplot2::ylab("Width of rolling window") +
    ggplot2::geom_text(ggplot2::aes_string(label = "labels"), hjust = -0.1, vjust = -0.1) +
    ggthemes::theme_tufte(base_size = 20)
  suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), 
                                   plot = p, width = 10, height = 7, dpi = 300))
  
  # Ignores chromosomes with too few genes for the rolling window
  to_keep <- rep(TRUE, nrow(chr_features))
  for (i in 1:nrow(chr_features)) {
    if (chr_features[i,"n_genes"] < chr_features[i,"roll_window"]) {
      to_keep[i] <- FALSE
    }
  }
  chr_features <- chr_features[to_keep,]
  
  # Identifies windows where chromosomal shifts occur
  noise_factor <- apply(guide_df[,2:ncol(guide_df)], 2, stats::sd, na.rm = TRUE)
  chr_shifted_genes <- list()
  for (i in 1:nrow(chr_features)) {
    chrom <- row.names(chr_features)[i]
    for (j in 2:ncol(guide_df)) {
      condition <- colnames(guide_df)[j]
      
      # ???
      th1_match <- .09 + noise_factor[condition] * .2
      th2_match <- .03 + noise_factor[condition] * .1
      th3_match <- .04 + noise_factor[condition] * .2
      th4_match <- .03 + noise_factor[condition] * .2
      chrom_guides <- which(df$chr %in% chrom)
      x <- df$start_loc[chrom_guides]
      y <- guide_df[chrom_guides, condition]
      y <- y[order(x)]
      x <- x[order(x)]
      
      # Takes a running mean across genomic coordinates, with windows
      # sized relative to the gene density for the chromosome
      rollwindow <- chr_features[chrom, "roll_window"]
      rmean <- rep(NA, (length(x) - rollwindow))
      for (k in 1:length(rmean)) {
        y_rw <- y[k:(rollwindow + k - 1)]
        y_q <- stats::quantile(y_rw, probs = c(.05, .95), na.rm = TRUE)
        y_rw <- y_rw[y_rw > y_q[[1]] & y_rw < y_q[[2]]]
        rmean[k] <- mean(y_rw, na.rm = TRUE)
      }
      
      # 2. d running mean
      d_rmean <- rmean[1:(length(rmean)-rollwindow)] - rmean[(1+rollwindow):length(rmean)]
      
      # 3.2 maxima
      counter <- 0
      x <- 1
      
      while(length(x) > counter) {
        counter <- counter + 1
        if (counter == 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, maxORmin = "max", pi = y)
        } else if (counter > 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, tb = x, maxORmin = "max", pi = y)
        }
      }
      x_maxima <- x
      
      # 3.1 minima
      counter <- 0
      x <- 1
      
      while(length(x) > counter) {
        counter <- counter + 1
        if (counter == 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, maxORmin = "min", pi = y)
        } else if (counter > 1) {
          x <- fragment_transition(d_rmean, rmean, th1_match, th2_match, tb = x, maxORmin = "min", pi = y)
        }
      }
      x_minima <- x
      
      # Identifies all shifted stretches
      chr_shifted_genes <- define_fragments(chr_shifted_genes, 
                                            x_max = x_maxima, x_min = x_minima,
                                            th3 = th3_match, th4 = th4_match,
                                            chromOI = chrom, condition = condition, 
                                            x_ref = rmean, pi = y,
                                            b = rollwindow, chrAnno = names(y))
      
      # # Plots shifted stretches
      # x <- df$start_loc[chrom_guides]
      # y <- guide_df[chrom_guides, condition]
      # 
      # plot(y[order(x)], 
      #      ylim = range(c(y, rmean, d_rmean), na.rm = TRUE), 
      #      cex = 0.4, pch = 16,
      #      col = ifelse(names(y[order(x)]) %in% unlist(chr_shifted_genes[[condition]][[chrom]]), "#4daf4a", "grey"),
      #      main = condition, bty = "n", 
      #      xlab = paste("Guide RNA on", chrom), 
      #      ylab = "Guide RNA Pi-score")
      # abline(h = c(0, th1_match, -th1_match), lty=c(2,3,3), col = c("black", "blue", "blue"))
      # 
      # for (i in 1:length(chr_shifted_genes[[condition]][[chrom]])) {
      #   abline(h = mean(y[chr_shifted_genes[[condition]][[chrom]][[i]]], na.rm = TRUE), lty = 2, col="orange")
      # }
      # 
      # lines(rmean, lwd = 1)
      # lines(d_rmean, lwd = 1, col = "blue")
      # abline(v = c(x_minima, x_minima+rollwindow, x_maxima, x_maxima+rollwindow), lty = 3, col = "red")
      # 
      # plot(x, y, ylim = range(c(y, rmean, d_rmean), na.rm = TRUE), cex = 0.4, pch = 16,
      #      col = ifelse(names(y) %in% unlist(chr_shifted_genes[[condition]][[chrom]]), "#4daf4a", "grey"),
      #      main = condition, bty = "n", xlab = paste("Location of guide RNA on", chrom), ylab = "Guide RNA Pi-score")
      # abline(h = 0, lty = 2)
      # 
      # for (i in 1:length(chr_shifted_genes[[condition]][[chrom]])) {
      #   abline(h = mean(y[chr_shifted_genes[[condition]][[chrom]][[i]]], na.rm = TRUE), lty = 2, col = "orange")
      # }
      # par(par.shift)
      # dev.off()
    }
  }
  
  # Applies chromosomal correction to shifted genes if any exist
  if (length(unlist(chr_shifted_genes)) > 0) {
    for (condition in colnames(guide_df)[2:dim(guide_df)]) {
      for (chrom in names(chr_shifted_genes[[condition]])) {
        if (length(unlist(chrom)) > 0) {
          for (gene in chr_shifted_genes[[condition]][[chrom]]) {
            goi <- chr_shifted_genes[[condition]][[chrom]][[gene]]
            goi <- which(guide_df$gene %in% goi) 
            guide_df[goi, condition] <- df[goi, condition] - mean(df[goi, condition], na.rm = TRUE)
          }
        }
      }
    } 
  }
  
  # Returns corrected scores
  return(guide_df)
}


#' Plots drug response for scored data.
#' 
#' Pretty-plots response for chemogenomic (e.g. for directly comparing drug 
#' response to DMSO response). Assumes that data was scored by 
#' \code{score_drugs_vs_control} and significant effects were called by 
#' \code{call_drug_hits}.
#' 
#' @param scores Dataframe of scores returned from \code{call_drug_hits}.
#' @param control_name Name of control passed to \code{call_drug_hits}, or 
#'   NULL for data scored by \code{score_drugs_vs_controls} (default NULL).
#' @param condition_name Name of condition passed to \code{call_drug_hits}.
#' @param output_folder Folder to output plots to. 
#' @param loess If true and data was loess-normalized, plots loess null model instead
#'   (default TRUE).
#' @param neg_type Label for significant effects with a negative differential effect
#'   passed to \code{call_drug_hits} (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   passed to \code{call_drug_hits} (default "Positive").
#' @param label_fdr_threshold Threshold below which to plot gene labels for significant
#'   hits, or NULL to plot without labels (default NULL).
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
plot_drug_response <- function(scores, control_name = NULL, 
                               condition_name = NULL, output_folder = NULL,
                               loess = TRUE, neg_type = "Negative", 
                               pos_type = "Positive", label_fdr_threshold = NULL,
                               plot_type = "png") {
  
  # Gets variables depending on scoring type
  plot_file <- paste0(condition_name, "_vs_", control_name, "_scatter.", plot_type)
  control_mean_col <- paste0("mean_", control_name)
  condition_mean_col <- paste0("mean_", condition_name)
  response_col <- paste0("effect_type_", condition_name)
  diff_col <- paste0("differential_", condition_name, "_vs_", control_name) 
  fdr_col <- paste0("fdr_", condition_name, "_vs_", control_name) 
  control_label <- paste0(control_name, " LFC")
  if (is.null(control_name)) {
    plot_file <- paste0(condition_name, "_vs_controls_scatter.", plot_type)
    control_mean_col <- paste0("mean_controls_", condition_name)
    diff_col <- paste0("differential_", condition_name, "_vs_controls") 
    fdr_col <- paste0("fdr_", condition_name, "_vs_controls")
    control_label <- paste0("Weighted control LFC")
  }
  
  # Adds backticks to column names if they contain special characters
  control_mean_col_aes <- paste0("`", control_mean_col, "`")
  condition_mean_col_aes <- paste0("`", condition_mean_col, "`")
  response_col_aes <- paste0("`", response_col, "`")
  diff_col_aes <- paste0("`", diff_col, "`")
  fdr_col_aes <- paste0("`", fdr_col, "`")
  response_col_aes <- paste0("`", response_col, "`")
  control_label_aes <- paste0("`", control_label, "`")
  
  # Sets factors to plot significant effects last
  scores$sort_col <- abs(scores[[diff_col]])
  scores <- dplyr::arrange(scores, "sort_col")
  scores[[response_col]] <- forcats::fct_inorder(scores[[response_col]])
  
  # Gets subset dataframe for plotting labels
  subset_scores <- NULL
  if (!is.null(label_fdr_threshold)) {
    subset_scores <- scores[scores[[fdr_col]] < label_fdr_threshold &
                              scores[[response_col]] != "None",] 
  }
  
  # Manually sets colors for plot
  colors <- c("None" = "Gray", neg_type = "Black", pos_type = "Black")
  names(colors)[2] <- neg_type
  names(colors)[3] <- pos_type
  fill <- c("None" = "Gray", neg_type = "Blue", pos_type = "Yellow")
  names(fill)[2] <- neg_type
  names(fill)[3] <- pos_type

  # Builds basic plot
  p <- ggplot2::ggplot(scores, ggplot2::aes_string(x = control_mean_col_aes, 
                                                   y = condition_mean_col_aes)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray")
  
  # Appends choice of null model to plot
  if (loess) {
  } else {
    p <- p + 
      ggplot2::geom_abline(slope = 1, intercept = 0, size = 1.5, alpha = 0.5, color = "Black")
  }
  
  # Appends each layer to plot individually
  point_levels <- levels(scores[[response_col]])
  layer1 <- scores[scores[[response_col]] == point_levels[1],]
  layer2 <- scores[scores[[response_col]] == point_levels[2],]
  layer3 <- scores[scores[[response_col]] == point_levels[3],]
  p <- p +
    ggplot2::geom_point(data = layer1, 
                        ggplot2::aes_string(color = response_col_aes, fill = response_col_aes), 
                        shape = 21, alpha = 0.7) +
    ggplot2::geom_point(data = layer2, 
                        ggplot2::aes_string(color = response_col_aes, fill = response_col_aes), 
                        shape = 21, alpha = 0.7) +
    ggplot2::geom_point(data = layer3, 
                        ggplot2::aes_string(color = response_col_aes, fill = response_col_aes), 
                        shape = 21, alpha = 0.7)
  
  # Appends labels to plot
  if (!is.null(label_fdr_threshold)) {
    p <- p +
      ggrepel::geom_text_repel(data = subset_scores,
                               ggplot2::aes_string(x = control_mean_col_aes, 
                                                   y = condition_mean_col_aes,
                                                   label = "gene"),
                               color = "black")
  }

  # Finishes plot
  p <- p + 
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = fill) +
    ggplot2::xlab(control_label) +
    ggplot2::ylab(paste0(condition_name, " LFC")) +
    ggplot2::labs(fill = "Significant response") +
    ggplot2::guides(color = "none", size = "none") +
    ggthemes::theme_tufte(base_size = 20) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = "Black", size = 16),
                   axis.text.y = ggplot2::element_text(color = "Black", size = 16),
                   legend.text = ggplot2::element_text(size = 16))
  
  # Saves to file
  suppressWarnings(ggplot2::ggsave(file.path(output_folder, plot_file), 
                                   width = 10, height = 7, dpi = 300))
}

#' Plot guide-level residuals for all hits
#' 
#' Plots guide-level residuals for each called hit and outputs plots to a given folder. 
#' Works for data returned from \code{call_drug_hits}.
#' 
#' @param scores Dataframe of scores returned from \code{call_drug_hits}.
#' @param residuals Residuals returned with the return_residuals argument set to true
#'   from \code{call_drug_hits}.
#' @param control_name Name of control passed to \code{call_drug_hits}.
#' @param condition_name Name of condition passed to \code{call_drug_hits}.
#' @param output_folder Folder to output plots to. 
#' @param neg_type Label for significant effects with a negative differential effect
#'   passed to \code{call_drug_hits} (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   passed to \code{call_drug_hits} (default "Positive").
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export 
plot_drug_residuals <- function(scores, residuals, control_name, 
                                condition_name, output_folder,
                                neg_type = "Negative", pos_type = "Positive", 
                                plot_type = "png") {
  
  # Checks plot type and converts to lowercase
  plot_type <- tolower(plot_type)
  if (plot_type != "png" & plot_type != "pdf") {
    stop("plot_type must be either png or pdf")
  }
  
  # Makes output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  # Gets top hits
  response_col <- paste0("effect_type_", condition_name)
  control_col <- paste0("mean_", control_name)
  condition_col <- paste0("mean_", condition_name)
  diff_col <- paste0("differential_", condition_name, "_vs_", control_name)
  scores <- scores[scores[[response_col]] != "None",]
  residuals <- residuals[residuals$n %in% as.numeric(rownames(scores)),]
  residuals$lfc <- residuals[[condition_col]] - residuals[[control_col]]
  
  # Returns if scores have no hits
  if (nrow(scores) == 0) {
    cat(paste("No significant hits for", control_name, "vs", condition_name, "guide-level plotting\n"))
  } else {
    
    # Gets ranking of top hits
    neg_order <- order(scores[[diff_col]])
    scores$neg_rank <- NA
    scores$pos_rank <- NA
    scores$neg_rank[neg_order] <- 1:nrow(scores)
    scores$pos_rank[neg_order] <- nrow(scores):1
    
    # Makes LFC plots for all top hits
    for (i in unique(residuals$n)) {
      
      # Gets data and gene names
      df <- residuals[residuals$n == i,]
      ind <- which(as.numeric(rownames(scores)) == i)
      gene <- scores$gene[ind]
      x_label <- paste0("Guides")
      y_label <- paste0("Average differential LFC across replicates")
      
      # Adds ID column for plotting
      df$ID <- paste("Guide", 1:nrow(df))
      
      # Plots data
      p <- ggplot2::ggplot(df) +
        ggplot2::xlab(x_label) +
        ggplot2::ylab(y_label) +
        ggplot2::geom_bar(ggplot2::aes_string(x = "ID", y = "lfc"), stat = "identity", color = "Black", 
                          fill = ggplot2::alpha(c("gray30"), 1)) +
        ggplot2::geom_hline(yintercept = 1, linetype = 2, size = 1, alpha = 0.75, color = "Yellow") +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 0.75, color = "Gray") +
        ggplot2::geom_hline(yintercept = -1, linetype = 2, size = 1, alpha = 0.75, color = "Blue") +
        ggplot2::coord_flip() +
        ggthemes::theme_tufte(base_size = 20)
      
      # Gets type and rank of effect
      effect <- ""
      rank <- 0
      effect_type <- scores[[response_col]][ind]
      if (effect_type == neg_type) {
        effect <- "neg"
        rank <- scores$neg_rank[ind]
      } else {
        effect <- "pos"
        rank <- scores$pos_rank[ind]
      }
      
      # Saves to file
      file_name <- paste0(effect, "_", rank, "_", gene, ".", plot_type)
      suppressWarnings(ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300))
    }
  }
}

#' Plot guide-level residuals for a given gene.
#' 
#' Plots guide-level residuals for a given gene and returns the plot. Works for data 
#' returned from \code{call_drug_hits}.
#' 
#' @param scores Dataframe of scores returned from \code{call_drug_hits}.
#' @param residuals Residuals returned with the return_residuals argument set to true
#'   from \code{call_drug_hits}.
#' @param gene Gene name for guides to plot.
#' @param control_name Name of control passed to \code{call_drug_hits}.
#' @param condition_name Name of condition passed to \code{call_drug_hits}.
#' @return A ggplot object.
#' @export 
plot_gene_residuals <- function(scores, residuals, gene, control_name, 
                                condition_name) {
  
  # Gets guide-level residuals for the given gene
  ind <- which(scores$gene == gene)
  df <- residuals[residuals$n == ind,]
  x_label <- paste0("Guides")
  y_label <- paste0("Average differential LFC across replicates")
  
  # Computes residuals
  control_col <- paste0("mean_", control_name)
  condition_col <- paste0("mean_", condition_name)
  df$lfc <- df[[condition_col]] - df[[control_col]]
  
  # Adds ID column for plotting
  df$ID <- paste("Guide", 1:nrow(df))
  
  # Plots data and returns plot
  p <- ggplot2::ggplot(df) +
    ggplot2::xlab(x_label) +
    ggplot2::ylab(y_label) +
    ggplot2::geom_bar(ggplot2::aes_string(x = "ID", y = "lfc"), stat = "identity", color = "Black", 
                      fill = ggplot2::alpha(c("gray30"), 1)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, size = 1, alpha = 0.75, color = "Yellow") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 0.75, color = "Gray") +
    ggplot2::geom_hline(yintercept = -1, linetype = 2, size = 1, alpha = 0.75, color = "Blue") +
    ggplot2::coord_flip() +
    ggthemes::theme_tufte(base_size = 20)
  return(p)
}

#' Plot guide-level LFCs for a given gene.
#' 
#' Plots guide-level LFCs for a given gene and returns the plot. Works for data 
#' returned from \code{call_drug_hits}.
#' 
#' @param scores Dataframe of scores returned from \code{call_drug_hits}.
#' @param residuals Residuals returned with the return_residuals argument set to true
#'   from \code{call_drug_hits}.
#' @param gene Gene name for guides to plot.
#' @param control_name Name of control passed to \code{call_drug_hits}.
#' @param condition_name Name of condition passed to \code{call_drug_hits}.
#' @return A ggplot object.
#' @export 
plot_gene_lfc <- function(scores, residuals, gene, control_name, 
                          condition_name) {
  
  # Gets guide-level residuals for the given gene
  ind <- which(scores$gene == gene)
  df <- residuals[residuals$n == ind,]
  x_label <- paste0("Guides")
  y_label <- paste0("Average LFC across replicates")
  control_col <- paste0("mean_", control_name)
  condition_col <- paste0("mean_", condition_name)
  
  # Adds ID column for plotting
  df$ID <- paste("Guide", 1:nrow(df))
  
  # Reshapes to long format for plotting
  df <- stats::reshape(df, direction = "long", varying = c(control_col, condition_col), 
                       idvar = "ID", sep = "", timevar = "screen", v.names = "screen_LFC",
                       times = c(control_name, condition_name))
  
  # Plots data
  p <- ggplot2::ggplot(df) +
    ggplot2::xlab(x_label) +
    ggplot2::ylab(y_label) +
    ggplot2::geom_bar(ggplot2::aes_string(x = "ID", y = "screen_LFC"), stat = "identity", color = "Black", 
                      fill = ggplot2::alpha(c("gray30"), 1)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, size = 1, alpha = 0.75, color = "Yellow") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 0.75, color = "Gray") +
    ggplot2::geom_hline(yintercept = -1, linetype = 2, size = 1, alpha = 0.75, color = "Blue") +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(. ~ screen) +
    ggthemes::theme_tufte(base_size = 20, base_family = "sans") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, color = "black"))
  return(p)
}

#' Makes scree plot
#' 
#' Outputs scree plot for a dataset to help select a number of principal components to
#' remove. 
#' 
#' @param df LFC dataframe.
#' @param cols Numerical column names to normalize with PCA. 
#' @param scale Whether or not to scale replicates before extracting principal 
#'   components (default FALSE).
#' @param na_behavior Whether to replace NAs with row (per-guide) mean values or to
#'   omit NAs, as either "mean_replace" or "omit" (default "mean_replace").
#' @param exclude_screens A list of screen names to exclude, e.g. for replicates with
#'   mostly NA values (default NULL).
#' @return Scree plot as a ggplot2 object.
#' @export 
plot_scree <- function(df, cols, scale = FALSE, na_behavior = "mean_replace",
                       exclude_screens = NULL) {
  
  # Replaces NAs with row means
  temp <- data.matrix(df[,cols])
  if (na_behavior == "mean_replace") {
    for (i in 1:nrow(temp)) {
      na_ind <- is.na(temp[i,])
      if (any(na_ind)) {
        if (!all(na_ind)) {
          temp[i,na_ind] <- mean(temp[i,], na.rm = TRUE)
        } else {
          temp[i,] <- mean(temp, na.rm = TRUE)
        }
      }
    }
  } else if (na_behavior == "omit") {
    temp <- stats::na.omit(temp)
  } else {
    stop("na_behavior must be either 'mean_replace' or 'omit'")
  }
  
  # Makes scree plot
  n_components <- length(cols) - 1
  pca <- stats::prcomp(temp, center = TRUE, scale. = scale)
  var_df <- data.frame(PC = paste0("PC", 1:length(pca$sdev)),
                       var = 100 * (pca$sdev ^ 2) / sum(pca$sdev ^ 2))
  var_df$PC <- factor(var_df$PC, levels = var_df$PC)
  p <- ggplot2::ggplot(var_df, ggplot2::aes_string(x = "PC", y = "var"))+
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("Principal component") +
    ggplot2::ylab("% variance explained") +
    ggthemes::theme_tufte(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))
  
  # Returns plot
  return(p)
}
