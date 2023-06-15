######
# PRE-PROCESSING CODE
######

#' Normalizes reads for the given screens
#' 
#' Log2 and depth-normalizes reads between a given list of columns and a given
#' column of an earlier timepoint (e.g. T0) specified in the "normalize_name"
#' entry of each screen. Screens with NULL for their "normalize_name" entry
#' are log2 and depth-normalized, but not normalized to earlier timepoints.
#' If a screen to normalize against has multiple replicates, those replicates
#' are averaged before normalization. Multiple replicates for a screen being
#' normalized, however, are normalized separately against the provided early
#' timepoint.
#' 
#' @param df Reads dataframe.
#' @param screens List of screens generated with \code{add_screens}. 
#' @param filter_names List of screen names to filter based on read counts by. 
#' @param cf1 Scaling factor (default 1e6).
#' @param cf2 Pseudocount (default 1).
#' @param min_reads Minimum number of reads to keep (default 30, anything
#'   below this value will be filtered out).
#' @param max_reads Maximum number of reads to keep (default 10000, anything
#'   above this value will be filtered out).
#' @param nonessential_norm Whether or not to normalize each screen against its
#'   population of core non-essential genes, as defined by Traver et al. 2015 
#'   (default FALSE).
#' @param replace_NA Whether or not to replace NA and NULL values in non-T0 screens 
#'   with 0's after filtering out T0 guides with too few or NA readcounts 
#'   (default TRUE).
#' @return Normalized dataframe.
#' @export 
normalize_screens <- function(df, screens, filter_names = NULL, cf1 = 1e6, cf2 = 1, 
                              min_reads = 30, max_reads = 10000, nonessential_norm = TRUE,
                              replace_NA = TRUE) {
  
  # Checks for input errors
  check_screen_params(df, screens)
  
  # Loads gene standard from internal data
  nonessentials <- traver_nonessentials
  
  # Ensures that screens are only filtered by valid screen names
  if (length(filter_names) > 0) {
    all_names <- names(screens)
    names_to_keep <- rep(TRUE, length(filter_names))
    for (i in 1:length(filter_names)) {
      name <- filter_names[i]
      if (!(name %in% all_names)) {
        cat(paste("WARNING: screen", name, "not found, data not filtered by this screen\n"))
        names_to_keep[i] <- FALSE
      }
    }
    filter_names <- filter_names[names_to_keep] 
  }
  
  # Flags guides with too few or too many read counts if there are valid names to filter by
  to_remove <- rep(FALSE, nrow(df))
  sum_low <- 0
  sum_high <- 0
  sum_na <- 0
  if (length(filter_names) > 0) {
    filter_cols <- sapply(screens[filter_names], "[[", "replicates")
    for (col in filter_cols) {
      to_remove[df[,col] < min_reads] <- TRUE
    }
    sum_low <- sum(to_remove, na.rm = TRUE)
    for (col in filter_cols) {
      to_remove[df[,col] > max_reads] <- TRUE
    }
    sum_high <- sum(to_remove, na.rm = TRUE) - sum_low
    for (col in filter_cols) {
      to_remove[is.na(df[,col])] <- TRUE
    }
    sum_na <- sum(to_remove, na.rm = TRUE) - (sum_high + sum_low)
    removed_guides_ind <- which(to_remove) 
  }
  
  # Log2 and depth-normalizes every screen
  for (screen in screens) {
    for (col in screen[["replicates"]]) {
      df[,col] <- normalize_reads(df[,col], cf1, cf2)
    }
  }
  
  # Replaces NAs in non-T0 screens with 0 values if specified. Also affects NULL values
  if (replace_NA) {
    for (screen in screens) {
      for (col in screen[["replicates"]]) {
        na_ind <- !stats::complete.cases(df[,col])
        df[na_ind, col] <- 0
      }
    }
  }
  
  # Normalizes specified screens to earlier timepoints
  new_df <- df
  for (screen in screens) {
    normalize_name <- screen[["normalize_name"]]
    if (!is.null(normalize_name)) {
      if (normalize_name %in% all_names) {
        for (col in screen[["replicates"]]) {
          rep_cols <- screens[[normalize_name]][["replicates"]]
          rep_norm <- df[,rep_cols]
          if (length(rep_cols) > 1) {
            rep_norm <- rowMeans(rep_norm, na.rm = TRUE)
          }
          new_df[,col] <- df[,col] - rep_norm
        }
      } else {
        cat(paste("WARNING: screen", normalize_name, "not found.\n"))
      }
    }
  }
  
  # Normalizes all LFCs against per-screen non-essential median LFCs if specified
  nonessential_ind <- new_df$gene %in% nonessentials
  if (sum(nonessential_ind) > 0 & nonessential_norm) {
    for (screen in screens) {
      for (col in screen[["replicates"]]) {
        nonessential_vals <- as.numeric(unlist(new_df[nonessential_ind, col]))
        nonessential_median <- stats::median(nonessential_vals, na.rm = TRUE)
        new_df[,col] <- new_df[,col] - nonessential_median
      }
    }
  } else {
    cat(paste("Skipping LFC normalization to median per-screen non-essential gene LFC\n"))
  }
  
  # Removes flagged guides if applicable
  if (length(filter_names) > 0) {
    new_df <- new_df[!to_remove,]
    cat(paste("Excluded a total of", sum_low, "guides for low t0 representation\n"))
    cat(paste("Excluded a total of", sum_high, "guides for high t0 representation\n"))
    cat(paste("Excluded a total of", sum_na, "guides with NA t0 values\n"))
  } else {
    cat(paste("Filtering skipped because no valid screens were specified\n"))
  }
  
  # Explicitly returns filtered, processed data
  return(new_df)
}

#' PCA-normalizes residual LFCs for the given screens
#' 
#' PCA-normalizes residual LFCs for a dataset by extracting a given number of principal
#' components from the dataset, projecting the dataset to those components, and 
#' subtracting the projected dataset from the original dataset. After performing 
#' PCA-normalization, consider re-centering LFCs by the mean of non-essential genes. 
#' 
#' @param df LFC dataframe.
#' @param cols Numerical column names to normalize with PCA. 
#' @param n_components Vector containing indices of principal components to remove from data. 
#' @param scale Whether or not to scale replicates before extracting principal 
#'   components (default FALSE).
#' @param na_behavior Whether to replace NAs with row (per-guide) mean values or to
#'   omit NAs, as either "mean_replace" or "omit" (default "mean_replace").
#' @param exclude_screens A list of screen names to exclude, e.g. for replicates with
#'   mostly NA values (default NULL).
#' @return PCA-normalized dataframe.
#' @export 
pca_screens <- function(df, cols, n_components, scale = FALSE, na_behavior = "mean_replace",
                        exclude_screens = NULL) {
  
  # Checks number of components to remove
  valid_ind <- 1:length(cols)
  if (length(n_components) > length(cols)) {
    stop("ERROR: Specified too many PCs relative to number of columns in data\n")
  } else if (!(all(n_components %in% valid_ind))) {
    stop("ERROR: Specified invalid PCs to remove from data \n")
  }
  
  # Replaces NAs with row means
  na_mask <- is.na(df[,cols])
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
  
  # PCA-normalizes data
  pca <- stats::prcomp(temp, center = TRUE, scale. = scale)
  real_v <- pca$rotation[,n_components]
  projected <- temp %*% real_v %*% t(real_v)
  df[,cols] <- temp - projected
  df[,cols][na_mask] <- NA
  
  # Returns PCA-normalized data
  return(df)
}

#' Filters guides with too few read counts.
#' 
#' Filters guides out with too few read counts from a given reads dataframe
#' and a given set of columns (e.g. T0 columns).
#' 
#' @param df Reads dataframe.
#' @param cols Columns to filter by.
#' @param min_reads Minimum number of reads to keep (anything below
#'   this value will be filtered out).
#' @param max_reads Maximum number of reads to keep (anything above
#'   this value will be filtered out).
#' @return Filtered dataframe.
filter_reads <- function(df, cols, min_reads = 30, max_reads = 10000) {
  to_remove <- rep(FALSE, nrow(df))
  for (col in cols) {
    to_remove[df[,col] < min_reads] <- TRUE
  }
  sum_low <- sum(to_remove)
  for (col in cols) {
    to_remove[df[,col] > max_reads] <- TRUE
  }
  sum_high <- sum(to_remove) - sum_low
  removed_guides_ind <- which(to_remove)
  df <- df[!to_remove,]
  cat(paste("Excluded a total of", sum_low, "guides for low t0 representation\n"))
  cat(paste("Excluded a total of", sum_high, "guides for low t0 representation\n"))
  return(df)
}

#' Computes essential gene recovery AUC.
#' 
#' Computes area under the curve for ROC curves that measure how well each technical replicate
#' recovers signal for essential-targeting guides. Only computes AUC for guides that target 
#' essential genes twice, guides that target two different essential genes, or guides that 
#' target an essential gene and an intergenic region.
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}. 
#' @param negative_controls List of negative control genes to append to default list of
#'   non-essential genes (default NULL).
#' @param append_to_negatives Whether to append the negative controls specified in
#'   negative_controls to the list of gold-standard human non-essential genes if TRUE or to
#'   replace the list of non-essential genes if FALSE (default FALSE).
#' @return Returns a dataframe with three columns for replicate name, essential AUC relative 
#'   to all other genes, and essential AUC relative to a specified set of non-essentials.
essential_lfc_qc <- function(df, screens, negative_controls = NULL, 
                             append_to_negatives = FALSE) {
  
  # Checks that the given gene_col is in the data
  if (!("gene" %in% colnames(df))) {
    stop(paste("gene name column 'gene' not in df"))
  }
  
  # Loads gene standards from internal data
  essentials <- traver_core_essentials
  nonessentials <- traver_nonessentials
  
  # Adds genes to nonessentials list if specified
  if (!is.null(negative_controls)) {
    if (append_to_negatives) {
      nonessentials <- c(nonessentials, negative_controls) 
    } else {
      nonessentials <- negative_controls
    }
  }
  
  # Gets indices of essential and non-essential guides in data
  essential_ind <- df$gene %in% essentials
  nonessential_ind <- df$gene %in% nonessentials
  
  # Throws warning if too few genes in standards
  skip_nonessential <- FALSE
  if (sum(essential_ind) < 10) {
    warning(paste("too few essential-targeting guides in df, skipping all AUC computation"))
    return(NULL)
  }
  if (sum(nonessential_ind) < 10) {
    warning(paste("too few nonessential-targeting guides in df, skipping non-essential AUC computation"))
    skip_nonessential <- TRUE
  }
  
  # Gets PR curves for all essential genes and all technical replicates
  results <- data.frame(screen = NA,
                        replicate = NA, 
                        essential_AUC_all = NA,
                        essential_AUC_nonessential = NA,
                        essential_AUC_all_gene = NA,
                        essential_AUC_nonessential_gene = NA)
  counter <- 1
  for (screen_name in names(screens)) {
    screen <- screens[[screen_name]]
    for (rep in screen[["replicates"]]) {
      temp <- df[!is.na(df[[rep]]),]
      essential_ind_rep <- temp$gene %in% essentials
      nonessential_ind_rep <- NULL
      if (!skip_nonessential) {
        nonessential_ind_rep <- temp$gene %in% nonessentials
      }
      
      # Computes AUC relative to all guides
      auc1 <- NA
      auc2 <- NA
      if (sum(essential_ind_rep) < 10) {
        warning(paste("too few essential-targeting guides for replicate", rep, ", skipping AUC computation"))
      } else {
        roc <- PRROC::roc.curve(-temp[[rep]], weights.class0 = as.numeric(essential_ind_rep), curve = TRUE)
        auc1 <- roc$auc 
      }
      
      # Computes AUC relative to nonessential guides
      if (!skip_nonessential & sum(nonessential_ind_rep) > 10) {
        essential_rownames <- rownames(temp)[essential_ind_rep]
        temp <- temp[essential_ind_rep | nonessential_ind_rep,]
        temp_essential_ind <- rownames(temp) %in% essential_rownames
        roc <- PRROC::roc.curve(-temp[[rep]], weights.class0 = as.numeric(temp_essential_ind), curve = TRUE)
        auc2 <- roc$auc
      }
      
      # Merges guide-level data to gene-level data
      gene_df <- data.frame(gene = unique(temp$gene), 
                            val = NA)
      for (i in 1:nrow(gene_df)) {
        gene <- gene_df$gene[i]
        gene_df$val[i] <- mean(temp[temp$gene == gene, rep], na.rm = TRUE)
      }
      temp <- gene_df
      
      # Gets gene-level indices
      essential_ind_gene <- temp$gene %in% essentials
      nonessential_ind_gene <- NULL
      if (!skip_nonessential) {
        nonessential_ind_gene <- temp$gene %in% nonessentials
      }
      
      # Computes AUC relative to all genes
      auc3 <- NA
      auc4 <- NA
      if (sum(essential_ind_gene) < 10) {
        warning(paste("too few essential-targeting guides for replicate", rep, ", skipping AUC computation"))
      } else {
        roc <- PRROC::roc.curve(-temp$val, weights.class0 = as.numeric(essential_ind_gene), curve = TRUE)
        auc3 <- roc$auc 
      }
      
      # Computes AUC relative to nonessential genes
      if (!skip_nonessential & sum(nonessential_ind_gene) > 10) {
        essential_rownames <- rownames(temp)[essential_ind_gene]
        temp <- temp[essential_ind_gene | nonessential_ind_gene,]
        temp_essential_ind <- rownames(temp) %in% essential_rownames
        roc <- PRROC::roc.curve(-temp$val, weights.class0 = as.numeric(temp_essential_ind), curve = TRUE)
        auc4 <- roc$auc
      }
      results[counter,] <- c(screen_name, rep, auc1, auc2, auc3, auc4)
      counter <- counter + 1
    }
  }
  return(results)
}

#' Log-normalizes reads.
#' 
#' Log2-normalizes reads with a given pseudocount and scaling factor, and also
#' depth-normalizes the data. 
#' 
#' @param df List of read counts.
#' @param cf1 Scaling factor (default 1e6).
#' @param cf2 Pseudocount (default 1).
#' @return Log- and depth-normalized read counts.
#' @export
normalize_reads <- function(df, cf1 = 1e6, cf2 = 1) {
  log2((df / sum(df, na.rm = TRUE)) * cf1 + cf2)
}

