######
# SCORING CODE
######

#' @importFrom magrittr "%>%"

# Inner function to scale values between 0 and 1
scale_values <- function(x) {
  val <- (x-min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
}

#' Scores conditions against a single control.
#' Revised one-off score function with individual replicate-level loess fitting (loess smoothing Drug A vs DMSO A, B vs B, C vs C per guide)
#' 
#' Scores guides for any number of drug screens against a control screen
#' (e.g. for directly comparing drug response to DMSO response). After running 
#' this function, pass the resulting dataframe to \code{call_drug_hits} to 
#' call significant effects.
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}.
#' @param control_screen_name Name of a control screen to test condition screens against.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param control_genes List of control genes to remove, e.g. "luciferase" (default c("None", "")).
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ma_transform If true, M-A transforms data before running loess normalization. Only
#'   has an effect when loess = TRUE (default TRUE).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default "BY").
#' @param sd_scale_factor Factor to normalize SDs against for scaling. If NULL, this operation
#'   is not performed - this behavior is different for group scoring! (default NULL).
#' @param return_residuals If FALSE, returns NA instead of residuals dataframe (default TRUE).
#'   This is recommend if scoring large datasets and memory is a limitation.  
#' @param verbose If true, prints verbose output (default FALSE). 
#' @return A list containing two dataframes. The first entry, named "scored_data" in the list,
#'   contains scored data with separate columns given by the specified control and condition
#'   names. The second entry, named "residuals" in the list, is a dataframe containing control,
#'   condition and loess-normalized residuals for all guides.
#' @export
score_drugs_vs_control <- function(df, screens, control_screen_name, condition_screen_names, 
                                   control_genes = c("None", ""), min_guides = 3, test = "moderated-t", 
                                   loess = TRUE, ma_transform = TRUE, fdr_method = "BY",
                                   sd_scale_factor = NULL, return_residuals = TRUE, verbose = FALSE) {
  
  # Gets condition names and columns for any number of conditions
  if (verbose) {
    cat(paste("Preparing to score...\n"))
  }
  control_name <- control_screen_name
  control_cols <- screens[[control_name]][["replicates"]]
  condition_names <- c()
  condition_cols <- list()
  for (condition in condition_screen_names) {
    condition_names <- c(condition_names, condition)
    condition_cols[[condition]] <- screens[[condition]][["replicates"]]
  }
  
  # Removes control genes
  df <- df[!(df$gene %in% control_genes),]
  
  # Makes output dataframe
  unique_genes <- unique(df$gene)
  n_genes <- length(unique_genes)
  scores <- data.frame(gene = rep(NA, n_genes))
  
  #^^^^^^^^^^^^^^^^^^^^^^^^^changed^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  # Makes condition residual dataframes if necessary  ^^^^^^^^^^^^^^^^^^^ 
  max_guides <- -1
  condition_residuals <- list()
  if (test == "moderated-t") {
    # Gets max number of guides first
    max_guides <- max(table(df$gene))
    # Makes residual dataframes with columns equal to the max number of guides
    for (name in condition_names) {
      # Extract replicate names (add to col name) and count from condition_cols for condition name
      condition_reps=condition_cols[[name]]
      residual_df <- data.frame(matrix(nrow = n_genes, ncol = max_guides*length(condition_reps)))# ncol = max_guides*total_replicate
      #colnames(residual_df) <- paste0("guide_residual_", 1:max_guides,'_', rep(condition_reps,each=max_guides)) # version1
      colnames(residual_df) <- paste0("guide_residual_", rep(1:max_guides,each=length(condition_reps)),'_', condition_reps) # version2
      condition_residuals[[name]] <- residual_df
    }
  }
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^changed^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  # Pairwise LOESS smoothing
  # Computes loess-normalized residuals if specified
  loess_residuals <- NULL
  if (loess & test == "moderated-t") {
    loess_residuals <- data.frame(gene = df$gene) # get gene names from main dataframe
    for (name in condition_names) {
      condition_reps=condition_cols[[name]]
      for(rep_index in c(1:length(condition_reps))){ # pair-wise control-condition rep extracted with index
        control_values <- df[,paste0(control_cols[rep_index])]  #get control replicate scores from main dataframe
        rep = paste0(condition_reps[rep_index])
        condition_values <- df[,rep] # get condition replicate scores from main dataframe
        ##################loess qGI_utils.R#########################
        temp <- loess_MA(control_values, condition_values, ma_transform = ma_transform)
        #####################################################
        loess_residuals[[paste0("loess_residual_", name,'_',rep)]] <- temp[["residual"]]
        loess_residuals[[paste0("loess_predicted_", name,'_',rep)]] <- temp[["predicted"]]
      }
    }
  }
  
  # Appends additional columns for each condition
  new_cols <- c(paste0("n_", control_name), 
                paste0("mean_", control_name),
                paste0("variance_", control_name))
  for (name in condition_names) {
    new_cols <- c(new_cols, c(
      paste0("n_", name), 
      paste0("mean_", name),
      paste0("variance_", name),
      paste0("differential_", name, "_vs_", control_name),
      paste0("pval_", name, "_vs_", control_name),
      paste0("fdr_", name, "_vs_", control_name),
      paste0("significant_", name, "_vs_", control_name)
    ))
  }
  scores[new_cols] <- NA
  
  # Scores guides for each condition
  counter <- 1
  for (i in 1:n_genes) {
    # Gets gene names and control guide values across replicates and removes 
    # NaNs introduced by missing guides
    gene <- unique_genes[i]
    guide_vals <- df[df$gene == gene,]
    scores$gene[i] <- gene
    rep_mean_control <- rowMeans(data.frame(guide_vals[control_cols]), na.rm = TRUE)
    keep_ind <- !is.nan(rep_mean_control)
    rep_mean_control <- rep_mean_control[keep_ind]
    
    # Skips if too few guides
    if (length(rep_mean_control) < min_guides) {
      next
    }
    
    # Takes the mean across replicates for all conditions
    for (name in condition_names) {
      # Gets residual LFCs across replicates after removing NaNs
      rep_mean_condition <- rowMeans(data.frame(guide_vals[condition_cols[[name]]]), na.rm = TRUE)
      rep_mean_condition <- rep_mean_condition[keep_ind]
      rep_mean_condition[is.nan(rep_mean_condition)] <- NA
      diff <- rep_mean_condition - rep_mean_control
      
      # Stores gene-level stats
      scores[[paste0("n_", control_name)]][i] <- length(rep_mean_control)
      scores[[paste0("n_", name)]][i] <- length(rep_mean_condition)
      scores[[paste0("mean_", control_name)]][i] <- mean(rep_mean_control, na.rm = TRUE)
      scores[[paste0("mean_", name)]][i] <- mean(rep_mean_condition, na.rm = TRUE)
      scores[[paste0("variance_", control_name)]][i] <- stats::var(rep_mean_control, na.rm = TRUE)
      scores[[paste0("variance_", name)]][i] <- stats::var(rep_mean_condition, na.rm = TRUE)
      scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(diff, na.rm = TRUE)      
      
      #^^^^^^^^^^^^^^^^^^^^^^^^^changed^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      # Performs the specified type of testing or stores residuals for later testing
      if (test == "rank-sum") {
        scores[[paste0("pval_", name, "_vs_", control_name)]][i] <- 
          suppressWarnings(stats::wilcox.test(rep_mean_condition, rep_mean_control))$p.value
      } else if (test == "moderated-t") {
        loess_residual_rep = c()
        condition_reps = condition_cols[[name]]
        for(rep_index in c(1:length(condition_reps))){ 
          if (loess) {
            resid <- loess_residuals[[paste0("loess_residual_", name, '_', condition_reps[rep_index])]][loess_residuals$gene == gene]
            predicted <- loess_residuals[[paste0("loess_predicted_", name, '_', condition_reps[rep_index])]][loess_residuals$gene == gene]
            if (length(resid) < max_guides) { 
              resid <- c(resid, rep(NA, max_guides - length(resid))) 
            }
            condition_residuals[[name]][i, paste0("guide_residual_", 1:max_guides,'_', rep(condition_reps[rep_index], max_guides))] <- resid
            loess_residual_rep <- c(loess_residual_rep, mean(resid, na.rm = TRUE))
          }else{
            resid <- guide_vals[condition_reps[rep_index]] - guide_vals[control_cols[rep_index]]
            resid=resid[[condition_reps[rep_index]]] 
            if (length(resid) < max_guides) { resid <- c(resid, rep(NA, max_guides - length(resid))) } 
            condition_residuals[[name]][i, paste0("guide_residual_", 1:max_guides,'_', rep(condition_reps[rep_index], max_guides))] <- resid 
          }
        }
        if (loess) {
          # Update differential LFC values
          scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(loess_residual_rep, na.rm = TRUE)
        }
      }
    }
    counter <- counter + length(rep_mean_control)
  }
  
  # Scores condition response with moderated t-test
  if (test == "moderated-t") {
    for (name in condition_names) {
      #block <- rep(1:length(condition_cols[[name]]),each=max_guides) #version1 #group same replicates together -- within block guides
      #block <- rep(1:max_guides, length(condition_cols[[name]])) #version1_flipped #group same guides together -- within block replicates
      block <- rep(1:max_guides, each=length(condition_cols[[name]])) #version2 #group same guides together -- within block replicates
      dupcor <- limma::duplicateCorrelation(condition_residuals[[name]], block=block)
      ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]], design=NULL, block=block, correlation=dupcor$consensus))
      p_val <- ebayes_fit$p.value[,1]
      scores[[paste0("pval_", name, "_vs_", control_name)]] <- p_val
    }   
  }
  
  # Scales moderate effects in top and bottom 10% of data to de-emphasize those. 
  # The mean to divide SD values by is a pre-computed scalar
  if (!is.null(sd_scale_factor)) {
    for (name in condition_names) {
      resid <- condition_residuals[[name]]
      lfc_range <- stats::quantile(resid, probs = c(0.1, 0.9), na.rm = TRUE)
      target_sd <- stats::sd(resid[resid > lfc_range[1] & resid < lfc_range[2]], na.rm = TRUE)
      target_sd <- target_sd / sd_scale_factor
      condition_residuals[[name]] <- resid / target_sd
      mean_residuals <- rowMeans(condition_residuals[[name]])
      scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(diff, na.rm = TRUE)
    } 
  }
  
  # Computes FDRs
  for (name in condition_names) {
    scores[[paste0("fdr_", name, "_vs_", control_name)]] <- 
      stats::p.adjust(scores[[paste0("pval_", name, "_vs_", control_name)]], method = fdr_method)
  }
  
  # Removes extra zero row from loess-normalized residuals - no need for the new method
  # if (loess) {
  #   loess_residuals <- loess_residuals[1:(nrow(loess_residuals) - 1),]
  # }
  
  # Explicitly returns scored data
  output <- list()
  output[["scored_data"]] <- scores
  if (return_residuals & loess) {
    output[["residuals"]] <- loess_residuals
  } else if (return_residuals) {
    cat("WARNING: returning residuals is currently only supported with loess-normalization enabled\n")
    output[["residuals"]] <- NA
  } else {
    output[["residuals"]] <- NA
  }
  return(output)
}

#' Scores conditions against multiple controls.
#' 
#' Scores guides for any number of drug screens against multiple control screens
#' (e.g. for directly comparing drug response to DMSO response). After running 
#' this function, pass the resulting dataframe to \code{call_drug_hits} to 
#' call significant effects.
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}.
#' @param control_screen_names A list of control screen names to test the condition
#'   screens against.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param matched_controls A list of control screen names with each screen corresponding
#'   to the condition at the same position in condition_screen_names. 
#' @param output_folder Folder to output scree plot to.
#' @param control_genes List of control genes to remove, e.g. "luciferase" (default c("None", "")).
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ma_transform If true, M-A transforms data before running loess normalization. Only
#'   has an effect when loess = TRUE (default TRUE).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default "BY").
#' @param weight_method If "linear" weights non-matched control screens equally. If exponential,
#'   instead weights non-matched control screens according to a decreasing exponential function
#'   based on the similarity of control screens to the given condition screen.
#' @param matched_fraction The weight given to the matched control as a fraction of 1, where
#'   all non-matched controls receive a total weight equal to 1 - matched_fraction (default 0.75).
#' @param sd_scale_factor Factor to normalize SDs against for scaling. If NULL, the mean is computed 
#'   across guide-level residuals, otherwise the given scalar is used instead (default NULL).
#' @param n_components Vector containing indices of principal components to remove from data. 
#' @param chromosomal_correction If TRUE, corrects chromosomal shifts by down-weighting qGI
#'   scores for shifted regions (default FALSE).
#' @param return_residuals If FALSE, returns NA instead of residuals dataframe (default TRUE).
#'   This is recommend if scoring large datasets and memory is a limitation.  
#' @param intermediate_file Path to file to save intermediate output to, or NULL to not write
#'   any intermediate output (default NULL).
#' @param load_intermediate If TRUE, loads data from intermediate file instead of running 
#'   earlier steps in the scoring pipeline (default FALSE).
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @param verbose If true, prints verbose output (default FALSE). 
#' @return A list containing two dataframes. The first entry, named "scored_data" in the list,
#'   contains scored data with separate columns given by the specified control and condition
#'   names. The second entry, named "residuals" in the list, is a dataframe containing control,
#'   condition and loess-normalized residuals for all guides.
#' @export
score_drugs_vs_controls <- function(df, screens, control_screen_names, condition_screen_names, 
                                    matched_controls, output_folder, control_genes = c("None", ""), 
                                    min_guides = 2, loess = TRUE, ma_transform = TRUE, fdr_method = "BY", 
                                    weight_method = "exp", matched_fraction = 0.75, sd_scale_factor = NULL,
                                    n_components = 0, chromosomal_correction = FALSE, return_residuals = TRUE, 
                                    intermediate_file = NULL, load_intermediate = FALSE, plot_type = "png", 
                                    verbose = FALSE) {
  
  # Disables dplyr warnings
  options(dplyr.summarise.inform = FALSE)
  
  # Gets condition names and columns for all specified conditions and controls
  if (verbose) {
    cat(paste("Preparing to score...\n"))
  }
  control_names <- c()
  control_cols <- list()
  condition_names <- c()
  condition_cols <- list()
  all_condition_cols <- c()
  for (control in control_screen_names) {
    control_names <- c(control_names, control)
    control_cols[[control]] <- screens[[control]][["replicates"]]
  }
  for (condition in condition_screen_names) {
    condition_names <- c(condition_names, condition)
    condition_cols[[condition]] <- screens[[condition]][["replicates"]]
    all_condition_cols <- c(all_condition_cols, screens[[condition]][["replicates"]])
  }
  
  # Removes control genes
  df <- df[!(df$gene %in% control_genes),]
  
  # Sorts dataframe by alphabetical order of genes
  unique_genes <- sort(unique(df$gene))
  df <- df[order(match(df$gene, unique_genes)),]
  
  # Makes output dataframe
  n_genes <- length(unique_genes)
  n_control <- length(control_names)
  max_guides <- max(table(df$gene))
  scores <- data.frame(gene = unique_genes)
  
  # Makes residual dataframes with columns equal to the max number of guides
  condition_residuals <- list()
  for (name in condition_names) {
    residual_df <- data.frame(matrix(nrow = n_genes, ncol = max_guides*n_control))
    for (i in 1:max_guides)  {
      label <- paste0("guide_", i, "_residual_", 1:n_control)
      col_start <- 1 + ((i - 1) * n_control)
      col_end <- i * n_control
      colnames(residual_df)[col_start:col_end] <- label
    }
    condition_residuals[[name]] <- residual_df
  }
  
  # Takes replicate averages for all controls and conditions
  control_df <- data.frame(gene = df$gene)
  condition_df <- data.frame(gene = df$gene)
  for (control in names(control_cols)) {
    reps <- df[,control_cols[[control]]]
    control_df[[control]] <- rowMeans(reps, na.rm = TRUE)
  }
  for (condition in names(condition_cols)) {
    reps <- df[,condition_cols[[condition]]]
    condition_df[[condition]] <- rowMeans(reps, na.rm = TRUE)
  }
  
  # Ensures that NaNs are converted to NAs
  control_df[is.na(control_df)] <- NA
  condition_df[is.na(condition_df)] <- NA
  
  # Takes residuals between each control and each condition
  if (verbose & !load_intermediate) {
    cat(paste("Taking residuals between all condition-control LFC pairs...\n"))
  } else if (verbose & load_intermediate) {
    cat(paste("Skipping residual computation to load saved residuals...\n"))
  }
  residual_df <- data.frame(gene = control_df$gene)
  if (!load_intermediate) {
    for (condition in names(condition_cols)) {
      for (control in names(control_cols)) {
        col <- paste0(condition, "_vs_", control)
        if (!loess) {
          residual_df[[col]] <- condition_df[[condition]] - control_df[[control]]
        } else {
          
          # Performs loess-normalization if specified
          temp <- loess_MA(control_df[[control]], condition_df[[condition]], ma_transform = ma_transform)
          residual_df[[col]] <- temp[["residual"]]
        }
      }
    } 
  }
  
  # Weights matched controls to contribute 75% of the effect size, and all other controls
  # to contribute the other 25% of the effect size
  if (verbose & !load_intermediate) {
    cat(paste("Weighting controls for computing qGI scores...\n"))
  } else if (verbose & load_intermediate) {
    cat(paste("Skipping weighting to load saved weights...\n"))
  }
  
  # # Optionally loads weights if intermediate files are not loaded
  # save(control_df, condition_df, control_names, condition_names, matched_controls, weight_method, 
  #      matched_fraction, file = intermediate_file)
  control_scores <- NULL
  if (!load_intermediate) {
    control_scores <- compute_control_effects(control_df, control_names, method = "linear",
                                              loess = loess, ma_transform = ma_transform)
    weights <- compute_control_weights(control_scores, condition_df, control_names, condition_names,
                                       matched_controls, method = weight_method,
                                       matched_fraction = matched_fraction)
  }
  
  # Stores basic stats and fills residual dataframes for moderated t-testing
  if (!load_intermediate) {
    for (i in 1:length(condition_names)) {
      
      # Gets number of guides for each condition
      condition <- condition_names[i]
      if (verbose) {
        cat(paste("Storing gene-level stats for", condition, "\n"))
      }
      n_guides <- as.integer(ncol(condition_residuals[[name]]) / n_control)
      
      # Gets weighted control LFCs - importantly, we prevent rows of all NAs from
      # being returned as 0 and instead return them as NA with another call to
      # rowSums. See akrun's SO answer for more information: 
      # https://stackoverflow.com/questions/38544325/rowsums-na-na-gives-0
      temp <- as.matrix(control_df[,2:ncol(control_df)])
      temp <- temp %*% diag(weights[i,])
      temp <- rowSums(temp, na.rm = TRUE) * NA^!rowSums(!is.na(temp))
      temp <- data.frame(gene = control_df$gene, lfc = temp)
      temp <- temp %>%
        dplyr::group_by(gene) %>%
        stats::na.omit() %>%
        dplyr::summarise(mean = mean(lfc, na.rm = TRUE),
                         var = stats::var(lfc, na.rm = TRUE),
                         n = dplyr::n())
      
      # Ensures that all genes are represented before appending to scores,
      # with scores that didn't make it through instead given NA values
      if (!all(scores$gene %in% temp$gene)) {
        to_append <- scores$gene[!(scores$gene %in% temp$gene)]
        for (gene in to_append) {
          temp[nrow(temp) + 1,] <- list(gene, NA, NA, 0)
        }
        temp <- temp[order(temp$gene),]
      }
      scores[[paste0("n_controls_", condition)]] <- temp$n
      scores[[paste0("mean_controls_", condition)]] <- temp$mean
      scores[[paste0("variance_controls_", condition)]] <- temp$var
      
      # Gets un-normalized condition LFCs
      temp <- data.frame(gene = condition_df$gene, lfc = condition_df[,condition])
      temp <- temp %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(mean = mean(lfc, na.rm = TRUE),
                         var = stats::var(lfc, na.rm = TRUE),
                         n = dplyr::n())
      scores[[paste0("n_", condition)]] <- temp$n
      scores[[paste0("mean_", condition)]] <- temp$mean
      scores[[paste0("variance_", condition)]] <- temp$var
      
      # Fills residual dataframe for all genes individually
      for (j in 1:length(unique_genes)) {
        
        # Gets indices of corresponding residual columns
        col_start <- 2 + ((i - 1) * n_control)
        col_end <- 1 + i * n_control
        
        # Fills residual dataframe
        gene <- unique_genes[j]
        gene_residual <- residual_df[residual_df$gene == gene, col_start:col_end]
        if (nrow(gene_residual) >= min_guides) {
          for (guide_ind in 1:nrow(gene_residual)) {
            residuals <- rep(NA, n_control)
            diff <- gene_residual[guide_ind,]
            residuals[1:length(diff)] <- diff
            residual_start <- 1 + ((guide_ind - 1) * n_control)
            residual_end <- guide_ind * n_control
            condition_residuals[[condition]][j, residual_start:residual_end] <- residuals
          }
        } 
      }
    }
  } else if (verbose) {
    cat(paste("Skipping storing gene-level stats to load saved stats...\n"))
  }
  
  # Scores each condition and scores corresponding p-values and FDRs
  block <- rep(1:max_guides, each = n_control)
  if (!load_intermediate) {
    for (condition in condition_names) {
      if (verbose) {
        cat(paste("Computing FDRs for", condition, "...\n"))
      }
      corfit <- limma::duplicateCorrelation(condition_residuals[[condition]], ndups = 1, block = block)
      consensus_cor <- corfit$consensus.correlation
      if (consensus_cor <= 0) {
        warning(paste("consensus correlation for", condition, "is negative"))
      }
      fit <- limma::lmFit(condition_residuals[[condition]], block = block, correlation = consensus_cor)
      ebayes_fit <- limma::eBayes(fit)
      condition_residuals[[condition]]$mean <- fit$Amean
      condition_residuals[[condition]]$pval <- ebayes_fit$p.value
      condition_residuals[[condition]]$fdr <- stats::p.adjust(ebayes_fit$p.value, method = fdr_method)
      scores[[paste0("pval_", condition, "_vs_controls")]] <- condition_residuals[[condition]]$pval
      scores[[paste0("fdr_", condition, "_vs_controls")]] <- condition_residuals[[condition]]$fdr
    }
  } else if (verbose) {
    cat(paste("Skipping computing FDRs to load saved FDRs...\n"))
  }
  
  # Loads intermediate files if specified 
  if (load_intermediate) {
    if (file.exists(intermediate_file)) {
      cat(paste("Loading intermediate output...\n"))
      load(intermediate_file)
    } else {
      cat(paste("ERROR: If load_intermediate is TRUE, intermediate_file must specify path to existing intermediate output\n"))
    }
  }
  
  # Writes intermediate output to file if specified
  if (!is.null(intermediate_file) & !load_intermediate) {
    if (verbose) {
      cat(paste0("Saving intermediate output to ", intermediate_file, "...\n"))
    }
    save(scores, residual_df, control_df, weights, condition_residuals, control_scores, file = intermediate_file)
  }
  
  ## The remaining sections of code compute qGI scores, which are more heavily processed than FDRs
  
  # Multiplies guide-level differential LFCs by weights
  # save(scores, control_df, condition_df, control_names, condition_names, matched_controls, weight_method, 
  #      residual_df, weights, matched_fraction, n_control, n_components, file = intermediate_file)
  for (i in 1:nrow(weights)) {
    condition <- rownames(weights)[i]
    for (j in 1:ncol(weights)) {
      control <- colnames(weights)[j]
      col <- paste0(condition, "_vs_", control)
      if (col %in% colnames(residual_df)) {
        residual_df[[col]] <- residual_df[[col]] * weights[condition, control]
      } else {
        if (verbose) {
          cat(paste("Comparison", col, "not valid\n"))
        }
      }
    }
  }
  
  # Sums weighted LFCs
  qGI <- data.frame(gene = control_df$gene)
  for (i in 1:length(condition_names)) {
    condition <- condition_names[i]
    col_start <- 2 + ((i - 1) * n_control)
    col_end <- 1 + i * n_control
    qGI[[condition]] <- rowSums(residual_df[,col_start:col_end], na.rm = TRUE)
  }
  
  # Scales moderate effects in top and bottom 10% of data to de-emphasize those and 
  # merges pre- and post-scaling SDs into a single dataframe to be returned. The mean
  # to divide SD values by is either the mean of all SDs or a pre-computed scalar
  pre_scaling_sd <- apply(qGI[,2:ncol(qGI)], 2, stats::sd)
  lfc_range <- apply(qGI[,2:ncol(qGI)], 2, stats::quantile, probs = c(0.1, 0.9), na.rm = TRUE)
  target_sd <- rep(NA, ncol(qGI))
  for (i in 2:ncol(qGI)) {
    target_sd[i] <- stats::sd(qGI[qGI[,i] > lfc_range[1,i-1] & qGI[,i] < lfc_range[2,i-1], i], na.rm = TRUE)
  }
  if (is.null(sd_scale_factor)) {
    sd_scale_factor <- mean(target_sd[2:length(target_sd)])
  }
  target_sd <- target_sd / sd_scale_factor
  for (i in 2:ncol(qGI)) {
    qGI[,i] <- qGI[,i] / target_sd[i]
  }
  post_scaling_sd <- apply(qGI[,2:ncol(qGI)], 2, stats::sd)
  sd_table <- data.frame(rbind(pre_scaling_sd, post_scaling_sd))
  sd_table$source[1] <- "pre_scaling"
  sd_table$source[2] <- "post_scaling"
  
  # Removes principal components from data if specified
  if (n_components > 0) {
    if (verbose) {
      cat(paste("Removing principal components...\n"))
    }
    numerical_cols <- colnames(qGI)[2:ncol(qGI)]
    p <- plot_scree(qGI, numerical_cols)
    plot_file <- file.path(output_folder, paste0("scree_plot.", plot_type))
    suppressWarnings(ggplot2::ggsave(plot_file, p, dpi = 300, width = 10, height = 7))
    qGI <- pca_screens(qGI, numerical_cols, n_components)
  }
  
  # Runs chromosomal correction if possible and enabled
  if (chromosomal_correction) {
    if (all(c("chr", "start_loc", "stop_loc")  %in% colnames(df))) {
      qGI <- correct_chromosomal_effects(df, qGI)
    }
  }
  
  # Summarizes guide-level scores to gene-level scores
  if (verbose) {
    cat(paste("Storing qGI scores...\n"))
  }
  for (i in 1:length(condition_names)) {
    condition <- condition_names[i]
    temp <- data.frame(gene = residual_df$gene, resid = qGI[[condition]])
    temp <- temp %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(val = mean(resid, na.rm = TRUE))
    scores[[paste0("differential_", condition, "_vs_controls")]] <- temp$val
  }
  
  # Removes genes with too few control observations
  ind <- which(grepl("n_controls", colnames(scores)))[1]
  scores <- scores[scores[,ind] >= min_guides,]
  
  # Explicitly returns scored data
  if (verbose) {
    cat(paste("Scoring finished!\n"))
  }
  output <- list()
  output[["scored_data"]] <- scores
  if (return_residuals) {
    output[["residuals"]] <- condition_residuals
  } else {
    output[["residuals"]] <- NA
  }
  output[["scored_controls"]] <- control_scores
  output[["sd_table"]] <- sd_table
  return(output)
}

#' Corrects for GIs shifted in specific chromosomes.
#' 
#' Run this within \code{score_drugs_vs_controls} to implement qGI chromosomal correction 
#' steps.
#' 
#' @param df LFC dataframe passed into \code{score_drugs_vs_controls} with chromosomal location 
#'   specified in a column named "chr", start coordinates specified in "start_loc" and stop 
#'   coordinates specified in "stop_loc." 
#' @param guide_df Dataframe of guide-level qGI scores.
#'   
#' @return Dataframe of chromosomal-corrected, guide-level qGI scores. 
correct_chromosomal_effects <- function(df, guide_df) {
  
  # Gets guide density for each relevant chromosome
  chrom <- names(sort(table(df$chr)))
  chr_features <- array(NA, dim = c(length(chrom), 4),
                        dimnames = list(chrom, c("length", "n_genes", "length_per", "roll_window")))
  for (i in 1:length(chrom)) {
    chr_features[i,"length"] <- max(df$stop_loc[df$chr == chrom[i]]) - min(df$start_loc[df$chr == chrom[i]])
    chr_features[i,"n_genes"] <- length(which(df$chr == chrom[i]))
  }
  chr_features[,3] <- chr_features[,2] / chr_features[,1]
  chr_features[,4] <- ceiling(log2(chr_features[,3] / stats::median(chr_features[,3]) + 2) * 150)
  
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

#' Call significant responses for scored chemogenomic data.
#' 
#' Run this to call significant responses for data returned from 
#' \code{score_drugs_vs_control} or \code{score_drugs_vs_controls}.
#' 
#' @param scores Dataframe returned from \code{score_drugs_vs_control} or 
#'   \code{score_drugs_vs_controls}.
#' @param control_screen_name Name of a control screen to test condition screens against, 
#'   or NULL to score data returned from \code{score_drugs_vs_controls}.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @return Dataframe of scored data with differential effects called as significant
#'   for the specified conditions. 
#' @export
call_drug_hits <- function(scores, control_screen_name = NULL, condition_screen_names = NULL,
                           fdr_threshold = 0.1, differential_threshold = 0,
                           neg_type = "Negative", pos_type = "Positive") {
  
  # Gets condition names and columns for any number of conditions
  control_name <- control_screen_name
  condition_names <- c()
  for (condition in condition_screen_names) {
    condition_names <- c(condition_names, condition)
  }
  
  # Checks type of scoring
  if (!(is.null(control_name))) {
    
    # Calls significant differences for each condition against the control
    for (name in condition_names) {
      scores[[paste0("significant_", name, "_vs_", control_name)]] <- 
        scores[[paste0("fdr_", name, "_vs_", control_name)]] < fdr_threshold
    }
    
    # Makes thresholded calls for significant negative and positive effects
    for (name in condition_names) {
      response_col <- paste0("effect_type_", name)
      scores[[response_col]] <- "None"
      diffs <- scores[[paste0("differential_", name, "_vs_", control_name)]]
      sig <- scores[[paste0("significant_", name, "_vs_", control_name)]]
      scores[[response_col]][sig & diffs < 0 & abs(diffs) > differential_threshold] <- neg_type
      scores[[response_col]][sig & diffs > 0 & abs(diffs) > differential_threshold] <- pos_type
    }
    
  } else {
    
    # Calls significant differences for each condition against the control
    for (name in condition_names) {
      scores[[paste0("significant_", name, "_vs_controls")]] <- 
        scores[[paste0("fdr_", name, "_vs_controls")]] < fdr_threshold
    }
    
    # Makes thresholded calls for significant negative and positive effects
    for (name in condition_names) {
      response_col <- paste0("effect_type_", name)
      scores[[response_col]] <- "None"
      diffs <- scores[[paste0("differential_", name, "_vs_controls")]]
      sig <- scores[[paste0("significant_", name, "_vs_controls")]]
      scores[[response_col]][sig & diffs < 0 & abs(diffs) > differential_threshold] <- neg_type
      scores[[response_col]][sig & diffs > 0 & abs(diffs) > differential_threshold] <- pos_type
    }
  }
  
  # Look into chromosomal correction by removing interactions on the X chromosome
  
  # Explicitly returns scored data
  return(scores)
}

#' Scores multiple drugs against multiple controls
#' 
#' Takes in an input .tsv file with two columns for "Screen" and "Control" and scores
#' all screens listed in "Screen" against their corresponding screens listed in
#' "Control." Outputs all files and plots in the specified folder. 
#' 
#' @param df LFC dataframe.
#' @param screens List of screens generated with \code{add_screens}.
#' @param batch_file Path to .tsv file mapping screens to their controls for scoring, with two 
#'   columns for "Screen" and "Control." Screens to score against dervied null-models with the
#'   combn scoring mode must have their respective control labeled as "combn." 
#' @param output_folder Folder to output scored data and plots to. 
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ma_transform If true, M-A transforms data before running loess normalization. Only
#'   has an effect when loess = TRUE (default TRUE).
#' @param control_genes List of control genes to remove, e.g. "luciferase" (default c("None", "")).
#' @param n_components Vector containing indices of principal components to remove from data.  
#' @param chromosomal_correction If TRUE, corrects chromosomal shifts by down-weighting qGI
#'   scores for shifted regions (default FALSE).
#' @param weight_method If "linear" weights non-matched control screens equally. If exponential,
#'   instead weights non-matched control screens according to a decreasing exponential function
#'   based on the similarity of control screens to the given condition screen.
#' @param matched_fraction The weight given to the matched control as a fraction of 1, where
#'   all non-matched controls receive a total weight equal to 1 - matched_fraction (default 0.75).
#' @param sd_scale_factor Factor to normalize SDs against for scaling. If NULL for group scoring,
#'   the mean is computed across guide-level residuals, otherwise the given scalar is used instead.
#'   If NULL for one-off scoring, this operation is not performed (default NULL).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default
#'   "BY")
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0.5).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @param label_fdr_threshold Threshold below which to plot gene labels for significant
#'   hits, or NULL to plot without labels (default NULL).
#' @param save_residuals If true, saves residuals for each screen to the output folder
#'   (default FALSE).
#' @param plot_residuals If true, plots residual effects for all top hits (default TRUE).
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @param verbose If true, prints verbose output (default FALSE). 
#' @export
score_drugs_batch <- function(df, screens, batch_file, output_folder, 
                              min_guides = 3, test = "moderated-t", 
                              loess = TRUE, ma_transform = TRUE,
                              control_genes = c("None", ""), n_components = 0, 
                              chromosomal_correction = FALSE, weight_method = "exp",
                              matched_fraction = 0.75, sd_scale_factor = NULL,
                              fdr_method = "BY", fdr_threshold = 0.1, 
                              differential_threshold = 0.5, neg_type = "Negative", 
                              pos_type = "Positive", label_fdr_threshold = NULL,
                              save_residuals = FALSE, plot_residuals = TRUE, 
                              plot_type = "png", verbose = FALSE) {
  
  # Creates output folder if nonexistent
  if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }
  
  # Checks format of batch file, decides which version of scoring to run, and loads it
  scoring_type <- "Group"
  first_file <- utils::read.table(file = batch_file, header = F, nrows = 1, sep = "\t", encoding = "UTF-8")
  batch <- NULL
  if (ncol(first_file) == 4) {
    check_group_file(batch_file, screens)
  } else if (ncol(first_file) == 2) {
    check_batch_file(batch_file, screens)
    scoring_type <- "Single"
  } else {
    stop(paste("file", batch_file, "must contain exactly 2 or 4 columns (see score_drugs_batch documentation)"))
  }
  batch <- utils::read.csv(batch_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, encoding = "UTF-8")
  
  # Scores guides for each batch
  if (scoring_type == "Single") {
    all_scores <- NULL
    for (i in 1:nrow(batch)) {
      
      # Makes output folders if nonexistent
      lfc_folder <- file.path(output_folder, "lfc")
      plot_folder <- file.path(output_folder, "plots")
      if (!dir.exists(lfc_folder)) { dir.create(lfc_folder) }
      if (!dir.exists(plot_folder)) { dir.create(plot_folder) }
      
      # Scores each screen separately
      condition <- batch[i,1]
      control <- batch[i,2]
      temp <- score_drugs_vs_control(df, screens, control, condition, test = test, 
                                     min_guides = min_guides, 
                                     loess = loess, 
                                     ma_transform = ma_transform, 
                                     control_genes = control_genes, 
                                     fdr_method = fdr_method, 
                                     sd_scale_factor = sd_scale_factor,
                                     verbose = verbose)
      scores <- temp[["scored_data"]]
      residuals <- temp[["residuals"]]
      scores <- call_drug_hits(scores, control, condition,
                               neg_type = neg_type, pos_type = pos_type,
                               fdr_threshold = fdr_threshold, 
                               differential_threshold = differential_threshold)
      plot_drug_response(scores, 
                         control_name = control, 
                         condition_name = condition, 
                         output_folder = plot_folder,
                         neg_type = neg_type, 
                         pos_type = pos_type,
                         plot_type = plot_type, 
                         label_fdr_threshold = label_fdr_threshold)
      if (plot_residuals) {
        if (!is.na(residuals)) {
          plot_drug_residuals(scores, residuals, control, condition, lfc_folder, 
                              neg_type = neg_type, pos_type = pos_type,
                              plot_type = plot_type) 
        } else {
          cat("WARNING: residuals are set to NA, skipping residual plots\n")
        }
      }
      if (save_residuals) {
        if (!is.na(residuals)) {
          residuals_file <- paste0(condition, "_vs_", control, "_residuals.tsv")
          utils::write.table(residuals, file.path(output_folder, residuals_file), sep = "\t",
                             row.names = FALSE, col.names = TRUE, quote = FALSE) 
        } else {
          cat("WARNING: residuals are set to NA, skipping writing residuals to file\n")
        }
      }
      if (is.null(all_scores)) {
        all_scores <- scores
      } else {
        all_scores <- cbind(all_scores, scores[,3:ncol(scores)])
      }
    }
    if (!is.null(all_scores)) {
      utils::write.table(all_scores, file.path(output_folder, "condition_gene_calls.tsv"), sep = "\t",
                         row.names = FALSE, col.names = TRUE, quote = FALSE) 
    } 
  } else if (scoring_type == "Group") {
    
    # Gets unique groups from batch file and scores each group separately
    groups <- unique(batch$Group)
    for (group in groups) {
      
      # Sets path to intermediate file for current group
      intermediate_file <- file.path(output_folder, group, "intermediate.rda")
      
      # Scores dataset separately depending on number of PCs to remove
      for (i in 1:length(n_components)) {
        components <- n_components[i]
        
        # Makes output folders if nonexistent for each principal component to remove
        group_folder <- file.path(output_folder, group, paste0("PC_", components))
        lfc_folder <- file.path(group_folder, "lfc")
        plot_folder <- file.path(group_folder, "plots")
        if (!dir.exists(group_folder)) { dir.create(group_folder, recursive = TRUE) }
        if (!dir.exists(lfc_folder)) { dir.create(lfc_folder) }
        if (!dir.exists(plot_folder)) { dir.create(plot_folder) }
        
        # Stores intermediate output for first run and loads for subsequent runs
        load_intermediate <- TRUE
        if (i == 1) {
          load_intermediate <- FALSE
        }
        
        # Scores data for the current group
        group_ind <- batch$Group %in% group
        controls <- batch$Screen[group_ind & batch$Type == "control"]
        conditions <- batch$Screen[group_ind & batch$Type == "condition"]
        matched_controls <- batch$Control[group_ind & batch$Type == "condition"]
        temp <- score_drugs_vs_controls(df, screens, controls, conditions, matched_controls,
                                        plot_folder, 
                                        min_guides = min_guides, 
                                        loess = loess, 
                                        ma_transform = ma_transform, 
                                        control_genes = control_genes, 
                                        fdr_method = fdr_method,
                                        weight_method = weight_method,
                                        matched_fraction = matched_fraction,
                                        sd_scale_factor = sd_scale_factor,
                                        n_components = components,
                                        chromosomal_correction = chromosomal_correction,
                                        return_residuals = FALSE,
                                        intermediate_file = intermediate_file,
                                        load_intermediate = load_intermediate,
                                        verbose = verbose)
        scores <- temp[["scored_data"]]
        control_scores <- temp[["scored_controls"]]
        sd_table <- temp[["sd_table"]]
        scores <- call_drug_hits(scores, NULL, conditions,
                                 neg_type = neg_type, pos_type = pos_type,
                                 fdr_threshold = fdr_threshold, 
                                 differential_threshold = differential_threshold)
        for (condition in conditions) {
          plot_drug_response(scores, 
                             control_name = NULL, 
                             condition_name = condition, 
                             output_folder = plot_folder, 
                             neg_type = neg_type, 
                             pos_type = pos_type, 
                             plot_type = plot_type,
                             label_fdr_threshold = label_fdr_threshold) 
        }
        utils::write.table(scores, file.path(group_folder, "gene_calls.tsv"), sep = "\t",
                           row.names = FALSE, col.names = TRUE, quote = FALSE) 
        utils::write.table(control_scores, file.path(group_folder, "control_guide_effects.tsv"), sep = "\t",
                           row.names = FALSE, col.names = TRUE, quote = FALSE) 
        utils::write.table(sd_table, file.path(group_folder, "sd_table.tsv"), sep = "\t",
                           row.names = FALSE, col.names = TRUE, quote = FALSE) 
      }
    }
  }
}