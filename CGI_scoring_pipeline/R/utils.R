######
# UTILITY FUNCTIONS
######

#' Adds a list of screens from a sample table file.
#' 
#' For a sample table formatted as a .tsv file with columns for screen, replicates,
#' and the name of a screen to compute log fold-changes against, adds each screen
#' to a list and returns that list.
#' 
#' @param table_file A tab-separated sample table with three columns, described above: 
#'   Screen, Replicates, and NormalizeTo. Screen is the name of the screen, replicates 
#'   are each technical replicate separated by semicolons, and NormalizeTo is either 
#'   the screen to normalize against or NA if unnecessary (e.g. for T0 screens).
#' @return A named list corresponding to provided screen names, where each sub-list 
#'   contains a list of the replicate columns (in "replicates") and the screen to 
#'   normalize against (in "normalize_name").
#' @export
add_screens_from_table <- function(table_file) {
  table <- utils::read.csv(table_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                           encoding = "UTF-8")
  
  # Checks columns of table
  cols <- c("Screen", "Replicates", "NormalizeTo")
  for (col in cols) {
    if (!(col %in% colnames(table))) {
      stop(paste("column", col, "not in sample table"))
    } else {
      if (col != "NormalizeTo") {
        if (any(is.na(table[[col]])))
        stop(paste("column", col, "contains NA values"))
      }
    }
  }
  
  # Adds all screens in table to list and returns
  screens <- NULL
  for (i in 1:nrow(table)) {
    name <- table[["Screen"]][i]
    replicates <- unlist(strsplit(table[["Replicates"]][i], ";"))
    normalize_name <- table[["NormalizeTo"]][i]
    if (is.na(normalize_name)) {
      normalize_name <- NULL
    }
    if (i == 1) {
      screens <- add_screen(name = name, replicates = replicates, normalize_name = normalize_name)
    } else {
      screens <- add_screen(screens, name, replicates, normalize_name)
    }
  }
  return(screens)
}

#' Adds a new screen to a list of screens
#' 
#' Makes a list containing info for a given screen and optionally appends it to a given list. 
#' 
#' @param screen_list An existing list of screens (optional, default NULL).
#' @param name A name for the screen (e.g. "RPE1_T18").
#' @param replicates A list of columns containing replicate data for the given screen.
#' @param normalize_name A name of another screen to normalize data in this screen to
#'   (e.g. "RPE1_T0", optional, default NULL). 
#' @param target_coverage Target coverage of the given screen (optional, default 100).
#' @return A named list corresponding to provided screen names, where each sub-list 
#'   contains a list of the replicate columns (in "replicates") and the screen to 
#'   normalize against (in "normalize_name").
#' @export
add_screen <- function(screen_list = NULL, name, replicates, 
                       normalize_name = NULL, target_coverage = 200) {
  
  # Checks arguments
  if (is.na(name) | !is.character(name)) {
    stop("name must be a string")
  } else if (any(is.na(replicates)) | length(replicates) < 1 | !is.character(replicates)) {
    stop("replicates must be a list of one or more column names containing technical replicate readcounts")
  } else if (any(is.na(target_coverage)) | is.null(target_coverage) | is.character(target_coverage)) {
    stop("target_coverage must be a numeric value (e.g. 100 for 100x coverage)")
  } else if (name %in% names(screen_list)) {
    stop(paste("screen", name, "already contained in screens - duplicate screen names are not allowed"))
  }
  
  # Adds arguments to list
  screen <- list()
  screen[["replicates"]] <- replicates
  screen[["normalize_name"]] <- normalize_name
  screen[["target_coverage"]] <- target_coverage
  if (!is.null(screen_list)) {
    screen_list[[name]] <- screen
  } else {
    screen_list <- list()
    screen_list[[name]] <- screen
  }
  return(screen_list)
}

#' Removes a screen from a list of screens
#' 
#' Removes a screen with a given name from a list of screens. If multiple 
#' screens match the provided name, they are all removed.
#' 
#' @param screen_list An existing list of screens.
#' @param name A name for the screen to remove (e.g. "RPE1_T18").
#' @return \code{screen_list} with the named screen removed.
#' @export
remove_screen <- function(screen_list, name) {
  
  # Checks input
  if (length(name) > 1) {
    stop("name must be a single screen name contained in screen_list")
  } else if (is.name(name) | is.null(name)) {
    stop("name must be a single screen name contained in screen_list")
  }
  
  # Checks that screen is in screen_list
  ind <- names(screen_list) == name
  if (sum(ind)) {
    screen_list <- screen_list[!ind]
  } else {
    warning(paste("screen", name, "not found in screen_list"))
  }
  
  # Re-formats screen_list as empty list instead of empty named list
  if (length(screen_list) == 0) {
    screen_list <- list()
    message("No screens remaining in list. Returning empty list")
  }
  return(screen_list)
}

#' Checks replicates
#' 
#' Checks to make sure that all replicates contained in a screens object returned from 
#' \code{add_screens_from_table} or \code{add_screen} are contained in a given dataframe.
#' Returns four vectors in a list, where the first contains columns missing from the 
#' screens but contained in the dataframe, the second contains columns contained in the
#' screens but missing from the dataframe, the third contains replicates duplicated
#' in multiple screens, and the last contains duplicated screen names.
#' 
#' @param df Reads or LFC dataframe.
#' @param screens List of screens created with \code{add_screens}.
#' @return List of up to four vectors, where "missing_from_screens" contains all screens 
#'  in df but not screens, "missing_from_df" contains all screens in screens but not df,
#'  "dupe_reps" contains names of all technical replicates mapped to multiple screens
#'  in screens, and "dupe_screens" contains duplicated screen names.
#' @export
check_replicates <- function(df, screens) {
  
  # Gets all screens in each category
  all_reps <- c()
  missing_screens <- c()
  missing_df <- c()
  dupe_reps <- c()
  dupe_screens <- c()
  for (screen in screens) {
    reps <- screen$replicates
    for (rep in reps) {
      if (!(rep %in% colnames(df))) {
        missing_df <- c(missing_df, rep)
      }
      if (rep %in% all_reps) {
        dupe_reps <- c(dupe_reps, rep)
      }
      all_reps <- c(all_reps, rep)
    }
  }
  missing_screens <- colnames(df)[!(colnames(df) %in% all_reps)]
  
  # Gets any duplicate screens
  screen_table <- table(names(screens))
  if (any(screen_table) > 1) {
    dupe_screens <- names(screen_table[screen_table > 1])
  }
  
  # Returns three vectors of screens in a list
  output <- list()
  output[["missing_from_screens"]] <- missing_screens
  output[["missing_from_df"]] <- missing_df
  output[["dupe_reps"]] <- dupe_reps
  output[["dupe_screens"]] <- dupe_screens
  return(output)
}

#' Checks input parameters
#' 
#' Checks to make sure that the given screens match the given dataframe.
#' 
#' @param df Reads or LFC dataframe.
#' @param screens List of screens created with \code{add_screens}.
#' @return TRUE.
#' @keywords internal
check_screen_params <- function(df, screens) {
  
  # Gets all replicates
  reps <- c()
  for (screen in screens) {
    reps <- c(reps, screen[["replicates"]])
  }
  
  # Checks input
  if (!all(reps %in% colnames(df))) {
    ind <- which(!(reps %in% colnames(df)))
    rep_name <- reps[ind[1]]
    screen_name <- ""
    for (name in names(screens)) {
      if (rep_name %in% screens[[name]][["replicates"]]) {
        screen_name <- name
      }
    }
    stop(paste("replicate", rep_name, "not in df, remove screen", screen_name, "with remove_screens"))
  }
  return(TRUE)
}

#' Checks a batch scoring file
#' 
#' Checks to make sure that the contents of the .tsv file and its formats are appropriate for
#' running batch scoring functions.
#' 
#' @param batch_file Path to .tsv file mapping screens to their controls for scoring, with two 
#'   columns for "Screen" and "Control." Screens to score against derived null-models with the
#'   combn scoring mode must have their respective control labeled as "combn." 
#' @return TRUE.
#' @keywords internal
check_batch_file <- function(batch_file, screens) {
  
  # Checks if file exists and is a .tsv file
  ext <- tools::file_ext(batch_file)
  if (ext != "tsv") {
    stop(paste("file", batch_file, "not a .tsv file"))
  }
  if (!file.exists(batch_file)) {
    stop(paste("file", batch_file, "does not exist at the specified path"))
  }
  
  # Loads file and checks its format
  df <- utils::read.csv(batch_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, encoding = "UTF-8")
  if (colnames(df)[1] != "Screen" | colnames(df)[2] != "Control" | length(colnames(df)) > 2) {
    stop(paste("file", batch_file, "does not contain exactly two columns labeled Screen and Control"))
  }
  
  # Checks that all screens are represented in the screens object
  for (screen in df$Screen) {
    if (!(screen %in% names(screens))) {
      stop(paste("screen", screen, "not in screens"))
    }
  }
  for (screen in df$Control) {
    if (!(screen %in% names(screens))) {
      stop(paste("screen", screen, "not in screens"))
    }
  }
  
  return(TRUE)
}

#' Checks a group file
#' 
#' Checks to make sure that the contents of the .tsv file and its formats are appropriate for
#' running group scoring functions.
#' 
#' @param group_file Path to .tsv file mapping screens to their controls for scoring, with four 
#'   columns for "Screen", "Control", "Group" and "Type." The "Group" column can contain any 
#'   meaningful label for each group, whereas the "Type" column must be labeled either "control" 
#'   or "condition." For each unique group label, all conditions in that group are scored
#'   against all controls.
#' @return TRUE.
#' @keywords internal
check_group_file <- function(group_file, screens) {
  
  # Checks if file exists and is a .tsv file
  ext <- tools::file_ext(group_file)
  if (ext != "tsv") {
    stop(paste("file", group_file, "not a .tsv file"))
  }
  if (!file.exists(group_file)) {
    stop(paste("file", group_file, "does not exist at the specified path"))
  }
  
  # Loads file and checks its format
  df <- utils::read.csv(group_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, encoding = "UTF-8")
  if (colnames(df)[1] != "Screen" | colnames(df)[2] != "Control" | colnames(df)[3] != "Group" | 
      colnames(df)[4] != "Type" | length(colnames(df)) > 4) {
    stop(paste("file", group_file, "does not contain exactly four columns labeled Screen, Control, Group and Type"))
  }
  
  # Checks that all screens are represented in the screens object
  for (screen in df$Screen) {
    if (!(screen %in% names(screens))) {
      stop(paste("screen", screen, "not in screens"))
    }
  }
  for (type in df$Type) {
    if (type != "control" & type != "condition") {
      stop(paste("type", type, "must be either 'control' or 'condition'"))
    }
  }
  
  return(TRUE)
}
