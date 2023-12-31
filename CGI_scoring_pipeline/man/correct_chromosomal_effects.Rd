% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.R
\name{correct_chromosomal_effects}
\alias{correct_chromosomal_effects}
\title{Corrects for GIs shifted in specific chromosomes.}
\usage{
correct_chromosomal_effects(df, guide_df)
}
\arguments{
\item{df}{LFC dataframe passed into \code{score_drugs_vs_controls} with chromosomal location 
specified in a column named "chr", start coordinates specified in "start_loc" and stop 
coordinates specified in "stop_loc."}

\item{guide_df}{Dataframe of guide-level qGI scores.}
}
\value{
Dataframe of chromosomal-corrected, guide-level qGI scores.
}
\description{
Run this within \code{score_drugs_vs_controls} to implement qGI chromosomal correction 
steps.
}
