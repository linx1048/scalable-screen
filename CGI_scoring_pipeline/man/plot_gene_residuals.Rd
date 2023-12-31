% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_gene_residuals}
\alias{plot_gene_residuals}
\title{Plot guide-level residuals for a given gene.}
\usage{
plot_gene_residuals(scores, residuals, gene, control_name, condition_name)
}
\arguments{
\item{scores}{Dataframe of scores returned from \code{call_drug_hits}.}

\item{residuals}{Residuals returned with the return_residuals argument set to true
from \code{call_drug_hits}.}

\item{gene}{Gene name for guides to plot.}

\item{control_name}{Name of control passed to \code{call_drug_hits}.}

\item{condition_name}{Name of condition passed to \code{call_drug_hits}.}
}
\value{
A ggplot object.
}
\description{
Plots guide-level residuals for a given gene and returns the plot. Works for data 
returned from \code{call_drug_hits}.
}
