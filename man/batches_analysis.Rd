% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flowSim_similarity.R
\name{batches_analysis}
\alias{batches_analysis}
\title{batches_analysis}
\usage{
batches_analysis(all_groups, n_cores, final_df, path_dir, thr_score,
  progress_bar = T, nboots)
}
\arguments{
\item{all_groups}{All batches names to analyze}

\item{n_cores}{Number to use}

\item{final_df}{Dataframe of features reduced}

\item{path_dir}{Path to input directory. Mandatory.}

\item{thr_score}{Threshold for filtering the similarity scores}

\item{progress_bar}{If False,disable progress bar.}

\item{nboots}{Number of bootstrap simulations for the  density identity test. Default to 6}
}
\value{
List of results batch analysis
}
\description{
function to analyze each batch. A batch is a fragment of the original dataset. Internal function.
}
\keyword{Internal}
