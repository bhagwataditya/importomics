% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_read_omics.R
\name{merge_sdata}
\alias{merge_sdata}
\alias{merge_fdata}
\title{Merge sample/feature data}
\usage{
merge_sdata(
  object,
  dt,
  by.x = "sample_id",
  by.y = names(dt)[1],
  subgroupvar = NULL,
  verbose = TRUE
)

merge_fdata(
  object,
  dt,
  by.x = "feature_id",
  by.y = names(dt)[1],
  fnamevar = NULL,
  verbose = TRUE
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{dt}{data.frame, data.table, DataFrame}

\item{by.x}{object mergevar}

\item{by.y}{df mergevar}

\item{subgroupvar}{subgroup svar}

\item{verbose}{TRUE/FALSE}

\item{fnamevar}{featurename fvar}
}
\value{
SummarizedExperiment
}
\description{
Merge sample/feature data
}
\examples{
require(magrittr)
file <- download_data('halama18.metabolon.xlsx')
object <- read_metabolon(file, plot=FALSE)
object \%<>\% merge_sdata( data.frame(sample_id = object$sample_id,
                                    number = seq_along(object$sample_id)))
head(sdata(object))
}
