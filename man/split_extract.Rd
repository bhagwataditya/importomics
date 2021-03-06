% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_contrast.R
\name{nfactors}
\alias{nfactors}
\alias{split_extract}
\title{stri_split and extract}
\usage{
nfactors(x, sep = guess_sep(x))

split_extract(x, i, sep = guess_sep(x))
}
\arguments{
\item{x}{string}

\item{sep}{string}

\item{i}{integer}
}
\value{
character
}
\description{
stri_split and extract
}
\examples{
require(magrittr)
file <- download_data('halama18.metabolon.xlsx')
object <- read_metabolon(file, plot=FALSE)
x <- object$sample_id[1:5]
nfactors(x)
split_extract(x, 1:2)
split_extract(x, seq_len(nfactors(x)-1))
split_extract(x, nfactors(x))
}
