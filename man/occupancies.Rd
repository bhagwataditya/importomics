% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_read_proteingroups.R
\name{occupancies}
\alias{occupancies}
\alias{occupancies,SummarizedExperiment-method}
\alias{occupancies<-}
\alias{occupancies<-,SummarizedExperiment,matrix-method}
\alias{occupancies<-,SummarizedExperiment,numeric-method}
\title{Get/Set occupancies}
\usage{
occupancies(object)

\S4method{occupancies}{SummarizedExperiment}(object)

occupancies(object) <- value

\S4method{occupancies}{SummarizedExperiment,matrix}(object) <- value

\S4method{occupancies}{SummarizedExperiment,numeric}(object) <- value
}
\arguments{
\item{object}{SummarizedExperiment}

\item{value}{occupancy matrix (features x samples)}
}
\value{
occpuancy matrix (get) or updated object (set)
}
\description{
Get / Set phosphosite occupancies matrix
}
\examples{
file <- download_data('fukuda20.proteingroups.txt')
object <- read_proteingroups(file, plot=FALSE)
occupancies(object)
occupancies(object) <- exprs(object)
occupancies(object)[1:3, 1:3]
}
