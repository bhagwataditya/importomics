% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_read_somascan.R
\name{.read_somascan}
\alias{.read_somascan}
\alias{read_somascan}
\title{Read somascan}
\usage{
.read_somascan(
  file,
  fidvar = "SeqId",
  sidvar = "SampleId",
  subgroupvar = "SampleGroup"
)

read_somascan(
  file,
  fidvar = "SeqId",
  sidvar = "SampleId",
  subgroupvar = "SampleGroup",
  fname_var = "EntrezGeneSymbol",
  sample_type = "Sample",
  feature_type = "Protein",
  sample_quality = c("FLAG", "PASS"),
  feature_quality = c("FLAG", "PASS"),
  rm_na_svars = FALSE,
  rm_single_value_svars = FALSE,
  pca = FALSE,
  limma = FALSE,
  formula = NULL,
  block = NULL,
  contrastdefs = NULL,
  verbose = TRUE,
  plot = TRUE
)
}
\arguments{
\item{file}{*.adat file path (string)}

\item{fidvar}{featureid fvar (string)}

\item{sidvar}{sampleid svar (string)}

\item{subgroupvar}{subgroup svar (string)}

\item{fname_var}{featurename fvar (string)}

\item{sample_type}{subset of c('Sample','QC','Buffer','Calibrator')}

\item{feature_type}{subset of c('Protein',
'Hybridization Control Elution',
'Rat Protein')}

\item{sample_quality}{subset of c('PASS', 'FLAG', 'FAIL')}

\item{feature_quality}{subset of c('PASS', 'FLAG', 'FAIL')}

\item{rm_na_svars}{whether to rm NA svars}

\item{rm_single_value_svars}{whether to rm single value svars}

\item{pca}{whether to pca}

\item{limma}{whether to limma/contrast}

\item{formula}{design formula (using svars)}

\item{block}{block var}

\item{contrastdefs}{contrastdef vector/matrix/list}

\item{verbose}{whether to msg}

\item{plot}{whether to plot}
}
\value{
Summarizedexperiment
}
\description{
Read data from somascan adat file
}
\examples{
file <- download_data('atkin18.somascan.adat')
read_somascan(file, pca = TRUE, limma = TRUE, block = 'Subject_ID')
}
