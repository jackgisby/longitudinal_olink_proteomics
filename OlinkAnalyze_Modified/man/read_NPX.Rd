% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_npx_data.R
\name{read_NPX}
\alias{read_NPX}
\title{Function to read NPX data into long format}
\usage{
read_NPX(
  filename,
  sample_manifest = NULL,
  pheno = NULL,
  skip_mod = 0,
  panel = "NA",
  tab = 1,
  this_exp = "plasma"
)
}
\arguments{
\item{filename}{Path to file NPX Manager output file.}

\item{sample_manifest}{Per individual phenotypic data}

\item{pheno}{Per sample phenotypic data}

\item{skip_mod}{How many lines to skip in the excel file}

\item{panel}{Which panel is currently being read}

\item{tab}{Which excel tab the data is being read from}

\item{this_exp}{Whether the data is plasma or serum}
}
\value{
A tibble in long format.
}
\description{
Modifications have been made to allow the reading of our specific dataset and
add GeneIDs in addition to the olink assay names.
}
\examples{
\donttest{read_NPX("~/NPX data.xlsx")}
}
\keyword{NPX}
