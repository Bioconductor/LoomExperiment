\name{export-methods}

\docType{methods}

\alias{export}
\alias{export-methods}
\alias{export,LoomExperiment,LoomFile,ANY-method}

\title{Export LoomExperiment to loomFile}

\description{
    Exports a \code{LoomExperiment} to a \code{LoomFile}.
}

\usage{
\S4method{export}{LoomExperiment,LoomFile,ANY}(object, con,
    matrix=assayNames(object)[1], rownames_attr="rownames", colnames_attr="colnames")
}

\arguments{
  \item{object}{A LoomExperiment object to be exported}
  \item{con}{A character indicating the loom file to be created. File must have the \code{.loom} extension.}
  \item{matrix}{A matrix in which the column and rows for the \code{.loom} file will be derived.}
  \item{rownames_attr}{A character indicating the row names}
  \item{colnames_attr}{A character indicating the col names}
}

\value{
  An error code indiciating whether the operation was successful
}

\examples{

counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(assays = list(counts = counts))
scle <- SingleCellLoomExperiment(sce)
tempfile <- tempfile(fileext=".loom")
export(scle, tempfile)
}

\seealso{
  \code{\linkS4class{LoomExperiment}},
  \code{\linkS4class{LoomFile}},
}
