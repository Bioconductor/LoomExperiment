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
  \item{object}{A \code{LoomExperiment} object to be exported. File must have the \code{.loom} extension.}
  \item{con}{
    The connection from which data is loaded or to which data is saved. If this is a \code{character vector}, it is assumed to be a filename and a corresponding file connection is created and then closed after exporting the object. If a \code{RTLFile} derivative, the data is loaded from or saved to the underlying resource.
  }
  \item{matrix}{A \code{matrix} in which the column and rows for the Loom file will be derived. The default argument is dervied from the the first assay in the \code{LoomExperiment} object.}
  \item{rownames_attr}{A \code{character vector} indicating the name of the \code{rowData} attribute to represent the names \code{rownames} in the \code{LoomExperiment} object.}
  \item{colnames_attr}{A \code{character vector} indicating the name of the \code{colData} attribute to represent the names \code{colnames} in the \code{LoomExperiment} object.}
}

\value{
  An error code indiciating whether the operation was successful.
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
