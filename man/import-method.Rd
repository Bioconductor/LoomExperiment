\name{import-methods}

\docType{methods}

\alias{import}
\alias{import-methods}
\alias{import,LoomFile,ANY,ANY-method}

\title{Import LoomExperiment to loomFile}

\description{
    Imports a \code{LoomExperiment} from a \code{LoomFile}.
}

\usage{
\S4method{import}{LoomFile,ANY,ANY}(con, ...,
    type = c("SingleCellLoomExperiment", "SummarizedLoomExperiment", "RangedLoomExperiment"),
    rownames_attr=NULL, colnames_attr=NULL)
}

\arguments{
  \item{con}{A character indicating the loom file to be created. File must have the \code{.loom} extension.}
  \item{...}{Additional arguments}
  \item{type}{
    Either "SingleCellLoomExperiment", "SummarizedLoomExperiment", or "RangedLoomExperiment".
    This value decides what type of object that will be returned by \code{import}.  If left empty
    \code{import} will either determine what type of class should be used by the context of the
    file.  If it cannot be determined, the \code{LoomExperiment} type will default to
    \code{SingleCellLoomExperiment}.
  }
  \item{rownames_attr}{A \code{character vector} indicating the name of the \code{rowData} attribute to represent the names \code{rownames} in the \code{LoomExperiment} object.}
  \item{colnames_attr}{A \code{character vector} indicating the name of the \code{colData} attribute to represent the names \code{colnames} in the \code{LoomExperiment} object.}
}

\value{
  An object of class \code{LoomExperiment}
}

\examples{
l1_file <- system.file("extdata", "L1_DRG_20_example.loom", package = "LoomExperiment")
scle <- import(l1_file, type="SingleCellLoomExperiment")
scle
}

\seealso{
  \code{\linkS4class{LoomExperiment}},
  \code{\linkS4class{LoomFile}},
}
