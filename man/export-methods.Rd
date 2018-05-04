\name{export-methods}

\docType{methods}

\alias{export}
\alias{export-methods}
\alias{export,LoomExperiment,loomFile,"ANY"}

\title{Export LoomExperiment to loomFile}

\description{
    Exports a \code{LoomExperiment} to a \code{loomFile}.
}

\usage{
export(object, con, matrix=assayNames(object)[1], rownames_attr="rownames",
       colnames_attr="colnames")
}
\arguments{
  \item{object}{A LoomExperiment object to be exported}
  \item{con}{A character indicating the loom file to be created. File must have the \code{.loom} extension.}
  \item{matrix}{A matrix}
  \item{rownames_attr}{A character indicating the row names}
  \item{colnames_attr}{A character indicating the col names}
}

\seealso{
  \code{\linkS4class{LoomExperiment}}, \code{\linkS4class{loomFile}}
}
