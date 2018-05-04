\name{import-methods}

\docType{methods}

\alias{import}
\alias{import-methods}
\alias{import,LoomExperiment,loomFile,"ANY"}

\title{Export LoomExperiment to loomFile}

\description{
    Exports a \code{LoomExperiment} to a \code{loomFile}.
}

\usage{
import(con, ..., rownames_attr="rownames", colnames_attr="colnames")
}
\arguments{
  \item{con}{A character indicating the loom file to be created. File must have the \code{.loom} extension.}
  \item{...}{Additional arguments}
  \item{rownames_attr}{A character indicating the row names}
  \item{colnames_attr}{A character indicating the col names}
}

\seealso{
  \code{\linkS4class{LoomExperiment}}, \code{\linkS4class{loomFile}}
}
