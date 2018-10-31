\name{LoomFile}
\alias{LoomFile}
\alias{LoomFile-class}

\docType{class}

\title{LoomFile objects}

\description{
    A \code{LoomFile} class represents a loom file based on the Linnarson Lab's
    \href{Loom File Format}{http://linnarssonlab.org/loompy/index.html}.  A loom
    file is encoded as an hdf5 file. A loom file consists of a main matrix,
    optional additional layers, a variable number of row and column annotations
    and sparse graph objects.  It is used to efficienctly stroe very large omics
    datasets.

    The \code{LoomFile} class extends the functionality of the
    \code{\link[rtracklayer:RTLFile-class]{RTLFile}} from the \code{rtracklayer}
    package.
}

\author{
Daniel Van Twisk
}
