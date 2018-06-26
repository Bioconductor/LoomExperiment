\name{LoomGraphs}
\alias{LoomGraphs}
\alias{LoomGraphs-class}

\docType{class}

\title{LoomGraphs class}
\description{A description of the LoomGraphs class for storing multiple LoomGraph objects.}

\usage{
LoomGraphs(...)
}

\arguments{
    \item{...}{Arguments to pass to the \code{\link[S4Vectors:SimpleList-class]{SimpleList}} constructor.}
}

\value{
A LoomGraphs object is returned from the constructor.
}

\details{
The LoomGraphs class extends the \code{\link[S4Vectors:SimpleList-class]{SimpleList}} class.
It is meant to store multiple \code{LoomGraph} objects and was created to be used in conjunction with the \code{LoomExperiment} class.
Its purpose is to store multiple \code{col_graph} or \code{row_graph} as specified by the \code{loom} file format.
Only \code{LoomGraph} obects may be stored in a LoomGraphs object.
}

\seealso{
\code{\link{LoomExperiment}},
\code{\link{LoomGraph}}
}

\author{
Daniel Van Twisk
}

\examples{
lg1 <- LoomGraph(a=c(1, 2, 3), b=c(3, 2, 1), w=c(4, 7, 8))
lg2 <- LoomGraph(a=c(3, 3, 1), b=c(3, 1, 2))
lgs <- LoomGraphs(lg1, lg2)
lgs
}

