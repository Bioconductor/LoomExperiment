\name{LoomGraphs}
\alias{LoomGraphs}
\alias{LoomGraphs-class}

\docType{class}

\title{LoomGraphs class}
\description{
The LoomGraphs class extends the \code{\link[S4Vectors:SimpleList-class]{SimpleList}} class.
It is meant to store multiple \code{LoomGraph} objects and was created to be used in conjunction with the \code{LoomExperiment} class.
Its purpose is to store multiple \code{col_graph} or \code{row_graph} as specified by the \code{loom} file format.
Only \code{LoomGraph} obects may be stored in a LoomGraphs object.
}

\usage{
LoomGraphs(...)
}

\arguments{
    \item{...}{\code{LoomGraph} objects.}
}

\value{
A LoomGraphs object is returned from the constructor.
}

\seealso{
\code{\link{LoomExperiment}},
\code{\link{LoomGraph}},
\code{\link[S4Vectors:SimpleList]{SimpleList}}
}

\examples{

## Construction
lg1 <- LoomGraph(a=c(1, 2, 3), b=c(3, 2, 1), w=c(4, 7, 8))
lg2 <- LoomGraph(a=c(3, 3, 1), b=c(3, 1, 2))
lgs <- LoomGraphs(lg1, lg2)
lgs

## Subsetting
lgs[c(1, 2)]
}

\author{
Daniel Van Twisk
}

