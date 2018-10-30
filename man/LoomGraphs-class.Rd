\name{LoomGraphs}

%Class
\alias{LoomGraphs}
\alias{LoomGraphs-class}
\alias{rbind,LoomGraphs-method}

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

\S4method{rbind}{LoomGraphs}(..., deparse.level=1)
}

\arguments{
    \item{...}{\code{LoomGraph} objects.}
    
    \item{deparse.level}{See '?base::cbind' for a description of this argument}
}

\value{
A LoomGraphs object is returned from the constructor.
}

\seealso{
\code{\link{LoomExperiment}},
\code{\link{LoomGraph}},
}

\examples{

## Construction
lg1 <- LoomGraph(c(1, 2, 3), c(3, 2, 1), weight=c(4, 7, 8))
lg2 <- LoomGraph(c(3, 3, 1), c(3, 1, 2))
lgs <- LoomGraphs(lg1, lg2)
lgs
}

\author{
Daniel Van Twisk
}

