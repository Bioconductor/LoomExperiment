\name{LoomGraph}
\alias{LoomGraph}
\alias{LoomGraph-class}

\docType{class}

\title{LoomGraph class}
\description{A description of the LoomGraph class for storing multiple LoomGraph objects.}

\usage{
LoomGraph(...)
}

\arguments{
    \item{...}{Arguments to pass to the \code{DataFrame} constructor.}
}

\value{
A LoomGraph object is returned from the constructor.
}

\details{
The LoomGraph class extends the \linkS4class{DataFrame} class.
It is meant to store graph information and was created to be used in conjunction with the \code{LoomExperiment} class.
Its purpose is to store multiple \code{col_graph} or \code{row_graph} as specified by the \code{loom} file format.
A LoomGraph object must contain two or three columns labeled \code{a}, \code{b}, and \code{w}.
Columns \code{a} and \code{b} indicate an edge between two vertices.
The \code{w} column indicates the weight of the corresponding edge and is optional.
All values in a LoomGraph must be \code{numeric}.
}

\seealso{
\code{\link{LoomExperiment}},
\code{\link{LoomGraphs}}
}

\author{
Daniel Van Twisk
}

\examples{
lg1 <- LoomGraph(a=c(1, 2, 3), b=(3, 2, 1), w=c(4, 7, 8))
lg1
}

