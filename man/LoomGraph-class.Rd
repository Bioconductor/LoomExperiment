\name{LoomGraph}

%Class
\alias{LoomGraph}
\alias{LoomGraph-class}
\alias{rbind,LoomGraph-method}
\alias{cbind,LoomGraph-method}

\docType{class}

\title{LoomGraph class}
\description{
The \code{LoomGraph} class extends the \code{\link[S4Vectors:Hits-class]{SelfHits}} class.
The \code{SelfHits} class represents a set of hits between a set of left node and right nodes.
Only the Hits are stored in a \code{SelfHits object}.
The \code{LoomGraph} class is meant to store graph information and was created to be used in conjunction with the \code{LoomExperiment} class.
Its purpose is to store either a \code{col_graph} or \code{row_graph} as specified by the \code{loom} file format.
Attributes \code{from} and \code{to} indicate an edge between two vertices.
The \code{w} column indicates the weight of the corresponding edge and is optional.
}

\usage{
LoomGraph(from, to, nnode=max(from, to), ..., weight=NULL)

\S4method{rbind}{LoomGraph}(..., deparse.level=1)
\S4method{cbind}{LoomGraph}(..., deparse.level=1)
}

\arguments{

  \item{...}{Arguments to pass to the \code{SelfHits} constructor.}

  \item{from}{A numeric vector of nodes indicating one side of the graph's edge.}

  \item{to}{A numeric vector of node indicating the second side of the graph's edge.}

  \item{nnode}{An integer indicating the maximum number of nodes in the graph.}

  \item{weight}{A numeric vector indicating the weight between the prospective edges.}

  \item{deparse.level}{See '?base::cbind' for a description of this argument}
}

\value{
A LoomGraph object is returned from the constructor.
}

\seealso{
\code{\link{LoomExperiment}},
\code{\link{LoomGraphs}},
\code{\link[S4Vectors:Hits-class]{SelfHits}}
}

\examples{
## Construction
a <- c(1, 2, 3)
b <- c(3, 2, 1)
w <- c(100, 10, 1)
df <- DataFrame(a, b, w)
lg <- as(df, "LoomGraph")
# OR
lg <- LoomGraph(a, b, weight=w)
lg

## Subsetting
lg[c(1, 2)]
lg[-c(2)]
}

\author{
Daniel Van Twisk
}
