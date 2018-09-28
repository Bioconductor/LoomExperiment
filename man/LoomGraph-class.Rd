\name{LoomGraph}

%Class
\alias{LoomGraph}
\alias{LoomGraph-class}

%SelectHits
\alias{selectHits}
\alias{selectHits,LoomGraph-method}
\alias{selectHits,LoomGraph,ANY-method}

%DropHits
\alias{dropHits}
\alias{dropHits<-}
\alias{dropHits,LoomGraph-method}
\alias{dropHits,LoomGraph,ANY-method}
\alias{dropHits<-,LoomGraph-method}
\alias{dropHits<-,LoomGraph,ANY,ANY-method}

\docType{class}

\title{LoomGraph class}
\description{
The LoomGraph class extends the code{\link[S4Vectors:Hits-class]{SelfHits}} class.
It is meant to store graph information and was created to be used in conjunction with the \code{LoomExperiment} class.
Its purpose is to store multiple \code{col_graph} or \code{row_graph} as specified by the \code{loom} file format.
Attributes \code{from} and \code{to} indicate an edge between two vertices.
The \code{w} column indicates the weight of the corresponding edge and is optional.
}

\usage{
LoomGraph(from, to, nnode=max(from, to), ..., weight=NULL)

## SelectHits
\S4method{selectHits}{LoomGraph}(x, i, ...)

## DropHits
\S4method{dropHits}{LoomGraph}(x, i, ...)
\S4method{dropHits}{LoomGraph}(x, i, ...) <- value

}

\arguments{

  \item{x}{A LoomGraph object}

  \item{i}{Indices specifying elements to subset, select elements, or drop elements from LoomGraph by}

  \item{...}{Arguments to pass to the \code{SelfHits} constructor.}

  \item{from}{A numeric vector of nodes indicating one side of the graph's edge.}

  \item{to}{A numeric vector of node indicating the second side of the graph's edge.}

  \item{nnode}{An integer indicating the maximum number of nodes in the graph.}

  \item{weight}{A numeric vector indicating the weight between the prospective edges.}

  \item{value}{For subsetting methods, a \code{SelfHits} to replace the selected entry.
               For dropHits methods, a \code{numeric} to replace selected nodes.}
}

\value{
A LoomGraph object is returned from the constructor.
}

\seealso{
\code{\link{LoomExperiment}},
\code{\link{LoomGraphs}},
\code{\link[S4Vectors:DataFrame]{DataFrame}}
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

## SelectHits
selectHits(lg, c(1))

## DropHits
dropHits(lg, c(2))
dropHits(lg, c(1)) <- c(7)
}

\author{
Daniel Van Twisk
}
