\name{LoomGraph}
\alias{LoomGraph}
\alias{LoomGraph-class}
\alias{[,LoomGraph,ANY,missing,ANY-method}
\alias{[<-,LoomGraph,ANY,missing,numeric-method}

\docType{class}

\title{LoomGraph class}
\description{A description of the LoomGraph class for storing multiple LoomGraph objects.}

\usage{
LoomGraph(...)

## Subsetting

\S4method{[}{LoomGraph,ANY,missing,ANY}(x, i, j, ..., drop=TRUE)
\S4method{[}{LoomGraph,ANY,missing,numeric}(x, i, j, ...) <- value
}

\arguments{

  \item{x}{A LoomGraph object}

  \item{i, j}{Indices specifying elements to subset LoomGraph by}

  \item{...}{Arguments to pass to the \code{DataFrame} constructor.}

  \item{drop}{
    For matrices and arrays.  If 'TRUE' the result is coerced to the lowest
    posible dimesnion.  This only woeks for extracting elements, not for the
    replacement.
  }

  \item{value}{A \code{LoomGraph} to replace the selected entry}
}

\value{
A LoomGraph object is returned from the constructor.
}

\details{
The LoomGraph class extends the code{\link[S4Vectors:DataFrame-class]{DataFrame}} class.
It is meant to store graph information and was created to be used in conjunction with the \code{LoomExperiment} class.
Its purpose is to store multiple \code{col_graph} or \code{row_graph} as specified by the \code{loom} file format.
A LoomGraph object must contain two or three columns labeled \code{a}, \code{b}, and \code{w}.
Columns \code{a} and \code{b} indicate an edge between two vertices.
The \code{w} column indicates the weight of the corresponding edge and is optional.
All values in a LoomGraph must be \code{numeric}.
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
lg <- LoomGraph(df)
# OR
lg <- as(df, "LoomGraph")
# OR
lg <- LoomGraph(a, b, w)
lg

## Subsetting
lg[c(1, 2)]
lg[-c(2)]
}

\author{
Daniel Van Twisk
}
