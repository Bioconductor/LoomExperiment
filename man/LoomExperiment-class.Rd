
\name{LoomExperiment}

\docType{class}

% Class
\alias{LoomExperiment}
\alias{LoomExperiment-class}
\alias{RangedLoomExperiment}
\alias{RangedLoomExperiment-class}
\alias{SingleCellLoomExperiment}
\alias{SingleCellLoomExperiment-class}

% Accessors
\alias{colGraphs}
\alias{rowGraphs}
\alias{colGraphs<-}
\alias{rowGraphs<-}
\alias{colGraphs,LoomExperiment-method}
\alias{rowGraphs,LoomExperiment-method}
\alias{colGraphs<-,LoomExperiment-method}
\alias{rowGraphs<-,LoomExperiment-method}

% Subsetting
\alias{[,LoomExperiment-method}
\alias{[,LoomExperiment,ANY-method}
\alias{[,LoomExperiment,ANY,ANY,ANY-method}
\alias{[,RangedLoomExperiment-method}
\alias{[,RangedLoomExperiment,ANY-method}
\alias{[,RangedLoomExperiment,ANY,ANY,ANY-method}
\alias{[,SingleCellLoomExperiment-method}
\alias{[,SingleCellLoomExperiment,ANY-method}
\alias{[,SingleCellLoomExperiment,ANY,ANY,ANY-method}
\alias{rbind,LoomExperiment-method}
\alias{rbind,SingleCellLoomExperiment-method}
\alias{cbind,LoomExperiment-method}
\alias{cbind,SingleCellLoomExperiment-method}

% Display
\alias{show,LoomExperiment-method}
\alias{show,RangedLoomExperiment-method}
\alias{show,SingleCellLoomExperiment-method}

\title{LoomExperiment, RangedLoomExperiment, and SingleCellLoomExperiment classes}

\description{
    The LoomExperiment family of classes is used as a bridge between Bioconductor's "Experiment" classes and the Linnarson Lab's \href{Loom file format}{http://linnarssonlab.org/loompy/index.html}.
    The family of LoomExperiment classes all inherit from the class LoomExperiment as well as their respectively named parent classes.
    The LoomExperiment class inheirts from \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
}

\usage{

## Constructor

LoomExperiment(..., colGraphs = LoomGraphs(), rowGraphs = LoomGraphs())
RangedLoomExperiment(..., colGraphs = LoomGraphs(), rowGraphs = LoomGraphs())
SingleCellLoomExperiment(..., colGraphs = LoomGraphs(), rowGraphs = LoomGraphs())

## Accessors

\S4method{colGraphs}{LoomExperiment}(x, ...)
\S4method{colGraphs}{LoomExperiment}(x, ...) <- value
\S4method{rowGraphs}{LoomExperiment}(x, ...)
\S4method{rowGraphs}{LoomExperiment}(x, ...) <- value

## Subsetting

\S4method{[}{LoomExperiment}(x, i, j, ..., drop=TRUE)

## Binding

\S4method{rbind}{LoomExperiment}(..., deparse.level=1)
\S4method{cbind}{LoomExperiment}(..., deparse.level=1)
}


\arguments{

  \item{x}{A LoomExperiment object}

  \item{colGraphs,rowGraphs}{
    LoomGraphs to be placed in either the colGraphs or rowGraphs slot respectively
  }

  \item{value}{
    For \code{colGraphs<-} and \code{rowGraphs<-}, \code{value} will be the replacement to the slot.
    For \code{dropHits<-}, indices to replace selected indices with.
  }

  \item{...}{
    For constructors, \code{...} will be passed on to the respective Experiment constructor.
  }

  \item{i,j}{
    For subsetting, indices specifying elements to subset LoomGraph by.
    For dropHits, numeric indicating the node number
  }

  \item{drop}{
    For matrices and arrays.  If 'TRUE' the result is coerced to the lowest 
    posible dimesnion.  This only woeks for extracting elements, not for the 
    replacement.
  }

  \item{deparse.level}{
    See '?base::cbind' for a description of this argument.
  }

}

\value{
  An object of class \code{LoomExperiment}
}

\details{
The LoomExperiment class is a virtual class meant to to act as an interface for other \code{"_LoomExperiment"} classes.
It contains two slots:
\describe{
\item{\code{colGraphs}:}{A LoomGraphs object containing \code{col_graph} data as specified by the loom format.}
\item{\code{rowGraphs}:}{A LoomGraphs object containing \code{row_graph} data as specified by the loom format.}
}

The intended use of this class is as an interface that allows various slots and operations necessary for subsequent \code{"_LoomExperiment"} classes to be defined.

The \code{colGraphs} and \code{rowGraphs} slot stores a \code{LoomGraphs} object that stores a graph of edges between vertices and possibly associated weights. These slots may be \code{NULL}.
}

\seealso{
\code{\link{SummarizedExperiment}},
\code{\link{RangedSummarizedExperiment}},
\code{\link{SingleCellExperiment}}
}

\examples{

## Construction
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(assays = list(counts = counts))
scle <- SingleCellLoomExperiment(sce)
# OR
scle <- SingleCellLoomExperiment(assays = list(counts = counts))
# OR
scle <- as(sce, "SingleCellLoomExperiment")
scle

## Get and replace rowGraphs and colGraphs
colGraphs(scle)
rowGraphs(scle)

a <- c(1, 2, 3)
b <- c(3, 2, 1)
w <- c(100, 10, 1)
lg <- LoomGraph(a, b, weight=w)
lgs <- LoomGraphs(lg, lg)
names(lgs) <- c('lg1', 'lg2')
lgs

colGraphs(scle) <- lgs
rowGraphs(scle) <- lgs

colGraphs(scle)
rowGraphs(scle)
colGraphs(scle)[[1]]
rowGraphs(scle)[[1]]

## Subsetting
scle2 <- scle[c(1, 3), 1:2]
colGraphs(scle2)[[1]]
rowGraphs(scle2)[[1]]
}

\author{
Daniel Van Twisk
}

