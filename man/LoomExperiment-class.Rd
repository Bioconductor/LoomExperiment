
\name{LoomExperiment}

\docType{class}

% Class
\alias{LoomExperiment}
\alias{LoomExperiment-class}
\alias{SummarizedLoomExperiment}
\alias{SummarizedLoomExperiment-class}
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
\alias{[,SummarizedLoomExperiment-method}
\alias{[,SummarizedLoomExperiment,ANY-method}
\alias{[,SummarizedLoomExperiment,ANY,ANY,ANY-method}
\alias{[,RangedLoomExperiment-method}
\alias{[,RangedLoomExperiment,ANY-method}
\alias{[,RangedLoomExperiment,ANY,ANY,ANY-method}
\alias{[,SingleCellLoomExperiment-method}
\alias{[,SingleCellLoomExperiment,ANY-method}
\alias{[,SingleCellLoomExperiment,ANY,ANY,ANY-method}

% Display
\alias{show,SummarizedLoomExperiment-method}
\alias{show,RangedLoomExperiment-method}
\alias{show,SingleCellLoomExperiment-method}

\title{LoomExperiment, SummarizedLoomExperiment, RangedLoomExperiment, and SingleCellLoomExperiment class}

\description{

    The LoomExperiment family of classes is used as a bridge between Bioconductor's Experiment classes
    and the Linnarson Lab's Loom file format.
    The family of LoomExperiment classes all inherit from the virtual class LoomExperiment as well as their
    respectively named parent classes.
}

\usage{

## Constructor

SummarizedLoomExperiment(..., colGraphs = LoomGraphs(), rowGraphs = LoomGraphs())
RangedLoomExperiment(..., colGraphs = LoomGraphs(), rowGraphs = LoomGraphs())
SingleCellLoomExperiment(..., colGraphs = LoomGraphs(), rowGraphs = LoomGraphs())

## Accessors

\S4method{colGraphs}{LoomExperiment}(x, ...)
\S4method{colGraphs}{LoomExperiment}(x, ...) <- value
\S4method{rowGraphs}{LoomExperiment}(x, ...)
\S4method{rowGraphs}{LoomExperiment}(x, ...) <- value

## Subsetting

\S4method{[}{SummarizedLoomExperiment}(x, i, j, ..., drop=TRUE)
\S4method{[}{RangedLoomExperiment}(x, i, j, ..., drop=TRUE)
\S4method{[}{SingleCellLoomExperiment}(x, i, j, ..., drop=TRUE)
}

\arguments{

  \item{x}{A LoomExperiment object}

  \item{colGraphs,rowGraphs}{
    LoomGraphs to be placed in either the colGraphs or rowGraphs slot respectively
  }

  \item{value}{
    For \code{colGraphs<-} and \code{rowGraphs<-}, \code{value} will be the replacement to the slot.
  }

  \item{...}{
    For constructors, \code{...} will be passed on to the respective Experiment constructor.
  }

  \item{i,j}{
    Indices specifying elements to subset LoomGraph by
  }

  \item{drop}{
    For matrices and arrays.  If 'TRUE' the result is coerced to the lowest 
    posible dimesnion.  This only woeks for extracting elements, not for the 
    replacement.
  }

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

\author{
Daniel Van Twisk
}

