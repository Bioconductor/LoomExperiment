
\name{LoomFile-class}
\docType{class}

%% Classes:
\alias{class:LoomFile}
\alias{LoomFile-class}
\alias{class:CompressedFile}
\alias{CompressedFile-class}
\alias{LoomFileList-class}

%% Accessors
\alias{path}
\alias{path,LoomFile-method}
\alias{resource}

%% Show:
\alias{show,LoomFile-method}

%% Related functions:
\alias{FileForFormat}

%% Coerce:
\alias{as.character,LoomFile-method}

\title{LoomFile objects}

\description{
  A \code{LoomFile} object is the base class for classes representing
  files accessible with \code{rtracklayer}. It wraps a resource (either a path,
  URL or connection).
}

\section{Accessor Methods}{
  In the code snippets below, \code{x} represents a \code{LoomFile}
  object.
  
  \describe{
    \item{}{
      \code{path(x)}:
      Gets the path, as a \code{character} vector, to the resource
      represented by the \code{LoomFile} object, if possible.
    }
    \item{}{
      \code{resource(x)}:
      Gets the low-level resource, either a character vector (a path or
      URL) or a connection.
    }
  }
}

\section{Coercion}{
  \describe{
    \item{}{
      \code{as.character(x)}:
      Returns the path of the file as a character vector.
    }
  }
}

\section{Related functions}{
  \describe{
    \item{}{
      \code{FileForFormat(path, format = file_ext(path))}:
      Determines the file type of \code{path} and returns
      a high-level file object such as BamFile, BEDFile,
      BigWigFile etc.. 
    }
  }
}

\author{Daniel Van Twisk}

\keyword{methods}
\keyword{classes}

