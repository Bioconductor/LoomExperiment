
setClassUnion("GenomicRanges_OR_GRangesList_OR_NULL", c("GenomicRanges", "GRangesList", "NULL"))

### =========================================================================
### LoomExperiment objects
### -------------------------------------------------------------------------
###

.LoomExperiment <- setClass("LoomExperiment",
    contains="SummarizedExperiment",
    representation(
        int_elementMetadata="DataFrame",
        int_colData="DataFrame",
        int_metadata="list",
        reducedDims="SimpleList",
        rowRanges="GenomicRanges_OR_GRangesList_OR_NULL",
        colGraphs="LoomGraphs",
        rowGraphs="LoomGraphs"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.LoomExperiment <- function(x)
{
    NULL
}

setValidity2("LoomExperiment", .valid.LoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_LoomExperiment <-
    function(se, colGraphs, rowGraphs)
{
    le <- .LoomExperiment(se, colGraphs=colGraphs, rowGraphs=rowGraphs)
    le <- BiocGenerics:::replaceSlots(le, rowRanges=rowRanges(se), check=FALSE)
    if (is(se, "SingleCellExperiment")) {
        le <- BiocGenerics:::replaceSlots(le, reducedDims=se@reducedDims, check=FALSE)
        le <- BiocGenerics:::replaceSlots(le, int_colData=se@int_colData, check=FALSE)
        le <- BiocGenerics:::replaceSlots(le, int_metadata=se@int_metadata, check=FALSE)
        le <- BiocGenerics:::replaceSlots(le, int_elementMetadata=se@int_elementMetadata, check=FALSE)
    }
    le
}

#'
#' The LoomExperiment representation class
#'
#' SummarizedExperiment-like class that facilitate the transition of 
#' SummarizedExperiment, RangedSummarizedExperiment, and
#' SingleCellExperiment objects to storage or retrieval in the loom format.
#'
#' @author Daniel Van Twisk
#'
#' @slot int_elementMetaData DataFrame. Mirrors int_elementMetadata from
#'  SingleCellExperiment.
#' @slot int_colData DataFrame. Mirrors int_colData from SingleCellExperiment.
#' @slot int_metadata List. Mirrors int_metadata from SingleCellExperiment.
#' @slot reducedDims SimpleList of matrices. Mirrors reducedDims from
#'  SingleCellExperiment.
#' @slot colGraphs LoomGraphs containing the loom file's colGraphs
#'  information. Optional.
#' @slot rowGraphs LoomGraphs containing the loom file's rowGraphs
#'  information. Optional.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment RangedSummarizedExperiment
#' 
#' @examples
#'
#'  library(SingleCellExperiment)
#'  counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
#'  sce <- SingleCellExperiment(assays = list(counts = counts))
#'  le <- LoomExperiment(sce)
#'  le
#' 
#' @rdname LoomExperiment-class
#' @export
LoomExperiment <- function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)[[1]]
    if (is(te, "SummarizedExperiment"))
        .new_LoomExperiment(te, colGraphs=colGraphs, rowGraphs=rowGraphs)
    else {
        se <- SummarizedExperiment(...)
        if(is(se, "RangedSummarizedExperiment"))
            se <- as(se, "SummarizedExperiment")
        .new_LoomExperiment(se, colGraphs=colGraphs, rowGraphs=rowGraphs)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_SummarizedExperiment_to_LoomExperiment <- function(from)
{
    .new_LoomExperiment(from,
                        rowGraphs=LoomGraphs(),
                        colGraphs=LoomGraphs())
}

#' @export
setAs("SummarizedExperiment", "LoomExperiment",
    .from_SummarizedExperiment_to_LoomExperiment
)

.from_RangedSummarizedExperiment_to_LoomExperiment <- function(from)
{
    .new_LoomExperiment(from,
                        rowGraphs=LoomGraphs(),
                        colGraphs=LoomGraphs())
}

#' @export
setAs("RangedSummarizedExperiment", "LoomExperiment",
    .from_RangedSummarizedExperiment_to_LoomExperiment
)

.from_LoomExperiment_to_RangedSummarizedExperiment <- function(from)
{
    SummarizedExperiment(assays=assays(from),
                         rowRanges=le@rowRanges,
                         colData=colData(le),
                         metadata=metadata(le))
}

#' @export
setAs("LoomExperiment", "RangedSummarizedExperiment",
    .from_LoomExperiment_to_RangedSummarizedExperiment
)

.from_LoomExperiment_to_SingleCellExperiment <- function(from)
{
    le <- SingleCellExperiment(assays=assays(from),
                               rowRanges=le@rowRanges,
                               colData=colData(le),
                               metadata=metadata(le),
                               reducedDims=le@reducedDims)
    le <- BiocGenerics:::replaceSlots(le, int_colData=se@int_colData, check=FALSE)
    le <- BiocGenerics:::replaceSlots(le, int_metadata=se@int_metadata, check=FALSE)
    le <- BiocGenerics:::replaceSlots(le, int_elementMetadata=se@int_elementMetadata, check=FALSE)
    le
}

#' @export
setAs("LoomExperiment", "SingleCellExperiment",
    .from_LoomExperiment_to_SingleCellExperiment
)

.from_SingleCellExperiment_to_LoomExperiment <- function(from)
{
    .new_LoomExperiment(from,
                        rowGraphs=LoomGraphs(),
                        colGraphs=LoomGraphs())
}

#' @export
setAs("SingleCellExperiment", "LoomExperiment",
    .from_SingleCellExperiment_to_LoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

.get.colGraphs <- function(x, ...)
{
    x@colGraphs
}

.replace.colGraphs <- function(x, ..., value)
{
    BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
}

.get.rowGraphs <- function(x, ...)
{
    x@rowGraphs
}

.replace.rowGraphs <- function(x, ..., value)
{
    BiocGenerics:::replaceSlots(x, rowGraphs=value, check=FALSE)
}

#'
#' @examples
#'  ## colGraphs of LoomExperiment object
#'  colGraphs(le)
#' 
#' @rdname LoomExperiment-class
#' @exportMethod colGraphs
setMethod("colGraphs", "LoomExperiment", .get.colGraphs)

#'
#' @examples
#'  ## replace colGraphs of LoomExperiment object
#'  lg <- LoomGraph(a=c(1, 2, 3), b=(3, 2, 1))
#'  colGraphs(le) <- lg
#' 
#' @rdname LoomExperiment-class
#' @exportMethod colGraphs<-
setReplaceMethod("colGraphs", "LoomExperiment", .replace.colGraphs)

#'
#' @examples
#'  ## rowGraphs of LoomExperiment object
#'  rowGraphs(le)
#' 
#' @rdname LoomExperiment-class
#' @exportMethod rowGraphs
setMethod("rowGraphs", "LoomExperiment", .get.rowGraphs)

#'
#' @examples
#'  ## replace rowGraphs of LoomExperiment object
#'  lg <- LoomGraph(a=c(1, 2, 3), b=(3, 2, 1))
#'  rowGraphs(le) <- lg
#' 
#' @rdname LoomExperiment-class
#' @exportMethod rowGraphs<-
setReplaceMethod("rowGraphs", "LoomExperiment", .replace.rowGraphs)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

.subset.LoomExperiment <- function(x, i, j, ...)
{
    if (!missing(i))
        rowGraphs(x) <- endoapply(rowGraphs(x), function(y) y[i,])
    if (!missing(j))
        colGraphs(x) <- endoapply(colGraphs(x), function(y) y[i,])
    callNextMethod()
}

#'
#'
#' @rdname LoomExperiment-class
#' @exportMethod [
setMethod("[", c("LoomExperiment", "ANY", "ANY"), .subset.LoomExperiment)

.show.LoomExperiment <- function(object)
{
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    callNextMethod()
    
    ## SingleCellExperiment Terms
    if(length(object@reducedDims))
        scat("reducedDimNames(%d): %s\n", names(object@reducedDims))
    else
        cat("reducedDims(0): NULL\n")

    ## LoomExperiment Terms
    if (length(object@rowGraphs) > 0)
        scat("rowGraphs(%d): %s\n", names(object@rowGraphs))
    else
        cat("rowGraphs(0): NULL\n")
    if (length(object@colGraphs) > 0)
        scat("colGraphs(%d): %s\n", names(object@colGraphs))
    else
        cat("colGraphs(0): NULL\n")
}

#' @export
setMethod("show", "LoomExperiment", .show.LoomExperiment)

