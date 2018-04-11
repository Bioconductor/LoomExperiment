setClassUnion("GenomicRanges_OR_GRangesList_OR_NULL", c("GenomicRanges", "GRangesList", "NULL"))

### =========================================================================
### LoomExperiment objects
### -------------------------------------------------------------------------
###

#' LoomExperiment
#'
#' A class that helps facilitate the transition of SummarizedExperiment objects
#' to .loom files and vise versa.
#'
#' @slot colGraphs A SimpleList containing the colGraphs information
#' @slot rowGraphs A SimpleList containing the rowGraphs information
#'
#' @author Daniel Van Twisk
#' @import SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export

setClass("LoomExperiment",
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
    le <- new("LoomExperiment", se, colGraphs=colGraphs, rowGraphs=rowGraphs)
    le <- BiocGenerics:::replaceSlots(le, rowRanges=rowRanges(se), check=FALSE)
    if (is(se, "SingleCellExperiment")) {
        le <- BiocGenerics:::replaceSlots(le, reducedDims=se@reducedDims, check=FALSE)
        le <- BiocGenerics:::replaceSlots(le, int_colData=se@int_colData, check=FALSE)
        le <- BiocGenerics:::replaceSlots(le, int_metadata=se@int_metadata, check=FALSE)
        le <- BiocGenerics:::replaceSlots(le, int_elementMetadata=se@int_elementMetadata, check=FALSE)
    }
    le
}

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

setAs("SummarizedExperiment", "LoomExperiment",
    .from_SummarizedExperiment_to_LoomExperiment
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

setMethod("colGraphs", "LoomExperiment", .get.colGraphs)

setReplaceMethod("colGraphs", "LoomExperiment", .replace.colGraphs)

setMethod("rowGraphs", "LoomExperiment", .get.rowGraphs)

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

setMethod("show", "LoomExperiment", .show.LoomExperiment)
