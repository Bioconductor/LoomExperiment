### =========================================================================
### RangedLoomExperiment objects
### -------------------------------------------------------------------------
###

#' RangedLoomExperiment
#'
#' A class that helps facilitate the transition of SummarizedExperiment objects
#' to .loom files and vise versa.
#'
#' @slot colGraphs A LoomGraphs object containing the colGraphs information
#' @slot rowGraphs A LoomGraphs containing the rowGraphs information
#'
#' @author Daniel Van Twisk
#' @import SummarizedExperiment
#' @importFrom SummarizedExperiment RangedSummarizedExperiment
#' @export

setClass("RangedLoomExperiment",
    contains="RangedSummarizedExperiment",
    representation(
        colGraphs="LoomGraphs",
        rowGraphs="LoomGraphs"
    ),
    prototype(
        colGraphs=NULL,
        rowGraphs=NULL
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.RangedLoomExperiment.Graphs <- function(x)
{
    NULL
}

.valid.RangedLoomExperiment <- function(x)
{
    .valid.RangedLoomExperiment.Graphs(x)
}

#setValidity2("RangedLoomExperiment", .valid.RangedLoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_RangedLoomExperiment <- function(assays, rowRanges, colData,
                                colGraphs, rowGraphs, metadata)
{
    elementMetadata <- S4Vectors:::make_zero_col_DataFrame(length(rowRanges))
    if (!is(assays, "Assays"))
        assays <- Assays(assays)
    new("RangedLoomExperiment", rowRanges=rowRanges,
                                elementMetadata=rowData,
                                colData=colData,
                                assays=assays,
                                colGraphs=colGraphs,
                                rowGraphs=rowGraphs,
                                metadata=as.list(metadata))
}

#' @export
LoomExperiment <- function(..., colGraphs=NULL, rowGraphs=NULL) {
    se <- SummarizedExperiment(...)
    if (!is(se, "RangedSummarizedExperiment")) {
        .new_LoomExperiment(assays=se@assays,
                            names=se@NAMES,
                            rowData=NULL,
                            colData=se@colData,
                            metadata=se@metadata,
                            colGraphs=colGraphs,
                            rowGraphs=rowGraphs)
    } else {
        .new_RangedLoomExperiment(assays=se@assays,
                                  rowRanges=se@rowRanges,
                                  colData=se@colData,
                                  metadata=se@metadata,
                                  colGraphs=colGraphs,
                                  rowGraphs=rowGraphs)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_RangedLoomExperiment_to_RangedSummarizedExperiment <- function(from)
{
    SummarizedExperiment:::.new_RangedSummarizedExperiment(from@assays,
                                   from@rowRanges,
                                   from@colData,
                                   from@elementMetadata)
}

setAs("RangedLoomExperiment", "RangedSummarizedExperiment",
    .from_RangedLoomExperiment_to_RangedSummarizedExperiment
)

.from_RangedExperiment_to_RangedLoomExperiment <- function(from)
{
    .new_RangedLoomExperiment(assays=from@assays,
                              rowRanges=from@rowRanges,
                              colData=from@colData,
                              metadata=from@elementMetadata,
                              colGraphs=NULL,
                              rowGraphs=NULL)
}

setAs("RangedSummarizedExperiment", "RangedLoomExperiment",
    .from_RangedSummarizedExperiment_to_RangedLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

setMethod("colGraphs", "RangedLoomExperiment",
    function(x, ...) x@colGraphs)

setReplaceMethod("colGraphs", "RangedLoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
    }
)

setMethod("rowGraphs", "RangedLoomExperiment",
    function(x, ...) x@rowGraphs)

setReplaceMethod("rowGraphs", "RangedLoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, rowGraphs=value, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

setMethod("[", c("RangedLoomExperiment", "ANY", "ANY"),
    function(x, i, j, ...)
{
    if (!missing(i))
        rowGraphs(x) <- endoapply(rowGraphs(x), function(y) y[i,])
    if (!missing(j))
        colGraphs(x) <- endoapply(colGraphs(x), function(y) y[i,])
    callNextMethod()
})

setMethod("show", "RangedLoomExperiment",
    function(object)
{
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    callNextMethod()
    if (length(object@rowGraphs) > 0)
        scat("rowGraphs(%d): %s\n", names(object@rowGraphs))
    else
        cat("rowGraphs(0): NULL\n")
    if (length(object@colGraphs) > 0)
        scat("colGraphs(%d): %s\n", names(object@colGraphs))
    else
        cat("colGraphs(0): NULL\n")
})
