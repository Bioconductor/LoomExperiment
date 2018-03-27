### =========================================================================
### SingleCellLoomExperiment objects
### -------------------------------------------------------------------------
###

#' SingleCellLoomExperiment
#'
#' A class that helps facilitate the transition of SummarizedExperiment objects
#' to .loom files and vise versa.
#'
#' @slot colGraphs A SimpleList containing the colGraphs information
#' @slot rowGraphs A SimpleList containing the rowGraphs information
#'
#' @author Daniel Van Twisk
#' @import SingleCellExperiment
#' @importFrom SignleCellExperiment SingleCellExperiment
#' @export

setClass("SingleCellLoomExperiment",
    contains="SingleCellExperiment",
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

.valid.SingleCellLoomExperiment.Graphs <- function(x)
{
    NULL
}

.valid.SingleCellLoomExperiment <- function(x)
{
    .valid.SingleCellLoomExperiment.Graphs(x)
}

#setValidity2("SingleCellLoomExperiment", .valid.SingleCellLoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

SingleCellLoomExperiment <- function(..., colGraphs=NULL, rowGraphs=NULL) {
    sce <- SingleCellExperiment(...)
    new("SingleCellLoomExperiment", sce, colGraphs=colGraphs,
        rowGraphs=rowGraphs)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_SingleCellLoomExperiment_to_SingleCellExperiment <- function(from)
{
    SingleCellExperiment:::.SingleCellExperiment(assays=from@assays,
                         rowData=NULL,
                         colData=from@colData,
                         metadata=from@elementMetadata)
}

setAs("SingleCellLoomExperiment", "SingleCellExperiment",
    .from_SingleCellLoomExperiment_to_SingleCellExperiment
)

.from_SingleCellExperiment_to_SingleCellLoomExperiment <- function(from)
{
    SingleCellLoomExperiment(assays=from@assays,
                         rowData=NULL,
                         colData=from@colData,
                         metadata=from@elementMetadata)
}

setAs("SingleCellExperiment", "SingleCellLoomExperiment",
    .from_SingleCellExperiment_to_SingleCellLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

setMethod("colGraphs", "SingleCellLoomExperiment",
    function(x, ...) x@colGraphs)

setReplaceMethod("colGraphs", "SingleCellLoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
    }
)

setMethod("rowGraphs", "SingleCellLoomExperiment",
    function(x, ...) x@rowGraphs)

setReplaceMethod("rowGraphs", "SingleCellLoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, rowGraphs=value, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

setMethod("[", c("SingleCellLoomExperiment", "ANY", "ANY"),
    function(x, i, j, ...)
{
    if (!missing(i))
        rowGraphs(x) <- endoapply(rowGraphs(x), function(y) y[i,])
    if (!missing(j))
        colGraphs(x) <- endoapply(colGraphs(x), function(y) y[i,])
    callNextMethod()
})

setMethod("show", "SingleCellLoomExperiment",
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
