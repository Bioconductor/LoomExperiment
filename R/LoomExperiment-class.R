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

.valid.LoomExperiment.Graphs <- function(x)
{
    NULL
}

.valid.LoomExperiment <- function(x)
{
    .valid.LoomExperiment.Graphs(x)
}

#setValidity2("LoomExperiment", .valid.LoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_LoomExperiment <- function(assays, names, rowData, colData,
                                colGraphs, rowGraphs, metadata)
{
    if (!is(assays, "Assays"))
        assays <- Assays(assays)
    if (is.null(rowData)) {
        if (is.null(names))
            nrow <- nrow(assays)
        else
            nrow <- length(names)
        rowData <- S4Vectors:::make_zero_col_DataFrame(nrow)
    } else {
        rownames(rowData) <- NULL
    }
    new("LoomExperiment", NAMES=names,
                                elementMetadata=rowData,
                                colData=colData,
                                assays=assays,
                                colGraphs=colGraphs,
                                rowGraphs=rowGraphs,
                                metadata=as.list(metadata))
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_LoomExperiment_to_SummarizedExperiment <- function(from)
{
    SummarizedExperiment(assays=from@assays,
                         rowData=NULL,
                         colData=from@colData,
                         metadata=from@elementMetadata)
}

setAs("LoomExperiment", "SummarizedExperiment",
    .from_LoomExperiment_to_SummarizedExperiment
)

.from_SummarizedExperiment_to_LoomExperiment <- function(from)
{
    LoomExperiment(assays=from@assays,
                         rowData=NULL,
                         colData=from@colData,
                         metadata=from@elementMetadata)
}

setAs("SummarizedExperiment", "LoomExperiment",
    .from_SummarizedExperiment_to_LoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

setGeneric("colGraphs", function(x, ...) standardGeneric("colGraphs"))

setMethod("colGraphs", "LoomExperiment",
    function(x, ...) x@colGraphs)

setGeneric("colGraphs<-", function(x, ..., value) standardGeneric("colGraphs<-"))

setReplaceMethod("colGraphs", "LoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
    }
)

setGeneric("rowGraphs", function(x, ...) standardGeneric("rowGraphs"))

setMethod("rowGraphs", "LoomExperiment",
    function(x, ...) x@rowGraphs)

setGeneric("rowGraphs<-", function(x, ..., value) standardGeneric("rowGraphs<-"))

setReplaceMethod("rowGraphs", "LoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, rowGraphs=value, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

setMethod("[", c("LoomExperiment", "ANY", "ANY"),
    function(x, i, j, ...)
{
    if (!missing(i))
        rowGraphs(x) <- endoapply(rowGraphs(x), function(y) y[i,])
    if (!missing(j))
        colGraphs(x) <- endoapply(colGraphs(x), function(y) y[i,])
    callNextMethod()
})

setMethod("show", "LoomExperiment",
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
