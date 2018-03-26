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

.new_SingleCellLoomExperiment <- function(assays, names, rowData, colData,
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
    new("SingleCellLoomExperiment", NAMES=names,
                                elementMetadata=rowData,
                                colData=colData,
                                assays=assays,
                                colGraphs=colGraphs,
                                rowGraphs=rowGraphs,
                                metadata=as.list(metadata))
}

#' @export
setGeneric("SingleCellLoomExperiment",
    function(assays, ...) standardGeneric("SingleCellLoomExperiment")
)

#' @export
setMethod("SingleCellLoomExperiment", "SimpleList",
   function(assays, rowData=NULL, rowRanges=GRangesList(), colData=DataFrame(),
            colGraphs=SimpleList(), rowGraphs=SimpleList(), metadata=list())
{
    if (missing(colData) && 0L != length(assays)) {
        assay <- assays[[1]]
        nms <- colnames(assay)
        colData <- DataFrame(x=seq_len(ncol(assay)), row.names=nms)[, FALSE]
    } else if (!missing(colData)) {
        colData <- as(colData, "DataFrame")
        if (is.null(rownames(colData)))
            rownames(colData) <- SummarizedExperiment:::.get_colnames_from_assays(assays)
    }
    ans_colnames <- rownames(colData)

    if (is.null(rowData)) {
        if (missing(rowRanges)) {
            ans_rownames <- SummarizedExperiment:::.get_rownames_from_assays(assays)
        } else {
            if (is.null(names(rowRanges)))
                names(rowRanges) <- SummarizedExperiment:::.get_rownames_from_assays(assays)
            ans_rownames <- names(rowRanges)
        }
    } else {
        if (!missing(rowRanges))
            stop("only one of 'rowData' and 'rowRanges' can be specified")
        if (is(rowData, "GenomicRanges_OR_GRangesList")) {
            rowRanges <- rowData
            if (is.null(names(rowRanges)))
                names(rowRanges) <- SummarizedExperiment:::.get_rownames_from_assays(assays)
            ans_rownames <- names(rowRanges)
        } else {
            rowData <- as(rowData, "DataFrame")
            ans_rownames <- rownames(rowData)
            if (is.null(ans_rownames))
                ans_rownames <- SummarizedExperiment:::.get_rownames_from_assays(assays)
        }
    }

    ## validate
#    ok <- vapply(assays, function(x) {
#        colnames <- colnames(x)
#        test <- is.null(colnames) || identical(colnames, ans_colnames)
#        if (!test)
#            stop("assay colnames() must be NULL or equal colData rownames()")
#
#        rownames <- rownames(x)
#        test <- test &&
#            is.null(rownames) || identical(rownames, ans_rownames)
#        if (!test) {
#            txt <- "assay rownames() must be NULL or equal rowData rownames() /
#                    rowRanges names()"
#            stop(paste(strwrap(txt, exdent=2), collapse="\n"))
#        }
#
#        test
#    }, logical(1))

    assays <- Assays(assays)

    if (missing(rowRanges) && !is(rowData, "GenomicRanges_OR_GRangesList")) {
        .new_SingleCellLoomExperiment(assays, ans_rownames, rowData, colData,
                            colGraphs, rowGraphs, metadata)
    } else {
        .new_RangedSummarizedExperiment(assays, rowRanges, colData, metadata)
    }
})

#' @export
setMethod("SingleCellLoomExperiment", "ANY",
    function(assays, ...)
{
    if (is.matrix(assays) && is.list(assays))
        assays <- list(assays)
    SingleCellLoomExperiment(assays, ...)
})

#' @export
setMethod("SingleCellLoomExperiment", "list",
    function(assays, ...)
{
    SingleCellLoomExperiment(do.call(SimpleList, assays), ...)
})

#' @export
setMethod("SingleCellLoomExperiment", "missing",
    function(assays, ...)
{
    SingleCellLoomExperiment(SimpleList(), ...)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_SingleCellLoomExperiment_to_SummarizedExperiment <- function(from)
{
    SummarizedExperiment(assays=from@assays,
                         rowData=NULL,
                         colData=from@colData,
                         metadata=from@elementMetadata)
}

setAs("SingleCellLoomExperiment", "SummarizedExperiment",
    .from_SingleCellLoomExperiment_to_SummarizedExperiment
)

.from_SummarizedExperiment_to_SingleCellLoomExperiment <- function(from)
{
    SingleCellLoomExperiment(assays=from@assays,
                         rowData=NULL,
                         colData=from@colData,
                         metadata=from@elementMetadata)
}

setAs("SummarizedExperiment", "SingleCellLoomExperiment",
    .from_SummarizedExperiment_to_SingleCellLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

setGeneric("colGraphs", function(x, ...) standardGeneric("colGraphs"))

setMethod("colGraphs", "SingleCellLoomExperiment",
    function(x, ...) x@colGraphs)

setGeneric("colGraphs<-", function(x, ..., value) standardGeneric("colGraphs<-"))

setReplaceMethod("colGraphs", "SingleCellLoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
    }
)

setGeneric("rowGraphs", function(x, ...) standardGeneric("rowGraphs"))

setMethod("rowGraphs", "SingleCellLoomExperiment",
    function(x, ...) x@rowGraphs)

setGeneric("rowGraphs<-", function(x, ..., value) standardGeneric("rowGraphs<-"))

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
