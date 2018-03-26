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

.new_RangedLoomExperiment <- function(assays, names, rowData, colData,
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
    new("RangedLoomExperiment", NAMES=names,
                                elementMetadata=rowData,
                                colData=colData,
                                assays=assays,
                                colGraphs=colGraphs,
                                rowGraphs=rowGraphs,
                                metadata=as.list(metadata))
}

#' @export
setGeneric("RangedLoomExperiment",
    function(assays, ...) standardGeneric("RangedLoomExperiment")
)

#' @export
setMethod("RangedLoomExperiment", "SimpleList",
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
        .new_RangedLoomExperiment(assays, ans_rownames, rowData, colData,
                            colGraphs, rowGraphs, metadata)
    } else {
        .new_RangedSummarizedExperiment(assays, rowRanges, colData, metadata)
    }
})

#' @export
setMethod("RangedLoomExperiment", "ANY",
    function(assays, ...)
{
    if (is.matrix(assays) && is.list(assays))
        assays <- list(assays)
    RangedLoomExperiment(assays, ...)
})

#' @export
setMethod("RangedLoomExperiment", "list",
    function(assays, ...)
{
    RangedLoomExperiment(do.call(SimpleList, assays), ...)
})

#' @export
setMethod("RangedLoomExperiment", "missing",
    function(assays, ...)
{
    RangedLoomExperiment(SimpleList(), ...)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_RangedLoomExperiment_to_SummarizedExperiment <- function(from)
{
    SummarizedExperiment(assays=from@assays,
                         rowData=NULL,
                         colData=from@colData,
                         metadata=from@elementMetadata)
}

setAs("RangedLoomExperiment", "SummarizedExperiment",
    .from_RangedLoomExperiment_to_SummarizedExperiment
)

.from_SummarizedExperiment_to_RangedLoomExperiment <- function(from)
{
    RangedLoomExperiment(assays=from@assays,
                         rowData=NULL,
                         colData=from@colData,
                         metadata=from@elementMetadata)
}

setAs("SummarizedExperiment", "RangedLoomExperiment",
    .from_SummarizedExperiment_to_RangedLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

setGeneric("colGraphs", function(x, ...) standardGeneric("colGraphs"))

setMethod("colGraphs", "RangedLoomExperiment",
    function(x, ...) x@colGraphs)

setGeneric("colGraphs<-", function(x, ..., value) standardGeneric("colGraphs<-"))

setReplaceMethod("colGraphs", "RangedLoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
    }
)

setGeneric("rowGraphs", function(x, ...) standardGeneric("rowGraphs"))

setMethod("rowGraphs", "RangedLoomExperiment",
    function(x, ...) x@rowGraphs)

setGeneric("rowGraphs<-", function(x, ..., value) standardGeneric("rowGraphs<-"))

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
