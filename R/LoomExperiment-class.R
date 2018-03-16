### =========================================================================
### LoomExperiment objects
### -------------------------------------------------------------------------
###

#' LoomExperiment
#'
#' A class that helps facilitate the transition of SummarizedExperiment objects
#' to .loom files and vise versa.
#'
#' @slot colGraph A SimpleList containing the colGraph information
#' @slot rowGraph A SimpleList containing the rowGraph information
#'
#' @author Daniel Van Twisk
#' @import SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export

setClass("LoomExperiment",
    contains="SummarizedExperiment",
    representation(
        colGraph="SimpleList",
        rowGraph="SimpleList"
    ),
    prototype(
        colGraph=SimpleList(),
        rowGraph=SimpleList()
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.LoomExperiment.checkGraph <- function(graph, type)
{
    if (length(graph) == 0L)
        return(NULL)
    test_type <- vapply(graph, function(x) is(x, "DataFrame"), logical(1))
    if(!all(test_type)) {
        txt <- sprintf("\n '%s' must only contain DataFrames", type)
    }
    len <- lengths(graph)
    if(length(len)==2) {
        if (all(names(len) == c('a', 'b'))) {
            txt <- sprintf("\n '%s' vectors must be named 'a' and 'b'", type)
            return(txt)
        }
    }
    if(length(len)==3) {
        if (all(names(len) == c('a', 'b', 'w'))) {
            txt <- sprintf(
                "\n '%s' vectors must be named 'a', 'b', and 'w'", type
            )
            return(txt)
        }
    }
    numerical <- vapply(graph, is.numeric, logical(1))
    if (all(numerical)) {
        txt <- sprintf(
            "\n '%s' vector must be numeric", type
        )
        return(txt)
    }
    NULL
}

.valid.LoomExperiment.Graphs <- function(x)
{
    .valid.LoomExperiment.checkGraph(x@colGraph, "colGraph")
    .valid.LoomExperiment.checkGraph(x@rowGraph, "rowGraph")
}

.valid.LoomExperiment <- function(x)
{
    .valid.LoomExperiment.Graphs(x)
}

setValidity2("LoomExperiment", .valid.LoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_LoomExperiment <- function(assays, names, rowData, colData,
                                colGraph, rowGraph, metadata)
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
                                colGraph=colGraph,
                                rowGraph=rowGraph,
                                metadata=as.list(metadata))
}

#' @export
setGeneric("LoomExperiment",
    function(assays, ...) standardGeneric("LoomExperiment")
)

#' @export
setMethod("LoomExperiment", "SimpleList",
   function(assays, rowData=NULL, rowRanges=GRangesList(), colData=DataFrame(),
            colGraph=SimpleList(), rowGraph=SimpleList(), metadata=list())
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
        .new_LoomExperiment(assays, ans_rownames, rowData, colData,
                            colGraph, rowGraph, metadata)
    } else {
        .new_RangedSummarizedExperiment(assays, rowRanges, colData, metadata)
    }
})

#' @export
setMethod("LoomExperiment", "ANY",
    function(assays, ...)
{
    if (is.matrix(assays) && is.list(assays))
        assays <- list(assays)
    LoomExperiment(assays, ...)
})

#' @export
setMethod("LoomExperiment", "list",
    function(assays, ...)
{
    LoomExperiment(do.call(SimpleList, assays), ...)
})

#' @export
setMethod("LoomExperiment", "missing",
    function(assays, ...)
{
    LoomExperiment(SimpleList(), ...)
})


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

setGeneric("colGraph", function(x, ...) standardGeneric("colGraph"))

setMethod("colGraph", "LoomExperiment",
    function(x, ...) x@colGraph)

setGeneric("colGraph<-", function(x, ..., value) standardGeneric("colGraph<-"))

setReplaceMethod("colGraph", "LoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, colGraph=value, check=FALSE)
    }
)

setGeneric("rowGraph", function(x, ...) standardGeneric("rowGraph"))

setMethod("rowGraph", "LoomExperiment",
    function(x, ...) x@rowGraph)

setGeneric("rowGraph<-", function(x, ..., value) standardGeneric("rowGraph<-"))

setReplaceMethod("rowGraph", "LoomExperiment",
    function(x, ..., value) {
        BiocGenerics:::replaceSlots(x, rowGraph=value, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

setMethod("[", c("LoomExperiment", "ANY", "missing"),
    function(x, i, j, ...)
{
    res <- lapply(colGraph(x), function(y) {
        y[!(!y[,"a"] %in% i | !y[,"b"] %in% i),]
    })
    res <- as(res, "SimpleList")
    LoomExperiment(assays=assays(le), rowData=rowData(le), colData=x@colData,
                   metadata=metadata(le), colGraph=res, rowGraph=rowGraph(le))
})

setReplaceMethod("[", c("LoomExperiment", "ANY", "missing"),
    function(x, i, j, ...)
{
    res <- lapply(colGraph(x), function(y) {
        y[!(!y[,"a"] %in% i | !y[,"b"] %in% i),]
    })
    colGraph(x) <- as(res, "SimpleList")
    #LoomExperiment(assays=x@assays, colData=x@colData, metadata=x@elementMetaData,
    #               colGraph=x@colGraph, rowGraph=x@rowGraph)
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
    if (length(object@rowGraph) > 0)
        scat("rowGraph(%d): %s\n", names(object@rowGraph))
    else
        cat("rowGraph(0): NULL\n")
    if (length(object@colGraph) > 0)
        scat("colGraph(%d): %s\n", names(object@colGraph))
    else
        cat("colGraph(0): NULL\n")
})
