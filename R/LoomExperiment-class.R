
### =========================================================================
### LoomExperiment objects
### -------------------------------------------------------------------------
###

#' @import SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
setClass('LoomExperiment',
    contains=c('SummarizedExperiment'),
    representation(
        colGraphs='LoomGraphs',
        rowGraphs='LoomGraphs'
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.Experiment <- function(x)
{
    clgs <- colGraphs(x)
    rlgs <- rowGraphs(x)
    cols <- c(0, seq_len(dim(x)[[2]]))
    rows <- c(0, seq_len(dim(x)[[1]]))
    ## Check that no 'a' or 'b' columns in LoomGraphs lie outside of dimensions.
    col_test <- lapply(clgs, function(lg) {
        txt <- 'All LoomGraph objects in LoomExperiment must reference a column in the LoomExperiment'
        test <- from(lg) %in% cols & to(lg) %in% cols
        if ((length(test) == 0 || !all(test)) && length(from(lg)) > 0)
            return(txt)
#        txt <- 'All LoomGraph objects in LoomExperiment colGraphs must have the same nnode as columns in the LoomExperiment'
#        if (nnode(lg) != nrow(x))
#            return(txt)
    })
    row_test <- lapply(rlgs, function(lg) {
        txt <- 'All LoomGraph objects in LoomExperiment must reference a row in the LoomExperiment'
        test <- from(lg) %in% rows & to(lg) %in% rows
        if ((length(test) == 0 || !all(test)) && length(from(lg)) > 0)
            return(txt)
#        txt <- 'All LoomGraph objects in LoomExperiment rowGraphs must have the same nnode as rows in the LoomExperiment'
#        if (nnode(lg) != nrow(x))
#            return(txt)
    })
    res <- list(col_test, row_test)
    unlist(res)
}

setValidity2('LoomExperiment', .valid.Experiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_LoomExperiment <-
    function(se, colGraphs, rowGraphs)
{
    new('LoomExperiment', se, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

#' @export
#' @importFrom methods as
LoomExperiment <- function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)
    if (length(te) > 0 && is(te[[1]], 'SingleCellExperiment')) {
        se <- te[[1]]
        .new_SingleCellLoomExperiment(se,
                            colGraphs=.change.nnode(colGraphs, ncol(se)),
                            rowGraphs=.change.nnode(rowGraphs, nrow(se)))
    }
    else if (length(te) > 0 && is(te[[1]], 'RangedSummarizedExperiment')) {
        se <- te[[1]]
        .new_RangedLoomExperiment(se,
                            colGraphs=.change.nnode(colGraphs, ncol(se)),
                            rowGraphs=.change.nnode(rowGraphs, nrow(se)))
    }
    else if (length(te) > 0 && is(te[[1]], 'SummarizedExperiment')) {
        se <- te[[1]]
        .new_LoomExperiment(se,
                            colGraphs=.change.nnode(colGraphs, ncol(se)),
                            rowGraphs=.change.nnode(rowGraphs, nrow(se)))
    }
    else {
        se <- SummarizedExperiment(...)
        if(is(se, 'RangedSummarizedExperiment'))
            se <- as(se, 'SummarizedExperiment')
        .new_LoomExperiment(se,
                            colGraphs=.change.nnode(colGraphs, ncol(se)),
                            rowGraphs=.change.nnode(rowGraphs, nrow(se)))
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

setAs('SummarizedExperiment', 'LoomExperiment',
    .from_SummarizedExperiment_to_LoomExperiment
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

.change.nnode <- function(x, nr)
{
    endoapply(x, function(y) {
        y@nLnode <- nr
        y@nRnode <- nr
        y
    })
}

.get.colGraphs <- function(x, ...)
{
    x@colGraphs
}

#' @export
setMethod('colGraphs', 'LoomExperiment', .get.colGraphs)

#' @importFrom methods validObject callNextMethod
.replace.colGraphs <- function(x, ..., value)
{
    #x <- BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
    x@colGraphs <- .change.nnode(value, ncol(x))
    validObject(x)
    x
}

#' @export
setReplaceMethod('colGraphs', 'LoomExperiment', .replace.colGraphs)

.get.rowGraphs <- function(x, ...)
{
    x@rowGraphs
}

#' @export
setMethod('rowGraphs', 'LoomExperiment', .get.rowGraphs)

.replace.rowGraphs <- function(x, ..., value)
{
    #x <- BiocGenerics:::replaceSlots(x, rowGraphs=value, check=FALSE)
    x@rowGraphs <- .change.nnode(value, nrow(x))
    validObject(x)
    x
}

#' @export
setReplaceMethod('rowGraphs', 'LoomExperiment', .replace.rowGraphs)

#' @export
setMethod('[', c('LoomExperiment', 'ANY', 'ANY'), .subset.LoomExperiment)

.loomSelectHits.LoomExperiment <- function(x, i, ...)
{
    rowGraphs(x) <- endoapply(rowGraphs(x), function(y) loomSelectHits(y, i))
    colGraphs(x) <- endoapply(colGraphs(x), function(y) loomSelectHits(y, i))
    x
}

#' @importFrom S4Vectors rbind nnode
#' @export
setMethod('rbind', 'LoomExperiment', .rbind.LoomExperiment)

#' @importFrom S4Vectors cbind
#' @export
setMethod('cbind', 'LoomExperiment', .cbind.LoomExperiment)

.loomDropHits.LoomExperiment <- function(x, i, ...)
{ 
    rowGraphs(x) <- endoapply(rowGraphs(x), function(y) loomDropHits(y, i))
    colGraphs(x) <- endoapply(colGraphs(x), function(y) loomDropHits(y, i))
    x
}

.loomDropHits.replace.LoomExperiment <- function(x, i, ..., value)
{
    rowGraphs(x) <- endoapply(rowGraphs(x), function(y){
        loomDropHits(y, i) <- value
        y
    })
    colGraphs(x) <- endoapply(colGraphs(x), function(y){
        loomDropHits(y, i) <- value
        y
    })
    x
}

#' @export
setMethod('show', 'LoomExperiment', .show.LoomExperiment)

