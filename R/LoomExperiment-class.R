
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
    })
    row_test <- lapply(rlgs, function(lg) {
        txt <- 'All LoomGraph objects in LoomExperiment must reference a row in the LoomExperiment'
        test <- from(lg) %in% rows & to(lg) %in% rows
        if ((length(test) == 0 || !all(test)) && length(from(lg)) > 0)
            return(txt)
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
    if (length(te) > 0 && is(te[[1]], 'SummarizedExperiment'))
        se <- te[[1]]
    else {
        se <- SummarizedExperiment(...)
        if(is(se, 'RangedSummarizedExperiment'))
            se <- as(se, 'SummarizedExperiment')
    }
    .new_LoomExperiment(se, colGraphs=colGraphs, rowGraphs=rowGraphs)
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
    x@colGraphs <- value
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
    x@rowGraphs <- value
    validObject(x)
    x
}

#' @export
setReplaceMethod('rowGraphs', 'LoomExperiment', .replace.rowGraphs)

#' @export
setMethod('[', c('LoomExperiment', 'ANY', 'ANY'), .subset.LoomExperiment)

.selectHits.LoomExperiment <- function(x, i, ...)
{
    rowGraphs(x) <- endoapply(rowGraphs(x), function(y) selectHits(y, i))
    colGraphs(x) <- endoapply(colGraphs(x), function(y) selectHits(y, i))
    x
}

#' @export
setMethod('selectHits', c('LoomExperiment', 'ANY'), .selectHits.LoomExperiment)

.dropHits.LoomExperiment <- function(x, i, ...)
{ 
    rowGraphs(x) <- endoapply(rowGraphs(x), function(y) dropHits(y, i))
    colGraphs(x) <- endoapply(colGraphs(x), function(y) dropHits(y, i))
    x
}

#' @export
setMethod('dropHits', c('LoomExperiment', 'ANY'), .dropHits.LoomExperiment)

.dropHits.replace.LoomExperiment <- function(x, i, ..., value)
{
    rowGraphs(x) <- endoapply(rowGraphs(x), function(y){
        dropHits(y, i) <- value
        y
    })
    colGraphs(x) <- endoapply(colGraphs(x), function(y){
        dropHits(y, i) <- value
        y
    })
    x
}

#' @export
setReplaceMethod('dropHits', c('LoomExperiment', 'ANY', 'ANY'),
    .dropHits.replace.LoomExperiment)

#' @export
setMethod('show', 'LoomExperiment', .show.LoomExperiment)

