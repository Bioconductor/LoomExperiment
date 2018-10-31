
### =========================================================================
### RangedLoomExperiment objects
### -------------------------------------------------------------------------
###

#' @import SummarizedExperiment
#' @export
setClass('RangedLoomExperiment',
    contains=c('RangedSummarizedExperiment', 'LoomExperiment')
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_RangedLoomExperiment <- function(se, colGraphs, rowGraphs)
{
    new('RangedLoomExperiment', se, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

#' @export
RangedLoomExperiment <-
    function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)
    if(length(te) > 0 && is(te[[1]], 'SummarizedExperiment'))
        rse <- te[[1]]
    else {
        rse <- SummarizedExperiment(...)
        if(!is(rse, 'RangedSummarizedExperiment'))
            rse <- as(rse, 'RangedSummarizedExperiment')
    }
    .new_RangedLoomExperiment(rse,
                              colGraphs=.change.nnode(colGraphs, ncol(rse)),
                              rowGraphs=.change.nnode(rowGraphs, nrow(rse)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_RangedSummarizedExperiment_to_RangedLoomExperiment <- function(from)
{
    .new_RangedLoomExperiment(from,
                              colGraphs=LoomGraphs(),
                              rowGraphs=LoomGraphs())
}

setAs('RangedSummarizedExperiment', 'RangedLoomExperiment',
    .from_RangedSummarizedExperiment_to_RangedLoomExperiment
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Methods.
###

#setMethod('[', 'RangedLoomExperiment', .subset.LoomExperiment)

