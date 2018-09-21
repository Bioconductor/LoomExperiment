
### =========================================================================
### SummarizedLoomExperiment objects
### -------------------------------------------------------------------------
###

#' @import SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
setClass('SummarizedLoomExperiment',
    contains=c('SummarizedExperiment', 'LoomExperiment')
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

setValidity2('SummarizedLoomExperiment', .valid.Experiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_SummarizedLoomExperiment <-
    function(se, colGraphs, rowGraphs)
{
    new('SummarizedLoomExperiment', se, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

#' @export
#' @importFrom methods as
SummarizedLoomExperiment <- function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)
    if (length(te) > 0 && is(te[[1]], 'SummarizedExperiment'))
        se <- te[[1]]
    else {
        se <- SummarizedExperiment(...)
        if(is(se, 'RangedSummarizedExperiment'))
            se <- as(se, 'SummarizedExperiment')
    }
    .new_SummarizedLoomExperiment(se, colGraphs=colGraphs, rowGraphs=rowGraphs)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_SummarizedExperiment_to_SummarizedLoomExperiment <- function(from)
{
    .new_SummarizedLoomExperiment(from,
                        rowGraphs=LoomGraphs(),
                        colGraphs=LoomGraphs())
}

setAs('SummarizedExperiment', 'SummarizedLoomExperiment',
    .from_SummarizedExperiment_to_SummarizedLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

#' @export
setMethod('[', c('SummarizedLoomExperiment', 'ANY', 'ANY'), .subset.LoomExperiment)

#' @export
setMethod('show', 'SummarizedLoomExperiment', .show.LoomExperiment)

