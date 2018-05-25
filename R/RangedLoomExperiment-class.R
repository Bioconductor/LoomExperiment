
### =========================================================================
### RangedLoomExperiment objects
### -------------------------------------------------------------------------
###

#' @import SummarizedExperiment
#' @export
setClass("RangedLoomExperiment",
    contains=c("RangedSummarizedExperiment", "LoomExperiment")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

setValidity2("RangedLoomExperiment", .valid.LoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_RangedLoomExperiment <- function(se, colGraphs, rowGraphs)
{
    new("RangedLoomExperiment", se, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

#' @export
RangedLoomExperiment <-
    function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)[[1]]
    if (is(te, "SummarizedExperiment"))
        .new_RangedLoomExperiment(te, colGraphs=colGraphs, rowGraphs=rowGraphs)
    else {    
        rse <- SummarizedExperiment(...)
        if(!is(rse, "RangedSummarizedExperiment"))
            rse <- as(rse, "RangedSummarizedExperiment")
        .new_RangedLoomExperiment(rse, colGraphs=colGraphs, rowGraphs=rowGraphs)
    }
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

setAs("RangedSummarizedExperiment", "RangedLoomExperiment",
    .from_RangedSummarizedExperiment_to_RangedLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

#' @export
setMethod("[", c("RangedLoomExperiment", "ANY", "ANY"), .subset.LoomExperiment)

#' @export
setMethod("show", "RangedLoomExperiment", .show.LoomExperiment)

