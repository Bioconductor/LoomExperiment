
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
#' @export
setClass("RangedLoomExperiment",
    contains=c("RangedSummarizedExperiment", "LoomExperiment")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.RangedLoomExperiment <- function(x)
{
    NULL
}

setValidity2("RangedLoomExperiment", .valid.RangedLoomExperiment)


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

setMethod("[", c("RangedLoomExperiment", "ANY", "ANY"), .subset.LoomExperiment)

setMethod("show", "RangedLoomExperiment", .show.LoomExperiment)

