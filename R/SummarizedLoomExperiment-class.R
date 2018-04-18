
### =========================================================================
### SummarizedLoomExperiment objects
### -------------------------------------------------------------------------
###

#' SummarizedLoomExperiment
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

setClass("SummarizedLoomExperiment",
    contains=c("SummarizedExperiment", "LoomExperiment")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.SummarizedLoomExperiment <- function(x)
{
    NULL
}

setValidity2("SummarizedLoomExperiment", .valid.SummarizedLoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_SummarizedLoomExperiment <-
    function(se, colGraphs, rowGraphs)
{
    new("SummarizedLoomExperiment", se, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

SummarizedLoomExperiment <- function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)[[1]]
    if (is(te, "SummarizedExperiment"))
        .new_SummarizedLoomExperiment(te, colGraphs=colGraphs, rowGraphs=rowGraphs)
    else {
        se <- SummarizedExperiment(...)
        if(is(se, "RangedSummarizedExperiment"))
            se <- as(se, "SummarizedExperiment")
        .new_SummarizedLoomExperiment(se, colGraphs=colGraphs, rowGraphs=rowGraphs)
    }
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

setAs("SummarizedExperiment", "SummarizedLoomExperiment",
    .from_SummarizedExperiment_to_SummarizedLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

setMethod("[", c("SummarizedLoomExperiment", "ANY", "ANY"), .subset.LoomExperiment)

setMethod("show", "SummarizedLoomExperiment", .show.LoomExperiment)

