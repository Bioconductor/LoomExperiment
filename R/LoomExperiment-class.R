### =========================================================================
### LoomExperiment objects
### -------------------------------------------------------------------------
###

#' LoomExperiment
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

setClass("LoomExperiment",
    contains="SummarizedExperiment",
    representation(
        colGraphs="LoomGraphs",
        rowGraphs="LoomGraphs"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

setValidity2("LoomExperiment", .valid.LoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_LoomExperiment <-
    function(se, colGraphs, rowGraphs)
{
    new("LoomExperiment", se, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

LoomExperiment <- function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)[[1]]
    if (is(te, "SummarizedExperiment"))
        .new_LoomExperiment(te, colGraphs=colGraphs, rowGraphs=rowGraphs)
    else {
        se <- SummarizedExperiment(...)
        if(is(se, "RangedSummarizedExperiment"))
            se <- as(se, "SummarizedExperiment")
        .new_LoomExperiment(se, colGraphs=colGraphs, rowGraphs=rowGraphs)
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

setAs("SummarizedExperiment", "LoomExperiment",
    .from_SummarizedExperiment_to_LoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

setMethod("colGraphs", "LoomExperiment", .get.colGraphs)

setReplaceMethod("colGraphs", "LoomExperiment", .replace.colGraphs)

setMethod("rowGraphs", "LoomExperiment", .get.rowGraphs)

setReplaceMethod("rowGraphs", "LoomExperiment", .replace.rowGraphs)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

setMethod("[", c("LoomExperiment", "ANY", "ANY"), .subset.LoomExperiment)

setMethod("show", "LoomExperiment", .show.LoomExperiment)
