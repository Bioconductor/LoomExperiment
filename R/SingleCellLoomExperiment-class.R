
### =========================================================================
### SingleCellLoomExperiment objects
### -------------------------------------------------------------------------
###

#' @import SingleCellExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
setClass("SingleCellLoomExperiment",
    contains=c("SingleCellExperiment", "LoomExperiment")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

setValidity2("SingleCellLoomExperiment", .valid.Experiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_SingleCellLoomExperiment <- function(sce, colGraphs, rowGraphs)
{
    new("SingleCellLoomExperiment", sce, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

#' @export
SingleCellLoomExperiment <-
    function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)
    if (length(te) > 0 && is(te[[1]], "SingleCellExperiment"))
        sce <- te[[1]]
    else
        sce <- SingleCellExperiment(...)
    .new_SingleCellLoomExperiment(sce, colGraphs=colGraphs, rowGraphs=rowGraphs)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_SingleCellExperiment_to_SingleCellLoomExperiment <- function(from)
{
    .new_SingleCellLoomExperiment(from,
                                  colGraphs=LoomGraphs(),
                                  rowGraphs=LoomGraphs())
}

setAs("SingleCellExperiment", "SingleCellLoomExperiment",
    .from_SingleCellExperiment_to_SingleCellLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

#' @export
setMethod("[", c("SingleCellLoomExperiment", "ANY", "ANY"), .subset.LoomExperiment)

#' @export
setMethod("show", "SingleCellLoomExperiment", .show.LoomExperiment)

