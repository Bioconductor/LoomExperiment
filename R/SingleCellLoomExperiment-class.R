
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

.valid.SingleCellLoomExperiment <- function(x)
{
    NULL
}

setValidity2("SingleCellLoomExperiment", .valid.SingleCellLoomExperiment)


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
    te <- list(...)[[1]]
    if (is(te, "SummarizedExperiment"))
        .new_SingleCellLoomExperiment(te, colGraphs=colGraphs, rowGraphs=rowGraphs)
    else {
        sce <- SingleCellExperiment(...)
        .new_SingleCellLoomExperiment(sce, colGraphs=colGraphs, rowGraphs=rowGraphs)
    }
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

