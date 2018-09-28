
### =========================================================================
### SingleCellLoomExperiment objects
### -------------------------------------------------------------------------
###

#' @import SingleCellExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
setClass('SingleCellLoomExperiment',
    contains=c('SingleCellExperiment', 'RangedLoomExperiment')
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_SingleCellLoomExperiment <- function(sce, colGraphs, rowGraphs)
{
    new('SingleCellLoomExperiment', sce, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

#' @export
SingleCellLoomExperiment <-
    function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
    te <- list(...)
    if (length(te) > 0 && is(te[[1]], 'SingleCellExperiment'))
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

setAs('SingleCellExperiment', 'SingleCellLoomExperiment',
    .from_SingleCellExperiment_to_SingleCellLoomExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Methods.
###

setMethod('[', 'SingleCellLoomExperiment', .subset.LoomExperiment)

setMethod('show', 'SingleCellLoomExperiment', .show.LoomExperiment)

