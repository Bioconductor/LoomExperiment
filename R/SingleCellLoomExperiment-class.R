### =========================================================================
### SingleCellLoomExperiment objects
### -------------------------------------------------------------------------
###

#' SingleCellLoomExperiment
#'
#' A class that helps facilitate the transition of SummarizedExperiment objects
#' to .loom files and vise versa.
#'
#' @slot colGraphs A SimpleList containing the colGraphs information
#' @slot rowGraphs A SimpleList containing the rowGraphs information
#'
#' @author Daniel Van Twisk
#' @import SingleCellExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export

setClass("SingleCellLoomExperiment",
    contains="SingleCellExperiment",
    representation(
        colGraphs="LoomGraphs",
        rowGraphs="LoomGraphs"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

setValidity2("SingleCellLoomExperiment", .valid.LoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_SingleCellLoomExperiment <- function(sce, colGraphs, rowGraphs)
{
    new("SingleCellLoomExperiment", sce, colGraphs=colGraphs, rowGraphs=rowGraphs)
}

SingleCellLoomExperiment <-
    function(..., colGraphs=LoomGraphs(), rowGraphs=LoomGraphs())
{
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
### Get and Replace methods.
###

setMethod("colGraphs", "SingleCellLoomExperiment", .get.colGraphs)

setReplaceMethod("colGraphs", "SingleCellLoomExperiment", .replace.colGraphs)

setMethod("rowGraphs", "SingleCellLoomExperiment", .get.rowGraphs)

setReplaceMethod("rowGraphs", "SingleCellLoomExperiment", .replace.rowGraphs)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellenious methods.
###

setMethod("[", c("SingleCellLoomExperiment", "ANY", "ANY"), .subset.LoomExperiment)

setMethod("show", "SingleCellLoomExperiment", .show.LoomExperiment)
