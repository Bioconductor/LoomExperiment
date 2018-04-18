
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
    representation(
        "VIRTUAL",
        colGraphs="LoomGraphs",
        rowGraphs="LoomGraphs"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.LoomExperiment <- function(x)
{
    NULL
}

setValidity2("LoomExperiment", .valid.LoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

.get.colGraphs <- function(x, ...)
{
    x@colGraphs
}

setMethod("colGraphs", "LoomExperiment", .get.colGraphs)

.replace.colGraphs <- function(x, ..., value)
{
    BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
}

setReplaceMethod("colGraphs", "LoomExperiment", .replace.colGraphs)

.get.rowGraphs <- function(x, ...)
{
    x@rowGraphs
}

setMethod("rowGraphs", "LoomExperiment", .get.rowGraphs)

.replace.rowGraphs <- function(x, ..., value)
{
    BiocGenerics:::replaceSlots(x, rowGraphs=value, check=FALSE)
}

setReplaceMethod("rowGraphs", "LoomExperiment", .replace.rowGraphs)

