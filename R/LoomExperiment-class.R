
### =========================================================================
### LoomExperiment objects
### -------------------------------------------------------------------------
###

#' @import SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export

setClass('LoomExperiment',
    representation(
        'VIRTUAL',
        colGraphs='LoomGraphs',
        rowGraphs='LoomGraphs'
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.LoomExperiment <- function(x)
{
    NULL
}

setValidity2('LoomExperiment', .valid.LoomExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get and Replace methods.
###

.get.colGraphs <- function(x, ...)
{
    x@colGraphs
}

#' @export
setMethod('colGraphs', 'LoomExperiment', .get.colGraphs)

#' @importFrom methods validObject callNextMethod
.replace.colGraphs <- function(x, ..., value)
{
    x <- BiocGenerics:::replaceSlots(x, colGraphs=value, check=FALSE)
    validObject(x)
    x
}

#' @export
setReplaceMethod('colGraphs', 'LoomExperiment', .replace.colGraphs)

.get.rowGraphs <- function(x, ...)
{
    x@rowGraphs
}

#' @export
setMethod('rowGraphs', 'LoomExperiment', .get.rowGraphs)

.replace.rowGraphs <- function(x, ..., value)
{
    x <- BiocGenerics:::replaceSlots(x, rowGraphs=value, check=FALSE)
    validObject(x)
    x
}

#' @export
setReplaceMethod('rowGraphs', 'LoomExperiment', .replace.rowGraphs)

