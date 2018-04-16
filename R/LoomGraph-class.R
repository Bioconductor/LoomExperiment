### =========================================================================
### LoomGraph Objects
### -------------------------------------------------------------------------
###

#' @export
.LoomGraph <- setClass("LoomGraph",
    contains = "DataFrame"
)

#' @export
.LoomGraphs <- setClass("LoomGraphs",
    contains = "SimpleList",
    prototype = prototype(
        elementType = "LoomGraph"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Vaidity
###

.valid.LoomGraph <- function(x) {
    if (isEmpty(x))
        return(NULL)
    numerals <- vapply(x, is.numeric, logical(1))
    if (!all(numerals)) {
        txt <- sprintf("\n A LoomGraph must only contain numeric elements")
        return(txt)
    }
    cols <- colnames(x)
    if (length(cols) == 2) {
        if(!all(cols == c('a', 'b'))) {
            txt <- sprintf("\n A LoomGraph with two columns must be named 'a' and 'b'")
            return(txt)
        }
    } else if (length(cols) == 3) {
        if(!all(cols == c('a', 'b', 'w'))) {
            txt <- sprintf("\n A LoomGraph with three columns must be named 'a', 'b', and 'w'")
            return(txt)
        }
    } else {
        txt <- sprintf("\n A LoomGraph must have two or three columns")
        return(txt)
    }
    NULL
}

setValidity2("LoomGraph", .valid.LoomGraph)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

#'
#' LoomGraph class for LoomExperiment
#'
#' @description The LoomGraph class is a DataFrame-like class that represents a
#'  weighted or unweighted graph from the loom file format.
#'
#' @details The LoomGraph object must have two or three columns. If the graph
#'  has two columns, it is an unweighted graph and must have columns named 'a'
#'  and 'b'. If the graph has three columns, it is treated as a weighted graph
#'  and must have column names 'a', 'b', and 'w'.  Each row in the LoomGraph
#'  represents an edge between two column/row numbers and may have a weight.
#'
#' @examples
#'  ## create new LoomGraph
#'  lg <- LoomGraph(a=c(1, 2, 3), b=c(3, 2, 1), w=c(5, 4, 1))
#'  lg
#'
#' @rdname LoomGraph-class
#' @export
LoomGraph <- function(...) {
    df <- DataFrame(...)
    .LoomGraph(df)
}

#'
#' LoomGraphs class for LoomExperiment
#'
#' @description The LoomGraphs class contain multiple LoomGraph objects to
#'  display either colGraph or rowGraph data from the loom file format.
#'
#' @examples
#'  ## create new LoomGraphs
#'  lg1 <- LoomGraph(a=c(1, 2, 3), b=c(3, 2, 1), w=c(3, 4, 5))
#'  lg2 <- LoomGraph(a=c(3, 5, 2), b=c(5, 3, 3))
#'  lgs <- LoomGraphs(lg1, lg2)
#'  lgs
#'
#' @rdname LoomGraphs-class
#' @export
LoomGraphs <- function(...) {
    sl <- SimpleList(...)
    .LoomGraphs(sl)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellanious methods
###

#'
#' @rdname LoomGraph-class
#' @exportMethod [
setMethod("[", c("LoomGraph", "ANY", "missing"),
    function(x, i, j, ..., drop=TRUE)          
{
    ii <- .convert_subset_index(i, rownames(x))
    x <- as(x, "DataFrame")
    x <- x[x[["a"]] %in% ii & x[["b"]] %in% ii,]  
    as(x, "LoomGraph")
})

#'
#' importFrom plyr mapvalues
#'
#' @rdname LoomGraphs-class
#' @exportMethod [
setReplaceMethod("[", c("LoomGraph", "ANY", "missing", "numeric"),
    function(x, i, j, ..., value)          
{
    ii <- .convert_subset_index(i, rownames(x))

    ld <- x@listData
    res <- lapply(ld[c('a', 'b')], function(x){
        mapvalues(x, ii, value)
    })
    ld[c('a', 'b')] <- res
    x@listData <- ld
})
