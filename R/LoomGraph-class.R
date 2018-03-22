### =========================================================================
### LoomGraph Objects
### -------------------------------------------------------------------------
###

#' @export
.LoomGraph <- setClass(
    "LoomGraph",
    contains = "DataFrame"
)

#' @export
.LoomGraphs <- setClass(
    "LoomGraphs",
    contains = "SimpleList",
    prototype = prototype(
        elementType = "LoomGraph"
    )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Vaidity
###

.LoomGraph.validity <- function(x) {
    if (isEmpty(x))
        return(NULL)
    numerals <- vapply(df, is.numeric, logical(1))
    if (!all(numerals)) {
        txt <- sprintf("\n LoomGraph must only contain numeric elements")
        return(txt)
    }
    cols <- colnames(x)
    if (length(cols) == 2) {
        if(!all(cols == c('a', 'b'))) {
            txt <- sprintf("\n LoomGraph with two columns must be named 'a' and 'b'")
            return(txt)
        }
    } else if (length(cols) == 3) {
        if(!all(cols == c('a', 'b'))) {
            txt <- sprintf("\n LoomGraph with two columns must be named 'a' and 'b'")
            return(txt)
        }
    } else {
        txt <- sprintf("\n LoomGraph must have two or three columns")
        return(txt)
    }
    NULL
}

#setValidity2("LoomGraph", .LoomGraph.validity)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

#setGeneric("LoomGraph", function(graph, ...), standardGeneric("LoomGraph"))

#setMethod("LoomGraph", "DataFrame",
#    function(graph, ...)
#{
#    as(lg, "LoomGraph")
#})

LoomGraph <- function(..., row.names = NULL, check.names = TRUE) {
    df <- DataFrame(..., row.names = row.names, check.names = check.names)
    as(df, "LoomGraph")
}

#setGeneric("LoomGraphs", function(graph, ...), standardGeneric("LoomGraphs"))

#setMethod("LoomGraph", "DataFrame",
#    function(graph, ...)
#{
#    lg <- .LoomGraph()
#    as(graph, "DataFrame")
#})

LoomGraphs <- function(...) {
    .LoomGraphs(list(...))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

#.from_LoomGraph_to_DataFrame <- function(from)
#{
    
#}

#setAs("LoomGraph", "DataFrame", .from_LoomGraph_to_DataFrame)

#.from_DataFrame_to_LoomGraph <- function(from)
#{
    
#}

#setAs("DataFrame", "LoomGraph", .from_DataFrame_to_LoomGraph)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Misc.
###

setMethod("[", c("LoomGraph", "ANY", "missing"),
    function(x, i, j, ..., drop=TRUE)          
{
    ii <- .convert_subset_index(i, rownames(x))
    subset(x, x[["a"]] %in% ii & x[["b"]] %in% ii)   
})

setReplaceMethod("[", c("LoomGraph", "ANY", "missing", "LoomGraph"),
    function(x, i, j, ..., value)          
{
    ii <- .convert_subset_index(i, rownames(x))

    rgx <- x[i,]

    ## update 'rgx' by
    ## - removing 'i' nodes from rgx
    ## - adding rgv to rgx

    rgx <- rbind(rgx, value)

    x@rownames <- rgx@rownames
    x@nrows <- rgx@nrows
    x@listData <- rgx@listData
    x
})
