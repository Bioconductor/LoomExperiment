### =========================================================================
### LoomGraph Objects
### -------------------------------------------------------------------------
###

#' @export
setClass("LoomGraph",
    contains = "DataFrame"
)

#' @export
setClass("LoomGraphs",
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
        if(!all(cols == c('a', 'b'))) {
            txt <- sprintf("\n A LoomGraph with two columns must be named 'a' and 'b'")
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

.new_LoomGraph <- function(df) {
    new("LoomGraph", df)
}

.new_LoomGraphs <- function(sl) {
    new("LoomGraphs", sl)
}

LoomGraph <- function(...) {
    df <- DataFrame(...)
    .new_LoomGraph(df)
}

LoomGraphs <- function(...) {
    sl <- SimpleList(...)
    .new_LoomGraphs(sl)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellanious methods
###

setMethod("[", c("LoomGraph", "ANY", "missing"),
    function(x, i, j, ..., drop=TRUE)          
{
    ii <- .convert_subset_index(i, rownames(x))
    subset(x, x[["a"]] %in% ii & x[["b"]] %in% ii)   
})

#' importFrom plyr mapvalues
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
