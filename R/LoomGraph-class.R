### =========================================================================
### LoomGraph Objects
### -------------------------------------------------------------------------
###

#' @export
setClass('LoomGraph',
    contains = 'DataFrame'
)

#' @export
setClass('LoomGraphs',
    contains = 'SimpleList',
    prototype = prototype(
        elementType = 'LoomGraph'
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Vaidity
###

#' @importFrom S4Vectors isEmpty
.valid.LoomGraph <- function(x) {
    if (isEmpty(x))
        return(NULL)
    integers <- vapply(x, is.numeric, logical(1))
    if (!all(integers)) {
        txt <- sprintf('\n A LoomGraph must only contain integer elements')
        return(txt)
    }
    cols <- colnames(x)
    if (length(cols) == 2) {
        if(!all(cols == c('a', 'b'))) {
            txt <- sprintf('\n A LoomGraph with two columns must be named "a" and "b"')
            return(txt)
        }
    } else if (length(cols) == 3) {
        if(!all(cols == c('a', 'b', 'w'))) {
            txt <- sprintf('\n A LoomGraph with three columns must be named "a", "b", and "w"')
            return(txt)
        }
    } else {
        txt <- sprintf('\n A LoomGraph must have two or three columns')
        return(txt)
    }
    NULL
}

setValidity2('LoomGraph', .valid.LoomGraph)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

#' @importFrom methods new
.new_LoomGraph <- function(df) {
    new('LoomGraph', df)
}

.new_LoomGraphs <- function(li) {
    new('LoomGraphs', listData = li)
}

#' @export
LoomGraph <- function(...) {
    df <- DataFrame(...)
    .new_LoomGraph(df)
}

#' @export
LoomGraphs <- function(...) {
    #sl <- SimpleList(...)
    list <- list(...)
    .new_LoomGraphs(list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellanious methods
###

#' @export
setMethod('[', c('LoomGraph', 'ANY', 'missing'),
    function(x, i, j, ..., drop=TRUE)          
{
    ii <- .convert_subset_index(i, rownames(x))
    subset(x, x[['a']] %in% ii & x[['b']] %in% ii)   
})

#' @export
#' @importFrom plyr mapvalues
setReplaceMethod('[', c('LoomGraph', 'ANY', 'missing', 'numeric'),
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
