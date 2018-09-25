### =========================================================================
### LoomGraph Objects
### -------------------------------------------------------------------------
###

#' @importFrom S4Vectors SelfHits
#' @export
setClass('LoomGraph',
    contains = 'SelfHits'
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
#    if (length(x) == 0)
#        return(NULL)
    mcol <- mcols(x)
    if (!is.integer(c(from(x), to(x)))) {
        txt <- sprintf('\n The nodes of a LoomGraph must be an integer')
        return(txt)
    }
    if (min(from(x), to(x)) < 0) {
        txt <- sprintf('\n The nodes of a LoomGraph must be non-negative')
        return(txt)
    }
    if (!is.null(mcol) && !all(names(mcol) == 'w')) {
        txt <- sprintf('\n A LoomGraph may only have one mcol named "w"')
        return(txt)
    }
    if (!is.null(w <- mcol$w) && !is.numeric(w)) {
        txt <- sprintf('\n The "w" mcol of a LoomGraph must numeric ')
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
LoomGraph <- function(from, to, nnode=max(from, to), ..., weight=NULL) {
    if (!is.numeric(c(from, to)))
        stop('"from" and "to" arguments to LoomGraph constructor  must be numeric')
    sh <- SelfHits(from=as.integer(from), to=as.integer(to), nnode=nnode, ...)
    mcols(sh)$w <- weight
    .new_LoomGraph(sh)
}

#' @export
LoomGraphs <- function(...) {
    list <- list(...)
    .new_LoomGraphs(list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_DataFrame_to_LoomGraph <- function(from)
{
    nam <- names(from)
    if (!all(nam %in% c('a', 'b', 'w')))
        stop('columns of DataFrame must be named "a" and "b" or "a", "b", "w"')
    a <- as.integer(from$a)
    b <- as.integer(from$b)
    lg <- LoomGraph(a, b, max(a,b))
    if('w' %in% nam)
        mcols(lg)$w <- as.numeric(from$w)
    lg
}

setAs('DataFrame', 'LoomGraph', .from_DataFrame_to_LoomGraph)

.from_LoomGraph_to_DataFrame <- function(from)
{
    a <- from(from)
    b <- to(from)
    if (is.null(w <- mcols(from)))
        DataFrame(a=a, b=b)
    else
        DataFrame(a=a, b=b, w)
}

setAs('LoomGraph', 'DataFrame', .from_LoomGraph_to_DataFrame)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellanious methods.
###

#' @export
setMethod('dropHits', c('LoomGraph', 'ANY'),
    function(x, i, ..., drop=TRUE)          
{
    ii <- .convert_subset_index(i, rownames(x))
    subset(x, x[['a']] %in% ii & x[['b']] %in% ii)   
})

