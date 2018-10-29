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
        txt <- sprintf('\n A LoomGraph may only have one metadata column named "w"')
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

#' @importFrom S4Vectors from to
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
setMethod('rbind', 'LoomGraph',
    function(..., deparse.level = 1)
{
    x <- c(...)
    nnode <- max(from(x), to(x))
    x@nLnode <- nnode
    x@nRnode <- nnode
    x
})

setMethod('rbind', 'LoomGraphs',
    function(..., deparse.level = 1)
{
    li <- list(...)
    res <- lapply(seq_along(li[[1]]), function(i) {
        lg <- lapply(li, function(x) {
            x[[i]]
        })
        do.call(rbind, lg)
    })
    do.call(LoomGraphs, res)
})

.correctHits <- function(x, i, decr)
{
    ifelse(x >= i, decr-1, decr)
}

setMethod('loomSelectHits', c('LoomGraph', 'ANY'),
    function(x, i, ...)          
{
    x <- x[from(x) %in% i & to(x) %in% i]
    from_top <- seq_len(max(i))
    i <- from_top[!from_top %in% i]

    from <- from(x)
    from_i <- rep(0, length(from))
    for (n in i)
        from_i <- .correctHits(from, n, from_i)
    from <- from + from_i
    to <- to(x)
    to_i <- rep(0, length(to))
    for (n in i)
        to_i <- .correctHits(to, n, to_i)
    to <- to + to_i

    nnode <- max(from, to)
    LoomGraph(from, to, nnode, weight = mcols(x)[[1]])
})

setMethod('loomDropHits', c('LoomGraph', 'ANY'),
    function(x, i, ...)          
{
    x <- x[!from(x) %in% i & !to(x) %in% i]

    from <- from(x)
    from_i <- rep(0, length(from))
    for (n in i)
        from_i <- .correctHits(from, n, from_i)
    from <- from + from_i
    to <- to(x)
    to_i <- rep(0, length(to))
    for (n in i)
        to_i <- .correctHits(to, n, to_i)
    to <- to + to_i

    nnode <- max(from, to)
    LoomGraph(from, to, nnode, weight = mcols(x)[[1]])
})

setReplaceMethod('loomDropHits', c('LoomGraph', 'ANY', 'ANY'),
    function(x, i, ..., value)
{
    from <- from(x)
    to <- to(x)
    from <- as.integer(replace(from, from %in% i, value))
    to <- as.integer(replace(to, to %in% i, value))

    nnode <- max(from, to)
    x@nLnode <- nnode
    x@nRnode <- nnode
    x@from <- from
    x@to <- to
    x
})

