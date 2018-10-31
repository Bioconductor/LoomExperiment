
#' @importFrom S4Vectors endoapply
.subset.LoomExperiment <- function(x, i, j, ...)
{
    rg <- rowGraphs(x)
    cg <- colGraphs(x)
    rowGraphs(x) <- LoomGraphs()
    colGraphs(x) <- LoomGraphs()
    x <- callNextMethod()
    if (!missing(i)) {
        if(all(i > 0))
            rowGraphs(x) <- .change.nnode(endoapply(rg, function(y) loomSelectHits(y, i)), nrow(x))
        else
            rowGraphs(x) <- .change.nnode(endoapply(rg, function(y) loomDropHits(y, abs(i))), nrow(x))
    }
    if (!missing(j)) {
        if(all(j > 0))
            colGraphs(x) <- .change.nnode(endoapply(cg, function(y) loomSelectHits(y, j)), ncol(x))
        else
            colGraphs(x) <- .change.nnode(endoapply(cg, function(y) loomDropHits(y, abs(j))), ncol(x))
    }
    x
}

.rbind.LoomExperiment <-
    function(..., deparse.level = 1)
{
    li <- list(...)
    rn <- names(rowGraphs(li[[1]]))

    clgs <- lapply(li, colGraphs)
    clgs <- do.call(c, clgs)

    rlgs <- lapply(li, rowGraphs)
    rlgs <- do.call(rbind, rlgs)
    if (is(rlgs, "matrix"))
        rlgs <- LoomGraphs()
    names(rlgs) <- rn
    x <- callNextMethod()
    rowGraphs(x) <- .change.nnode(rlgs, nrow(x))
    colGraphs(x) <- clgs
    x
}

.show.LoomExperiment <- function(object)
{
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse=' ')
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep='\n')
    }
    callNextMethod()
    if (length(object@rowGraphs) > 0) {
        if (is.null(names(object@rowGraphs)))
            cat(sprintf('rowGraphs(%d):\n', length(object@rowGraphs)))
        else
            scat('rowGraphs(%d): %s\n', names(object@rowGraphs))
    }
    else
        cat('rowGraphs(0): NULL\n')
    if (length(object@colGraphs) > 0) {
        if (is.null(names(object@rowGraphs)))
            cat(sprintf('colGraphs(%d):\n', length(object@rowGraphs)))
        else
            scat('colGraphs(%d): %s\n', names(object@colGraphs))
    }
    else
        cat('colGraphs(0): NULL\n')
}
