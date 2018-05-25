
## Miscellanious methods

.valid.LoomExperiment <- function(x)
{
    clgs <- colGraphs(x)
    rlgs <- rowGraphs(x)
    cols <- seq_len(dims(x)[[2]])
    rows <- seq_len(dims(x)[[1]])
    ## Check that no "a" or "b" columns in LoomGraphs lie outside of dimensions.
    txt <- "All LoomGraph objects in LoomExperiment reference a row in the LoomExperiment"
    for (lg in clgs) {
        if(!all(lg$a %in% rows & lg$b %in% rows))
            return(txt)
    }
    for (lg in rlgs) {
        if(!all(lg$a %in% cols & lg$b %in% cols))
            return(txt)
    }
    NULL
}

.subset.LoomExperiment <- function(x, i, j, ...)
{
    if (!missing(i))
        rowGraphs(x) <- endoapply(rowGraphs(x), function(y) y[i,])
    if (!missing(j))
        colGraphs(x) <- endoapply(colGraphs(x), function(y) y[i,])
    callNextMethod()
}

.show.LoomExperiment <- function(object)
{
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    callNextMethod()
    if (length(object@rowGraphs) > 0)
        scat("rowGraphs(%d): %s\n", names(object@rowGraphs))
    else
        cat("rowGraphs(0): NULL\n")
    if (length(object@colGraphs) > 0)
        scat("colGraphs(%d): %s\n", names(object@colGraphs))
    else
        cat("colGraphs(0): NULL\n")
}

