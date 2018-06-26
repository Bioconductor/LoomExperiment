
## Miscellanious methods

.valid.Experiment <- function(x)
{
    clgs <- colGraphs(x)
    rlgs <- rowGraphs(x)
    cols <- c(0, seq_len(dim(x)[[2]]))
    rows <- c(0, seq_len(dim(x)[[1]]))
    ## Check that no "a" or "b" columns in LoomGraphs lie outside of dimensions.
    col_test <- lapply(clgs, function(lg) {
        txt <- "All LoomGraph objects in LoomExperiment must reference a column in the LoomExperiment"
        test <- lg$a %in% cols & lg$b %in% cols
        if ((length(test) == 0 || !all(test)) && length(lg$a) > 0)
            return(txt)
    })
    row_test <- lapply(rlgs, function(lg) {
        txt <- "All LoomGraph objects in LoomExperiment must reference a row in the LoomExperiment"
        test <- lg$a %in% rows & lg$b %in% rows
        if ((length(test) == 0 || !all(test)) && length(lg$a) > 0)
            return(txt)
    })
    res <- list(col_test, row_test)
    unlist(res)
}

#' @importFrom S4Vectors endoapply
.subset.LoomExperiment <- function(x, i, j, ...)
{
    if (!missing(i))
        rowGraphs(x) <- endoapply(rowGraphs(x), function(y) y[i,])
    if (!missing(j))
        colGraphs(x) <- endoapply(colGraphs(x), function(y) y[,j])
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

