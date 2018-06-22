
#' @importFrom rhdf5 h5write
setMethod(".exportLoom", "matrix",
    function(object, con, name)
{
    object <- t(object)
    tryCatch({
        rhdf5::h5write(object, con, name)
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

#' @importFrom DelayedArray DelayedArray
#' @importFrom HDF5Array writeHDF5Array
setMethod(".exportLoom", "DelayedArray",
    function(object, con, name)
{
    HDF5Array::writeHDF5Array(t(object), con, name)
    0L
})

setMethod(".exportLoom", "vector",
    function(object, con, name, rowname_attr)
{
    object <- as.matrix(object)
    .exportLoom(object, con, name)
})

#' @importFrom rhdf5 h5write
setMethod(".exportLoom", "data.frame",
    function(object, con, name, rowname_attr)
{
    if (!is.null(rowname_attr))
        object[[rowname_attr]] <- rownames(object)

    is.factor <- vapply(object, is, logical(1), "factor")
    #if (any(is.factor))
        #warning(
        #    "'.exportLoom()' coerced 'factor' column(s) to character:\n  ",
        #    paste(sQuote(names(object)[is.factor]), collapse=", "),
        #    call. = FALSE
        #)
    object[is.factor] <- lapply(object[is.factor], as.character)

    names <- sprintf("/%s/%s", name, names(object))
    tryCatch({
        Map(rhdf5::h5write, object, names, MoreArgs = list(file = con))
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

setMethod(".exportLoom", "DataFrame",
    function(object, con, name, rowname_attr)
{
    object <- as.data.frame(object)
    .exportLoom(object, con, name, rowname_attr)
})

#' @import GenomicRanges
setMethod(".exportLoom", "GenomicRanges",
    function(object, con, name, rowname_attr)
{
    object <- as.data.frame(object)
    name <- paste0(name, "/GRanges")
    rhdf5::h5createGroup(con, name)
    names <- colnames(object)
    colnames(object) <- names
    .exportLoom(object, con, name, rowname_attr)
})

.get_empty_GRangesList_value <- function(type) {
    switch(type, "character" = '', "numeric" = 0, "integer" = 0, "double" = 0, '')
}

setMethod(".exportLoom", "GenomicRangesList",
    function(object, con, name, rowname_attr)
{
    lengths <- lengths(object)
    num <- length(object)
    max <- max(lengths)
    if (max == 0)
        max <- 1
    names <- names(object)
    if(is.null(names))
        names <- rep('', length(object))

    name <- paste0(name, "/GRangesList")
    rhdf5::h5createGroup(con, name)

    .exportLoom(lengths, con, paste0(name, "/lengths"), rowname_attr)
    .exportLoom(names, con, paste0(name, "/names"), rowname_attr)

    rownames <- unlist(lapply(object, function(x) rownames(as.data.frame(x))))

    df <- as.data.frame(object)
    df['rownames'] <- rownames

    names <- colnames(df)
    names <- names[!names %in% c("group")]
    
    na_types <- vapply(df, class, character(1))
    names(na_types) <- names

    dfs <- lapply(names, function(i) {
        val <- lapply(seq_len(num), function(idx) {
            na <- .get_empty_GRangesList_value(na_types[[i]])
            if(!idx %in% df$group)
                rep_len(na, max)
            else {
                temp <- df[df$group==idx,i]
                if(all(is.na(temp))) 
                    rep_len(na, max) 
                else {
                    temp <- rep_len(temp, max)
                    if(is(temp, "factor"))
                        temp <- as(temp, "character")
                    temp
                }
            }
        })
        do.call(rbind, val)
    })

    df_names <- paste0(name, '/', names)
    Map(.exportLoom, dfs, name = df_names, MoreArgs = list(con = con))
})

setMethod(".exportLoom", "LoomGraph", 
    function(object, con, name)
{
    rhdf5::h5createGroup(con, name)
    name <- paste0(name, '/', colnames(object))
    tryCatch({
        Map(rhdf5::h5write, object, name, MoreArgs = list(file = con))
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

setMethod(".exportLoom", "LoomGraphs", 
    function(object, con, name)
{
    rhdf5::h5createGroup(con, name)
    name <- paste0(name, '/', names(object))
    Map(.exportLoom, object, name = name, MoreArgs = list(con = con))
})

.exportLoom.LoomExperiment <-
        function(object, con,
             matrix = assayNames(object)[1],
             rownames_attr = "rownames", colnames_attr = "colnames")
{
    con <- path(con)

    stopifnot(
        !file.exists(con),
        is.character(matrix), length(matrix) == 1L, !is.na(matrix),
        matrix %in% assayNames(object),
        is.character(rownames_attr), length(rownames_attr) == 1L,
        !is.na(rownames_attr),
        is.character(colnames_attr), length(colnames_attr) == 1L,
        !is.na(colnames_attr)
    )
    
    if (!is.null(rownames(object)) && rownames_attr %in% names(rowData(object)))
        stop("'rownames_attr' must not be in names(rowData())")
    if (!is.null(colnames(object)) && colnames_attr %in% names(colData(object)))
        stop("'colnames_attr()' must not be in names(colData())")

    rhdf5::h5createFile(con)

    h5f <- H5Fopen(con) 
    tryCatch({
        rhdf5::h5writeAttribute(paste0("LoomExperiment-", as.character(
            packageVersion("LoomExperiment"))), name="CreatedWith", h5obj=h5f)
        rhdf5::h5writeAttribute(class(object), name="LoomExperiment-class",
            h5obj=h5f)
        Map(rhdf5::h5writeAttribute, metadata(object),
            name = names(metadata(object)), MoreArgs = list(h5obj = h5f))
    }, error = function(err) {
        warning(conditionMessage(err))
    }, finally = H5Fclose(h5f))

    assays <- assays(object, withDimnames = FALSE)
    layers <- setNames(paste0("/layers/", names(assays)), names(assays))
    layers[matrix] <- "/matrix"

    if (length(layers) > 1L)
        rhdf5::h5createGroup(con, "/layers")
    success <- unlist(Map(
        .exportLoom, assays, name = layers, MoreArgs = list(con = con)
    ))
    if (!all(success == 0L))
        stop(
            "'.exportLoom()' failed to write assay(s)\n  ",
            paste0(sQuote(names(layers)[success != 0]), collapse = ", ")
        )

    rhdf5::h5createGroup(con, "/col_attrs")
    rhdf5::h5createGroup(con, "/row_attrs")

    if (is(object, "SingleCellLoomExperiment")) {
        reducedDims_names <- "/col_attrs/reducedDims"
        rhdf5::h5createGroup(con, reducedDims_names)
        reducedDims_names <- paste0(reducedDims_names, '/', names(reducedDims(object)))
        if (length(reducedDims(object)) == 0)
            reducedDims_names <- character(0)
        Map(.exportLoom, reducedDims(object), name = reducedDims_names, MoreArgs = list(con = con))

        #rhdf5::h5createGroup(con, "/row_attrs/int_colData")
        #.exportLoom(object@int_colData, con, "row_attrs/int_colData", NULL)

        #rhdf5::h5createGroup(con, "/row_attrs/int_elementMetadata")
        #.exportLoom(object@int_elementMetadata, con, "row_attrs/int_elementMetadata", NULL)

        #int_metadata_names <- "/row_attrs/int_metadata"
        #rhdf5::h5createGroup(con, int_metadata_names)
        #int_metadata_names <- paste0(int_metadata_names, '/', names(object@int_metadata))
        #int_metadata <- lapply(object@int_metadata, as.character)
        #Map(rhdf5::h5write, obj = int_metadata, name = int_metadata_names, MoreArgs = list(file = con))
    }

    .exportLoom(colData(object), con, "col_attrs", colnames_attr)
    rowData <- rowData(object)
    if (is(object, "RangedSummarizedExperiment")) {
        rowRanges <- rowRanges(object)
        if (is(rowRanges, "GRangesList")) {
#            GRangesList_lengths <- lengths(rowRanges)
#            rowData <- cbind(rowData, GRangesList_lengths)
            .exportLoom(rowRanges, con, "row_attrs", rownames_attr)
        } else {
            .exportLoom(rowRanges, con, "row_attrs", rownames_attr)
        }
    }
    else
        .exportLoom(rowData, con, "row_attrs", rownames_attr)

    if (length(colGraphs(object)) > 0)
        .exportLoom(colGraphs(object), con, "col_graphs")

    if (length(rowGraphs(object)) > 0)
        .exportLoom(rowGraphs(object), con, "row_graphs")

    invisible(con)
}

#' @importFrom rhdf5 h5createGroup
#' @importFrom rtracklayer export
#' @export
setMethod("export", signature=c("LoomExperiment", "loomFile", "ANY"),
    .exportLoom.LoomExperiment)

