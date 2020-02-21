
#' @importFrom rhdf5 h5write
setMethod('.exportLoom', 'matrix',
    function(object, con, name)
{
    object <- t(object)
    tryCatch({
        rhdf5::h5write(object, con, name)
        0L
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

#' @importFrom DelayedArray DelayedArray
#' @importFrom HDF5Array writeHDF5Array
setMethod('.exportLoom', 'DelayedArray',
    function(object, con, name)
{
    HDF5Array::writeHDF5Array(t(object), con, name)
    0L
})

#' @importFrom Matrix t
setMethod('.exportLoom', 'dgCMatrix',
    function(object, con, name)
{
    HDF5Array::writeHDF5Array(Matrix::t(object), con, name)
    0L
})

setMethod('.exportLoom', 'vector',
    function(object, con, name, rowname_attr)
{
    object <- as.matrix(object)
    .exportLoom(object, con, name)
})

#' @importFrom rhdf5 h5write
setMethod('.exportLoom', 'data.frame',
    function(object, con, name, rowname_attr)
{
    rnames <- rownames(object)
    hasdfrownames <- identical(rownames(object),
        as.character(seq_along(rownames(object))))
    if (!is.null(rowname_attr) && !hasdfrownames)
        object[[rowname_attr]] <- rownames(object)

    is.factor <- vapply(object, is, logical(1), 'factor')
    object[is.factor] <- lapply(object[is.factor], as.character)

    names <- sprintf('/%s/%s', name, names(object))
    tryCatch({
        Map(rhdf5::h5write, object, names, MoreArgs = list(file = con))
        0L
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

setMethod('.exportLoom', 'DataFrame',
    function(object, con, name, rowname_attr)
{
    object <- as.data.frame(object)
    .exportLoom(object, con, name, rowname_attr)
})

#setMethod('.exportLoom', 'LoomGraph',
#    function(object, con, name, rowname_attr)
#{
#    object <- as(object, "DataFrame")
#    .exportLoom(object, con, name, rowname_attr)
#})

#' @import GenomicRanges
setMethod('.exportLoom', 'GenomicRanges',
    function(object, con, name, rowname_attr)
{
    object <- as.data.frame(object)
    name <- paste0(name, '/GRanges')
    rhdf5::h5createGroup(con, name)
    names <- colnames(object)
    colnames(object) <- names
    .exportLoom(object, con, name, rowname_attr)
})

.get_empty_GRangesList_value <- function(type) {
    switch(type, 'character' = '', 'numeric' = 0, 'integer' = 0, 'double' = 0, '')
}

setMethod('.exportLoom', 'GenomicRangesList',
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

    name <- paste0(name, '/GRangesList')
    rhdf5::h5createGroup(con, name)

    .exportLoom(lengths, con, paste0(name, '/lengths'), rowname_attr)
    .exportLoom(names, con, paste0(name, '/names'), rowname_attr)

    rownames <- unlist(lapply(object, function(x) rownames(as.data.frame(x))))

    df <- as.data.frame(object)
    df['rownames'] <- rownames

    names <- colnames(df)
    names <- names[!names %in% c('group')]

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
                    if(is(temp, 'factor'))
                        temp <- as(temp, 'character')
                    temp
                }
            }
        })
        do.call(rbind, val)
    })

    df_names <- paste0(name, '/', names)
    Map(.exportLoom, dfs, name = df_names, MoreArgs = list(con = con))
})

setMethod('.exportLoom', 'LoomGraph',
    function(object, con, name)
{
    rhdf5::h5createGroup(con, name)
    object <- as(object, "DataFrame")
    object <- as.matrix(object)
    object[,1] <- object[,1] - 1
    object[,2] <- object[,2] - 1
    object <- DataFrame(object)
    name <- paste0(name, '/', colnames(object))
    tryCatch({
        Map(rhdf5::h5write, object, name, MoreArgs = list(file = con))
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

setMethod('.exportLoom', 'LoomGraphs',
    function(object, con, name)
{
    rhdf5::h5createGroup(con, name)
    if (length(object) > 0) {
        name <- paste0(name, '/', names(object))
        Map(.exportLoom, object, name = name, MoreArgs = list(con = con))
    }
})

#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @importFrom rhdf5 H5Fclose H5Fopen
#' @importFrom rtracklayer path
#' @importFrom stats setNames
#' @importFrom utils packageVersion
.exportLoom.LoomExperiment <-
        function(object, con,
             matrix = assayNames(object)[1],
             rownames_attr = 'rownames', colnames_attr = 'colnames')
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
        stop('"rownames_attr" must not be in names(rowData())')
    if (!is.null(colnames(object)) && colnames_attr %in% names(colData(object)))
        stop('"colnames_attr" must not be in names(colData())')

    rhdf5::h5createFile(con)

    assays <- assays(object, withDimnames = FALSE)
    layers <- setNames(paste0('/layers/', names(assays)), names(assays))
    layers[matrix] <- '/matrix'

    if (length(layers) > 1L)
        rhdf5::h5createGroup(con, '/layers')
    success <- unlist(Map(
        .exportLoom, assays, name = layers, MoreArgs = list(con = con)
    ))
    if (!all(success == 0L))
        stop(
            '".exportLoom()" failed to write assay(s)\n  ',
            paste0(sQuote(names(layers)[success != 0]), collapse = ', ')
        )

    rhdf5::h5createGroup(con, '/col_attrs')
    rhdf5::h5createGroup(con, '/row_attrs')

    if (is(object, 'SingleCellLoomExperiment')) {
        rdo <- reducedDims(object)
        reducedDims_names <- paste0('/col_attrs/reducedDims_',
            names(rdo))
        lad <- seq_along(reducedDims_names)
        reducedDims_colnames <- paste0(reducedDims_names, "_colnames")
        reducedDims_rownames <- paste0(reducedDims_names, "_rownames")
        reducedDims_attr_names <- paste0('ReducedDimsName', lad)
        if (!length(rdo))
            reducedDims_names <- character(0)
        Map(.exportLoom, rdo, name = reducedDims_names, MoreArgs = list(con = con))

        rdcolnames <- lapply(rdo, colnames)
        rdrownames <- lapply(rdo, rownames)
        if (length(rdcolnames)) {
            reducedDims_attr_colnames <- paste0('ReducedDimsColNames', lad)
            Map(.exportLoom, lapply(rdo, colnames),
                name = reducedDims_colnames, MoreArgs = list(con = con))
        }
        if (any(!is.null(unlist(rdrownames)))) {
            reducedDims_attr_rownames <- paste0('ReducedDimsRowNames', lad)
            Map(.exportLoom, lapply(rdo, rownames),
                name = reducedDims_rownames, MoreArgs = list(con = con))
        }
    }

    .exportLoom(colData(object), con, 'col_attrs', colnames_attr)
    rowData <- rowData(object)
    if (is(object, 'RangedSummarizedExperiment') &&
        !all(lengths(rowRanges <- rowRanges(object)) == 0)) {
            if (is(rowRanges, 'GRangesList')) {
                .exportLoom(rowRanges, con, 'row_attrs', rownames_attr)
            } else {
                .exportLoom(rowRanges, con, 'row_attrs', rownames_attr)
        }
    }
    else
        .exportLoom(rowData, con, 'row_attrs', rownames_attr)

    h5f <- H5Fopen(con)
    tryCatch({
        rhdf5::h5writeAttribute('2.0.1', h5obj=h5f, name='LOOM_SPEC_VERSION')
        rhdf5::h5writeAttribute(paste0('LoomExperiment-', as.character(
            packageVersion('LoomExperiment'))), name='CreatedWith', h5obj=h5f)
        rhdf5::h5writeAttribute(class(object), name='LoomExperiment-class',
            h5obj=h5f)
        rhdf5::h5writeAttribute(matrix, h5obj=h5f, name='MatrixName')
        Map(rhdf5::h5writeAttribute, metadata(object),
            name = names(metadata(object)), MoreArgs = list(h5obj = h5f))
        Map(rhdf5::h5writeAttribute, reducedDims_names,
            name = reducedDims_attr_names, MoreArgs = list(h5obj = h5f))
        if (exists("reducedDims_attr_colnames"))
            Map(rhdf5::h5writeAttribute, reducedDims_colnames,
                name = reducedDims_attr_colnames, MoreArgs = list(h5obj = h5f))
        if (exists("reducedDims_attr_rownames"))
            Map(rhdf5::h5writeAttribute, reducedDims_rownames,
                name = reducedDims_attr_rownames, MoreArgs = list(h5obj = h5f))
    }, error = function(err) {
        warning(conditionMessage(err))
    }, finally = H5Fclose(h5f))

    .exportLoom(colGraphs(object), con, 'col_graphs')
    .exportLoom(rowGraphs(object), con, 'row_graphs')

    invisible(con)
}

#' @importFrom rhdf5 h5createGroup
#' @importFrom rtracklayer export
#' @export
setMethod('export', signature=c('LoomExperiment', 'LoomFile', 'ANY'),
    .exportLoom.LoomExperiment)

