
#' @importFrom HDF5Array HDF5Array
.importLoom_matrix <-
    function(con, name)
{
    t(HDF5Array::HDF5Array(con, name))
}

#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom rhdf5 h5read
.importLoom_DataFrame <-
    function(con, name, rowname)
{
    ls <- h5ls(con)
    names_ls <- ls[ls$group %in% name & ls$otype %in% 'H5I_DATASET', 'name']
    names <- paste0(name, '/', names_ls)

    df <- lapply(names, function(x) {
        rhdf5::h5read(con, x)
    })
    names(df) <- names_ls
    df <- DataFrame(df)
    #df <- DataFrame(rhdf5::h5read(con, name))
    df[] <- lapply(df, as.vector)
    if (!is.null(rowname)) {
        rownames(df) <- df[[rowname]]
        df <- df[, -match(rowname, colnames(df)), drop = FALSE]
    }
    if (nrow(df) == 0L)
        df <- NULL
    df
}

.importLoom_GRanges <-
    function(con, name, ls)
{
    name <- paste0(name, '/GRanges')
    ls <- h5ls(con)
    gr <- rhdf5::h5read(con, name)
    gr <- as.data.frame(gr)

    gr['seqnames'] <- as.character(gr[['seqnames']])
    gr['rownames'] <- as.character(gr[['rownames']])
    gr <- subset(gr, select = -c(rownames))

    GRanges(gr)
}

.importLoom_GRangesList <-
    function(con, name, ls)
{
    ls <- h5ls(con)
    grl <- rhdf5::h5read(con, name)
    grl <- grl[[1]]

    lengths <- as.vector(grl[['lengths']])
    names <- as.vector(grl[['names']])
    grl <- grl[!names(grl) %in% c('lengths', 'names')]

    offset <- 0

    final <- lapply(seq_along(lengths), function(idx) {
        len <- lengths[idx]
        if (len == 0) {
            offset <- offset + 1
            GRanges(NULL)
        } else {
            temp <- lapply(grl, `[`, seq_len(len), idx+offset)
            offset <- offset + len
            res <- do.call(cbind.data.frame, temp)
            rownames(res) <- as.character(res[['rownames']])
            group_name <- NULL
            res <- subset(res, select = -c(rownames, group_name))
            GRanges(res)
        }
    })

    GRangesList(final)
}

#' @importFrom rhdf5 h5ls h5readAttributes
#' @importFrom rtracklayer import
#' @importMethodsFrom rtracklayer import
#' @export
setMethod('import', 'LoomFile',
    function(con, ...,
             type = c('SingleCellLoomExperiment', 'LoomExperiment', 'RangedLoomExperiment'),
             rownames_attr = NULL, colnames_attr = NULL)
{
    con <- path(con)
    stopifnot(file.exists(con))

    ls <- rhdf5::h5ls(con)
    rowColnames <- ls[ls$group == '/row_attrs', 'name', drop=TRUE]
    colColnames <- ls[ls$group == '/col_attrs', 'name', drop=TRUE]
    if (missing(rownames_attr) && 'rownames' %in% rowColnames)
        rownames_attr <- 'rownames'
    if (missing(colnames_attr) && 'colnames' %in% colColnames)
        colnames_attr <- 'colnames'
    stopifnot(
        is.null(rownames_attr) || rownames_attr %in% rowColnames,
        is.null(colnames_attr) || colnames_attr %in% colColnames
    )

    metadata <- rhdf5::h5readAttributes(con, '/')

    assay <- .importLoom_matrix(con, '/matrix')
    layerNames <- ls[ls$group == '/layers', 'name', drop = TRUE]
    layers <- lapply(setNames(layerNames, layerNames), function(layer) {
        layer <- paste0('/layers/', layer)
        .importLoom_matrix(con, layer)
    })
    assays <- c(list(matrix = assay), layers)

    is_rangedloomexperiment <- any(grepl('GRanges', ls$name))
    is_singlecellloomexperiment <- nrow(ls[ls$name == 'reducedDims',]) > 0

    colData <- .importLoom_DataFrame(con, '/col_attrs', colnames_attr)

    if (is_rangedloomexperiment) {
        if (nrow(ls[grep('GRangesList', ls$name),]) > 0)
            rowData <- .importLoom_GRangesList(con, '/row_attrs', ls)
        else
            rowData <- .importLoom_GRanges(con, '/row_attrs', ls)
    } else {
        rowData <- .importLoom_DataFrame(con, '/row_attrs', rownames_attr)
    }

    if (is_singlecellloomexperiment) {
        reducedDims_names <- '/col_attrs/reducedDims'
        names <- ls[ls$group %in% reducedDims_names, 'name', drop=TRUE]

        reducedDims <- lapply(paste0(reducedDims_names, '/', names), function(x) {
            as.matrix(.importLoom_matrix(con, x))
        })
        names(reducedDims) <- names
        reducedDims <- SimpleList(reducedDims)
    }

    row_graphs <- ls[ls$group == '/row_edges', 'name', drop=TRUE]
    if (length(row_graphs) > 0)
        row_graphs_names <- '/row_edges/'
    else {
        row_graphs <- ls[ls$group == '/row_graphs', 'name', drop=TRUE]
        row_graphs_names <- '/row_graphs/'
    }
    col_graphs <- ls[ls$group == '/col_edges', 'name', drop=TRUE]
    if (length(col_graphs) > 0)
        col_graphs_names <- '/col_edges/'
    else {
        col_graphs <- ls[ls$group == '/col_graphs', 'name', drop=TRUE]
        col_graphs_names <- '/col_graphs/'
    }

    if (length(row_graphs) == 0)
        row_graphs <- LoomGraphs()
    if (length(col_graphs) == 0)
        col_graphs <- LoomGraphs()

    if (length(row_graphs) > 0) {
        row_graphs <- paste0(row_graphs_names, row_graphs)
        names(row_graphs) <- basename(row_graphs)

        row_graphs <- lapply(row_graphs, function(x) {
            LoomGraph(h5read(con, x))
        })

        row_graphs <- do.call('LoomGraphs', row_graphs)
    }

    if (length(col_graphs) > 0) {
        col_graphs <- paste0(col_graphs_names, col_graphs)
        names(col_graphs) <- basename(col_graphs)

        col_graphs <- lapply(col_graphs, function(x) {
            LoomGraph(h5read(con, x))
        })

        col_graphs <- do.call('LoomGraphs', col_graphs)
    }

    if (!missing(type)) { ## check if LoomExperiment class is specified
        type <- match.arg(type)
        le <- do.call(type, list(assays=assays, rowData=rowData, colData=colData,
                           rowGraphs=row_graphs, colGraphs=col_graphs))
        if (is_singlecellloomexperiment) {
            reducedDims(le) <- reducedDims
        }
    } else { ## discover
        if (is_singlecellloomexperiment) {
            le <- SingleCellLoomExperiment(assays, rowData = rowData, colData = colData,
                                           reducedDims = reducedDims,
                                           rowGraphs = row_graphs, colGraphs = col_graphs)
        } else if (is_rangedloomexperiment) {
            le <- RangedLoomExperiment(assays, rowData = rowData, colData = colData,
                                 rowGraphs = row_graphs, colGraphs = col_graphs)
        } else {
            le <- LoomExperiment(assays, rowData = rowData, colData = colData,
                                 rowGraphs = row_graphs, colGraphs = col_graphs)
        }
    }
    metadata(le) <- metadata
    le
})

