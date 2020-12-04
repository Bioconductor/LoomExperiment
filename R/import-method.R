
.importLoom_colchar <-
    function(con, name)
{
    name <- as.character(name)
    as.character(rhdf5::h5read(con, name))
}

#' @importFrom HDF5Array HDF5Array
.importLoom_matrix <-
    function(con, name)
{
    name <- as.character(name)
    t(HDF5Array::HDF5Array(con, name))
}

#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom rhdf5 h5read
.importLoom_DataFrame <-
    function(con, name, rowname, exclude)
{
    ls <- h5ls(con)
    names_ls <- ls[ls$group %in% name & ls$otype %in% 'H5I_DATASET', 'name']
    names <- paste0(name, '/', names_ls)
    cfa <- paste0(name, "/", "colnames_factor")
    has_factor <- cfa %in% names
    if (!missing(exclude)) {
        indices <- !names %in% c(exclude, cfa)
        names <- names[indices]
        names_ls <- names_ls[indices]
    }

    if (has_factor)
        isfactor <- as.logical(t(rhdf5::h5read(con, cfa))[1L, ])
    df <- lapply(names, function(x) {
        rhdf5::h5read(con, x)
    })
    names(df) <- names_ls
    if (has_factor) {
        df <- Map(function(x, y) { if (y) as.factor(x) else as.vector(x) },
            x = df, y = isfactor)
    } else
        df[] <- lapply(df, as.vector)
    df <- DataFrame(df)

    if (!is.null(rowname)) {
        rnames <- as.character(df[[rowname]])
        hasdfrownames <- identical(rnames, as.character(seq_along(rnames)))
        if (!hasdfrownames) {
            rownames(df) <- df[[rowname]]
        }
        df <- df[, -match(rowname, colnames(df)), drop = FALSE]
    }
    if (nrow(df) == 0L)
        df <- NULL
    df
}

#' @importFrom stringr str_replace
.importLoom_GRanges <-
    function(con, name, ls)
{
    #name <- paste0(name, '/GRanges')
    ls <- h5ls(con)
    labels <- ls$name
    labels <- labels[grep("GRanges_", labels)]
    names <- stringr::str_replace(labels, "GRanges_", "")
    labels <- paste0(name, '/', labels)

    gr <- Map(rhdf5::h5read, con, labels)
    names(gr) <- names
    gr <- as.data.frame(gr)
    gr['seqnames'] <- as.character(gr[['seqnames']])
    gr['rownames'] <- as.character(gr[['rownames']])
    gr <- subset(gr, select = -c(rownames))

    gr <- GRanges(gr)
    GRanges(gr)
}

.importLoom_GRangesList <-
    function(con, name, ls)
{
    ls <- h5ls(con)
    labels <- ls$name
    labels <- labels[grep("GRangesList_", labels)]
    names <- stringr::str_replace(labels, "GRangesList_", "")
    labels <- paste0(name, '/', labels)

    grl <- Map(rhdf5::h5read, con, labels)
    names(grl) <- names

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

.eq_dims <- function(asys, coldat , rowdat) {
    adims <- if (all(isEmpty(asys))) c(0L, 0L) else dim(asys[[1L]])
    zrows <- if (isEmpty(rowdat)) {
        adims[1]
    } else if (inherits(rowdat, "GenomicRanges_OR_GRangesList")) {
        length(rowdat)
    } else {
        nrow(rowdat)
    }
    crows <- if (isEmpty(coldat)) adims[2] else nrow(coldat)
    identical(adims, c(zrows, crows))
}

#' @importFrom rhdf5 h5ls h5readAttributes
#' @importFrom BiocIO import
#' @importFrom stringr str_split
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
    metadata_names <- names(metadata)
    idx <- grep("ReducedDims", metadata_names)
    reducedDims_names <- metadata_names[idx]

    assay <- .importLoom_matrix(con, '/matrix')
    layerNames <- ls[ls$group == '/layers', 'name', drop = TRUE]
    layers <- lapply(setNames(layerNames, layerNames), function(layer) {
        layer <- paste0('/layers/', layer)
        .importLoom_matrix(con, layer)
    })

    assay_matrix <- list(assay)
    matrix_name <- metadata$MatrixName
    if(is.null(matrix_name))
       matrix_name <- 'matrix'
    names(assay_matrix) <- matrix_name
    assays <- c(assay_matrix, layers)

    is_rangedloomexperiment <- any(grepl('GRanges', ls$name))
    is_singlecellloomexperiment <- any(grepl('reducedDims', ls$name))

    colData <- .importLoom_DataFrame(con, '/col_attrs', colnames_attr,
        unlist(metadata[reducedDims_names]))

    if (is_rangedloomexperiment) {
        if (any(grepl('GRangesList', ls$name)))
            rowData <- .importLoom_GRangesList(con, '/row_attrs', ls)
        else
            rowData <- .importLoom_GRanges(con, '/row_attrs', ls)
    } else {
        rowData <- .importLoom_DataFrame(con, '/row_attrs', rownames_attr)
    }

    if (is_singlecellloomexperiment) {
        if (length(reducedDims_names) == 0)
            reducedDims <- list()
        else {
            reducedDims <- grep("colnames", metadata[reducedDims_names],
                invert = TRUE, value = TRUE)
            reducedattrs <- strsplit(unlist(reducedDims), "_")
            names <- vapply(reducedattrs, `[[`, character(1L), 3L)
            rdimcols <- metadata[grep("ReducedDimsColNames", metadata_names)]
            reducedDims <- lapply(reducedDims, function(x)
                as.matrix(.importLoom_matrix(con, x))
            )
            if (length(rdimcols)) {
                withCols <- gsub("ColNames", "Name", names(rdimcols))
                reducedDims[withCols] <- Map(function(x, y) {
                    colnames(x) <- .importLoom_colchar(con, y)[seq_len(ncol(x))]
                    x
                },
                x = reducedDims[withCols],
                y = rdimcols)
            }
            names(reducedDims) <- names
            reducedDims <- SimpleList(reducedDims)
        }
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
            df <- DataFrame(h5read(con, x))
            df <- as.matrix(df)
            df[,1] <- df[,1] + 1L
            df[,2] <- df[,2] + 1L
            df <- DataFrame(df)
            as(df, "LoomGraph")
        })

        row_graphs <- do.call('LoomGraphs', row_graphs)
    }

    if (length(col_graphs) > 0) {
        col_graphs <- paste0(col_graphs_names, col_graphs)
        names(col_graphs) <- basename(col_graphs)

        col_graphs <- lapply(col_graphs, function(x) {
            df <- DataFrame(h5read(con, x))
            df <- as.matrix(df)
            df[,1] <- df[,1] + 1L
            df[,2] <- df[,2] + 1L
            df <- DataFrame(df)
            as(df, "LoomGraph")
        })

        col_graphs <- do.call('LoomGraphs', col_graphs)
    }

    ## if a transposition couldn't be performed on the array in export; do it here
    if (!.eq_dims(assays, colData, rowData))
        assays <- lapply(assays, t)

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

    metadata <- metadata[-c(idx)]
    metadata(le) <- metadata
    le
})

