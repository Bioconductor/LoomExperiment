
#' @importFrom HDF5Array HDF5Array
.importLoom_matrix <-
    function(con, name)
{
    t(HDF5Array::HDF5Array(con, name))
}

#' @importFrom rhdf5 h5read
.importLoom_DataFrame <-
    function(con, name, rowname)
{
    ls <- h5ls(con)
    names_ls <- ls[ls$group %in% name & ls$otype %in% "H5I_DATASET", "name"]
    names <- paste0(name, "/", names_ls)

    df <- lapply(names, function(x) {
        rhdf5::h5read(con, names)
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
    ls <- h5ls(con)
    #gr <- lapply(names, function(x) {
    #    rhdf5::h5read(con, paste0(name, '/', x))
    #})
    gr <- rhdf5::h5read(con, name)
    gr <- gr[[1]]
    gr <- as.data.frame(gr)

    #gr['seqnames'] <- as.character(gr[['seqnames']])
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
            temp <- lapply(grl, `[`, idx+offset, seq_len(len))
            offset <- offset + len
            res <- do.call(rbind, temp)
            res <- t(res)
            colnames(res) <- names(grl)
            rownames(res) <- res[['rownames']]
            res <- subset(res, select = -c(rownames))
            GRanges(res)
        }
    })

    GRangesList(final)
}

#' @importFrom rhdf5 h5ls h5readAttributes
#' @importFrom rtracklayer import
#' @importMethodsFrom rtracklayer import
#' @export
setMethod("import", "loomFile",
    function(con, ..., rownames_attr = NULL, colnames_attr = NULL)
{
    con <- path(con)
    stopifnot(file.exists(con))

    ls <- rhdf5::h5ls(con)
    rowColnames <- ls[ls$group == "/row_attrs", "name", drop=TRUE]
    colColnames <- ls[ls$group == "/col_attrs", "name", drop=TRUE]
    if (missing(rownames_attr) && "rownames" %in% rowColnames)
        rownames_attr <- "rownames"
    if (missing(colnames_attr) && "colnames" %in% colColnames)
        colnames_attr <- "colnames"
    stopifnot(
        is.null(rownames_attr) || rownames_attr %in% rowColnames,
        is.null(colnames_attr) || colnames_attr %in% colColnames
    )

    assay <- .importLoom_matrix(con, "/matrix")
    layerNames <- ls[ls$group == "/layers", "name", drop = TRUE]
    layers <- lapply(setNames(layerNames, layerNames), function(layer) {
        layer <- paste0("/layers/", layer)
        .importLoom_matrix(con, layer)
    })
    assays <- c(list(matrix = assay), layers)

    is_rangedloomexperiment <- nrow(ls[grep('GRanges', ls$name),]) > 0
    is_singlecellloomexperiment <- nrow(ls[ls$name == "reducedDims",]) > 0

    colData <- .importLoom_DataFrame(con, "/col_attrs", colnames_attr)

    if (is_rangedloomexperiment) {
        if (nrow(ls[grep('GRangesList', ls$name),]) > 0)
            rowData <- .importLoom_GRangesList(con, "/row_attrs", ls)
        else
            rowData <- .importLoom_GRanges(con, "/row_attrs", ls)
            
    } else {
        rowData <- .importLoom_DataFrame(con, "/row_attrs", rownames_attr)
    }

    if (is_singlecellloomexperiment) {
        reducedDims_names <- "/col_attrs/reducedDims"
        names <- ls[ls$group %in% reducedDims_names, "name", drop=TRUE]

#        reducedDims <- lapply(names, function(x) {
#            as.matrix(.importLoom_matrix(con, paste0(reducedDims_names, '/', x)))
#        })
        reducedDims <- lapply(paste0(reducedDims_names, '/', names), function(x) {
            #HDF5Array::HDF5Array(con, x)
            as.matrix(.importLoom_matrix(con, x))
        })
        names(reducedDims) <- names
        reducedDims <- SimpleList(reducedDims)

#        int_colData <- .importLoom_DataFrame(con, "/row_attrs/int_colData", rownames_attr)
#        int_elementMetadata <- .importLoom_DataFrame(con, "/row_attrs/int_elementMetadata", rownames_attr)
        
#        int_metadata_fields <- c("version", "spike_names", "size_factor_names")
#        int_metadata <- lapply(int_metadata_fields, function(x) {
#            int_metadata_names <- paste0("/singlecellexperiment/int_metadata/", x)
#            rhdf5::h5read(con, int_metadata_names)
#        })
#        int_metadata[["version"]] <- package_version(int_metadata[["version"]])
    }

    row_graphs <- ls[ls$group == "/row_edges", "name", drop=TRUE]
    if (length(row_graphs) > 0)
        row_graphs_names <- "/row_edges/"
    else {
        row_graphs <- ls[ls$group == "/row_graphs", "name", drop=TRUE]
        row_graphs_names <- "/row_graphs/"
    }
    col_graphs <- ls[ls$group == "/col_edges", "name", drop=TRUE]
    if (length(col_graphs) > 0)
        col_graphs_names <- "/col_edges/"
    else {
        col_graphs <- ls[ls$group == "/col_graphs", "name", drop=TRUE]
        col_graphs_names <- "/col_graphs/"
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

        row_graphs <- LoomGraphs(row_graphs) # as(row_graphs, "LoomGraphs")
    }

    if (length(col_graphs) > 0) {
        col_graphs <- paste0(col_graphs_names, col_graphs)
        names(col_graphs) <- basename(col_graphs)

        col_graphs <- lapply(col_graphs, function(x) {
            LoomGraph(h5read(con, x))
        })

        col_graphs <- LoomGraphs(col_graphs) # as(col_graphs, "LoomGraphs")
    }

    browser()

    if (is_singlecellloomexperiment) {
        le <- SingleCellLoomExperiment(assays, rowData = rowData, colData = colData,
                                       reducedDims = reducedDims,
                                       rowGraphs = row_graphs, colGraphs = col_graphs)
#        if (is.null(int_colData))
#            int_colData <- DataFrame(matrix(0, nrow(le), 0))
#        le@int_colData <- int_colData
#        if (is.null(int_elementMetadata))
#            int_elementMetadata <- DataFrame(matrix(0, nrow(le), 0))
#        le@int_elementMetadata <- int_elementMetadata
    } else if (is_rangedloomexperiment) {
        if (length(rowData) == 1)
            rowData <- rowData[[1]]
        le <- RangedLoomExperiment(assays, rowData = rowData, colData = colData,
                             rowGraphs = row_graphs, colGraphs = col_graphs)
    } else {
        le <- SummarizedLoomExperiment(assays, rowData = rowData, colData = colData,
                             rowGraphs = row_graphs, colGraphs = col_graphs)
    }
    metadata(le) <- rhdf5::h5readAttributes(con, "/")
    le
})

