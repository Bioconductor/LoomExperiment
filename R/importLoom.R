
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
    df <- DataFrame(rhdf5::h5read(con, name))
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
    function(con, name)
{
    ls <- h5ls(con)
    names <- ls[ls$group == name, "name", drop=TRUE]
    gr <- lapply(names, function(x) {
        rhdf5::h5read(con, paste0(name, '/', x))
    })
    names(gr) <- names
    gr <- as.data.frame(gr)
    gr['seqnames'] <- as.character(gr[['seqnames']])
    gr <- GRanges(gr)
    gr
}

.importLoom_GRangesList <-
    function(con, name)
{
    ls <- h5ls(con)
    names <- ls[grep('granges', ls$name), "name", drop=TRUE]
    if (names %in% "granges")
        return(.importLoom_GRanges(con, paste0(name, '/', names)))
    grl <- lapply(names, function(x) {
        .importLoom_GRanges(con, paste0(name, '/', x))
    })
    grl <- GRangesList(grl)
    grl
}

#' Function for importing .loom files
#'
#' @description
#'  a function for importing a \code{.loom} con as a \code{loomexperiment}.
#' @param con character(1) a character vector indicating the file path to the
#'  that is to be imported.
#' @param rownames_attr character indicating the row
#'  attributes in the .loom file that are to be designated as the LoomExperiment
#'  object's rownames.
#' @param rownames_attr character indicating the row
#'  attributes in the .loom con that are to be designated as the LoomExperiment
#'  object's rownames.
#' @return LoomExperiment contained the information from the .loom file.
#' @examples
#' test_con <- system.con(
#'      package="loomexperiment", "extdata", "example.loom"
#' )
#' le <- importLoom(test_con)
#' le
#' @export
#' @importFrom rhdf5 h5ls h5readAttributes
#' @importFrom rtracklayer import
#' @importMethodsFrom rtracklayer import
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

    is_rangedloomexperiment <- nrow(ls[grep('granges', ls$name),]) > 0
    is_singlecellloomexperiment <- nrow(ls[ls$name == "singlecellexperiment",]) > 0

    colData <- .importLoom_DataFrame(con, "col_attrs", colnames_attr)

    if (is_rangedloomexperiment) {
        rowData <- .importLoom_GRangesList(con, "/row_attrs")
    } else {
        rowData <- .importLoom_DataFrame(con, "/row_attrs", rownames_attr)
    }

    if (is_singlecellloomexperiment) {
        reducedDims_names <- "/singlecellexperiment/reducedDims"
        names <- ls[ls$group == reducedDims_names, "name", drop=TRUE]

        reducedDims <- lapply(names, function(x) {
            as.matrix(.importLoom_matrix(con, paste0(reducedDims_names, '/', x)))
        })
        names(reducedDims) <- names
        reducedDims <- SimpleList(reducedDims)

        int_colData <- .importLoom_DataFrame(con, "/singlecellexperiment/int_colData", rownames_attr)
        int_elementMetadata <- .importLoom_DataFrame(con, "/singlecellexperiment/int_elementMetadata", rownames_attr)
        
        int_metadata_fields <- c("version", "spike_names", "size_factor_names")
        int_metadata <- lapply(int_metadata_fields, function(x) {
            int_metadata_names <- paste0("/singlecellexperiment/int_metadata/", x)
            rhdf5::h5read(con, int_metadata_names)
        })
        int_metadata[["version"]] <- package_version(int_metadata[["version"]])
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

    if (is_singlecellloomexperiment) {
        if (length(rowData) == 1)
            rowData <- rowData[[1]]
        le <- SingleCellLoomExperiment(assays, rowData = rowData, colData = colData,
                                       reducedDims = reducedDims,
                                       rowGraphs = row_graphs, colGraphs = col_graphs)
        if (is.null(int_colData))
            int_colData <- DataFrame(matrix(0, nrow(le), 0))
        le@int_colData <- int_colData
        if (is.null(int_elementMetadata))
            int_elementMetadata <- DataFrame(matrix(0, nrow(se), 0))
        le@int_elementMetadata <- int_elementMetadata
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

