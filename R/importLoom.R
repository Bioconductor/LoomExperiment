#' @importFrom HDF5Array HDF5Array
.importLoom_matrix <-
    function(file, name)
{
    t(HDF5Array::HDF5Array(file, name))
}

#' @importFrom rhdf5 h5read
.importLoom_DataFrame <-
    function(file, name, rowname)
{
    df <- DataFrame(rhdf5::h5read(file, name))
    df[] <- lapply(df, as.vector)
    if (!is.null(rowname)) {
        rownames(df) <- df[[rowname]]
        df <- df[, -match(rowname, colnames(df)), drop = FALSE]
    }
    if (nrow(df) == 0L)
        df <- NULL
    df
}

#' Function for importing .loom files
#'
#' @description
#'  a function for importing a \code{.loom} file as a \code{loomexperiment}.
#' @param file character(1) a character vector indicating the file path to the
#'  that is to be imported.
#' @param rownames_attr character indicating the row
#'  attributes in the .loom file that are to be designated as the loomexperiment
#'  object's rownames.
#' @param rownames_attr character indicating the row
#'  attributes in the .loom file that are to be designated as the loomexperiment
#'  object's rownames.
#' @return loomexperiment contained the information from the .loom file.
#' @examples
#' test_file <- system.file(
#'      package="loomexperiment", "extdata", "example.loom"
#' )
#' le <- importLoom(test_file)
#' le
#' @importFrom rhdf5 h5ls h5readAttributes
#' @export
importLoom <-
    function(file, rownames_attr = NULL, colnames_attr = NULL)
{
    stopifnot(file.exists(file))

    ls <- rhdf5::h5ls(file)
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

    assay <- .importLoom_matrix(file, "/matrix")
    layerNames <- ls[ls$group == "/layers", "name", drop = TRUE]
    layers <- lapply(setNames(layerNames, layerNames), function(layer) {
        layer <- paste0("/layers/", layer)
        .importLoom_matrix(file, layer)
    })
    assays <- c(list(matrix = assay), layers)

    rowData <- .importLoom_DataFrame(file, "row_attrs", rownames_attr)
    colData <- .importLoom_DataFrame(file, "col_attrs", colnames_attr)

    row_edges <- ls[ls$group == "/row_edges", "name", drop=TRUE]
    col_edges <- ls[ls$group == "/col_edges", "name", drop=TRUE]

    if (length(row_edges) == 0)
        row_edges <- LoomGraphs()
    if (length(col_edges) == 0)
        col_edges <- LoomGraphs()

    if (length(row_edges) > 0) {
        row_edges <- paste0("/row_edges/", row_edges)
        names(row_edges) <- basename(row_edges)

        row_edges <- lapply(row_edges, function(x) {
            res <- as(h5read(file, x), "DataFrame")
            as(res, "LoomGraph")
        })

        row_edges <- .LoomGraphs(row_edges) # as(row_edges, "LoomGraphs")
    }

    if (length(col_edges) > 0) {
        col_edges <- paste0("/col_edges/", col_edges)
        names(col_edges) <- basename(col_edges)

        col_edges <- lapply(col_edges, function(x) {
            res <- as(h5read(file, x), "DataFrame")
            as(res, "LoomGraph")
        })

        col_edges <- .LoomGraphs(col_edges) # as(col_edges, "LoomGraphs")
    }

    le <- LoomExperiment(assays, rowData = rowData, colData = colData,
                         rowGraph = row_edges, colGraph = col_edges)
    metadata(le) <- rhdf5::h5readAttributes(file, "/")
    le
}
