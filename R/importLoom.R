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
setMethod("import.loom", signature=c("ANY"),
    function(con, ..., rownames_attr = NULL, colnames_attr = NULL)
#import_loom <- function(con, rownames_attr = NULL, colnames_attr = NULL)
{
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

    rowData <- .importLoom_DataFrame(con, "row_attrs", rownames_attr)
    colData <- .importLoom_DataFrame(con, "col_attrs", colnames_attr)

    row_graphs <- ls[ls$group == "/row_graphs", "name", drop=TRUE]
    col_graphs <- ls[ls$group == "/col_graphs", "name", drop=TRUE]

    if (length(row_graphs) == 0)
        row_graphs <- LoomGraphs()
    if (length(col_graphs) == 0)
        col_graphs <- LoomGraphs()

    if (length(row_graphs) > 0) {
        row_graphs <- paste0("/row_graphs/", row_graphs)
        names(row_graphs) <- basename(row_graphs)

        row_graphs <- lapply(row_graphs, function(x) {
            LoomGraph(h5read(con, x))
        })

        row_graphs <- LoomGraphs(row_graphs) # as(row_graphs, "LoomGraphs")
    }

    if (length(col_graphs) > 0) {
        col_graphs <- paste0("/col_graphs/", col_graphs)
        names(col_graphs) <- basename(col_graphs)

        col_graphs <- lapply(col_graphs, function(x) {
            LoomGraph(h5read(con, x))
        })

        col_graphs <- LoomGraphs(col_graphs) # as(col_graphs, "LoomGraphs")
    }

    le <- LoomExperiment(assays, rowData = rowData, colData = colData,
                         rowGraphs = row_graphs, colGraphs = col_graphs)
    metadata(le) <- rhdf5::h5readAttributes(con, "/")
    le
})
#}
