#' @importFrom HDF5Array HDF5Array
.import_loom_matrix <-
    function(file, name)
{
    t(HDF5Array::HDF5Array(file, name))
}

#' @importFrom rhdf5 h5read
.import_loom_DataFrame <-
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
#' le <- import_loom(test_file)
#' le
#' @importFrom rhdf5 h5ls h5readAttributes
#' @export
import_loom <-
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

    assay <- .import_loom_matrix(file, "/matrix")
    layerNames <- ls[ls$group == "/layers", "name", drop = TRUE]
    layers <- lapply(setNames(layerNames, layerNames), function(layer) {
        layer <- paste0("/layers/", layer)
        .import_loom_matrix(file, layer)
    })
    assays <- c(list(matrix = assay), layers)

    rowData <- .import_loom_DataFrame(file, "row_attrs", rownames_attr)
    colData <- .import_loom_DataFrame(file, "col_attrs", colnames_attr)

    le <- LoomExperiment(assays, rowData = rowData, colData = colData)
    metadata(le) <- rhdf5::h5readAttributes(file, "/")
    le
}

#' @export
setGeneric(
    "export_loom",
    function(object, file = tempfile(), ...) standardGeneric("export_loom"),
    signature = "object"
)

#' @export
#' @importFrom rhdf5 h5write
setMethod("export_loom", "matrix",
    function(object, file, name)
{
    object <- t(object)
    tryCatch({
        rhdf5::h5write(object, file, name)
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

#' @export
#' @importFrom DelayedArray DelayedArray
#' @importFrom HDF5Array writeHDF5Array
setMethod("export_loom", "DelayedArray",
    function(object, file, name)
{
    HDF5Array::writeHDF5Array(t(object), file, name)
    0L
})

#' @export
#' @importFrom rhdf5 h5write
setMethod("export_loom", "data.frame",
    function(object, file, name, rowname_attr)
{
    if (!is.null(rowname_attr))
        object[[rowname_attr]] <- rownames(object)

    is.factor <- vapply(object, is, logical(1), "factor")
    if (any(is.factor))
        warning(
            "'export_loom()' coerced 'factor' column(s) to character:\n  ",
            paste(sQuote(names(object)[is.factor]), collapse=", "),
            call. = FALSE
        )
    object[is.factor] <- lapply(object[is.factor], as.character)

    names <- sprintf("/%s/%s", name, names(object))
    tryCatch({
        #for (i in seq_along(names))
        #    rhdf5::h5write(object[[i]], names[i], file=file)
        Map(rhdf5::h5write, object, names, MoreArgs = list(file = file))
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

#' @export
setMethod("export_loom", "DataFrame",
    function(object, file, name, rowname_attr)
{
    object <- as.data.frame(object)
    export_loom(object, file, name, rowname_attr)
})

#' @export
#' @import GenomicRanges
setMethod("export_loom", "GenomicRanges",
    function(object, file, name, rowname_attr)
{
    object <- as.data.frame(rowRanges(object))
    export_loom(object, file, name, rowname_attr)
})    

#' @export
setMethod("export_loom", "GenomicRangesList",
    function(object, file, name, rowname_attr)
{
    warning(
        "'export_loom()' does not support '", class(object),
        "'; using mcols()",
        call. = FALSE
    )
    object <- mcols(object, use.names=TRUE)
    export_loom(object, file, name, rowname_attr)
})    

#' Method for exporting to .loom files
#'
#' @description
#'  A function for exporting a \code{LoomExperiment} object to a \code{.loom}
#'  file.
#' @param object LoomExperiment the object that is to be exported.
#' @param file character(1) indicating the file path to the
#'  that is to be imported.
#' @param matrix character(1) indicating the name of the matrix in the .loom
#'  file.
#' @param rownames_attr character indicating the row
#'  attributes in the .loom file that are to be designated as the LoomExperiment
#'  object's rownames.
#' @param rownames_attr character indicating the row
#'  attributes in the .loom file that are to be designated as the LoomExperiment
#'  object's rownames.
#' @return LoomExperiment contained the information from the .loom file.
#' @examples
#' temp_file <- tempfile(fileext=".h5")
#' export_loom(le, temp_file, rownames_attr="id", colnames_attr="id")
#' rhdf5::h5ls(temp_file)
#' rhdf5::h5dump(temp_file)
#' @export
#' @importFrom rhdf5 h5createGroup
setMethod("export_loom", "LoomExperiment",
    function(object, file,
             matrix = assayNames(object)[1],
             rownames_attr = "rownames", colnames_attr = "colnames")
{
    stopifnot(
        !file.exists(file),
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

    rhdf5::h5createFile(file)

    assays <- assays(object, withDimnames = FALSE)
    layers <- setNames(paste0("/layers/", names(assays)), names(assays))
    layers[matrix] <- "/matrix"

    if (length(layers) > 1L)
        rhdf5::h5createGroup(file, "/layers")
    success <- unlist(Map(
        export_loom, assays, layers, MoreArgs = list(file = file)
    ))
    if (!all(success == 0L))
        stop(
            "'export_loom()' failed to write assay(s)\n  ",
            paste0(sQuote(names(layers)[success != 0]), collapse = ", ")
        )

    rhdf5::h5createGroup(file, "/col_attrs")
    export_loom(colData(object), file, "col_attrs", colnames_attr)
    rhdf5::h5createGroup(file, "/row_attrs")
    if (is(object, "RangedSummarizedExperiment"))
        rowData <- rowRanges(object)
    else
        rowData <- rowData(object)
    export_loom(rowData, file, "row_attrs", rownames_attr)

#    rhdf5::h5createGroup(file, "col_graphs")
#    rhdf5::h5createGroup(file, "row_graphs")

    h5f <- H5Fopen(file)
 
    tryCatch({
        Map(rhdf5::h5writeAttribute, metadata(object),
            name = names(metadata(object)), MoreArgs = list(h5obj = h5f))
    }, error = function(err) {
        warning(conditionMessage(err))
    }, finally = H5Fclose(h5f))

    invisible(file)
})

