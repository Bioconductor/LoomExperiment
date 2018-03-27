#' @export
#' @importFrom rhdf5 h5write
setMethod("exportLoom", "matrix",
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
setMethod("exportLoom", "DelayedArray",
    function(object, file, name)
{
    HDF5Array::writeHDF5Array(t(object), file, name)
    0L
})

#' @export
#' @importFrom rhdf5 h5write
setMethod("exportLoom", "data.frame",
    function(object, file, name, rowname_attr)
{
    if (!is.null(rowname_attr))
        object[[rowname_attr]] <- rownames(object)

    is.factor <- vapply(object, is, logical(1), "factor")
    if (any(is.factor))
        warning(
            "'exportLoom()' coerced 'factor' column(s) to character:\n  ",
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
setMethod("exportLoom", "DataFrame",
    function(object, file, name, rowname_attr)
{
    object <- as.data.frame(object)
    exportLoom(object, file, name, rowname_attr)
})

#' @export
#' @import GenomicRanges
setMethod("exportLoom", "GenomicRanges",
    function(object, file, name, rowname_attr)
{
    object <- as.data.frame(rowRanges(object))
    exportLoom(object, file, name, rowname_attr)
})    

#' @export
setMethod("exportLoom", "GenomicRangesList",
    function(object, file, name, rowname_attr)
{
    warning(
        "'exportLoom()' does not support '", class(object),
        "'; using mcols()",
        call. = FALSE
    )
    object <- mcols(object, use.names=TRUE)
    exportLoom(object, file, name, rowname_attr)
})    

setMethod("exportLoom", "LoomGraph", 
    function(object, file, name, rowname_attr)
{
    rhdf5::h5createGroup(file, name)
    for (i in names(object))
        rhdf5::h5write(object[[i]], file, paste0(name, '/', i))
})

setMethod("exportLoom", "LoomGraphs", 
    function(object, file, name, rowname_attr)
{
    rhdf5::h5createGroup(file, name)
    name <- paste0(name, names(object))
    for (i in seq_len(length(object)))
        exportLoom(object[[i]], file, name[i], rowname_attr)
})

.exportLoom.Experiment <-
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

    h5f <- H5Fopen(file) 
    tryCatch({
        rhdf5::h5writeAttribute("LoomExperiment", name="createdWith", h5obj=h5f)
        Map(rhdf5::h5writeAttribute, metadata(object),
            name = names(metadata(object)), MoreArgs = list(h5obj = h5f))
    }, error = function(err) {
        warning(conditionMessage(err))
    }, finally = H5Fclose(h5f))

    assays <- assays(object, withDimnames = FALSE)
    layers <- setNames(paste0("/layers/", names(assays)), names(assays))
    layers[matrix] <- "/matrix"

    if (length(layers) > 1L)
        rhdf5::h5createGroup(file, "/layers")
    success <- unlist(Map(
        exportLoom, assays, layers, MoreArgs = list(file = file)
    ))
    if (!all(success == 0L))
        stop(
            "'exportLoom()' failed to write assay(s)\n  ",
            paste0(sQuote(names(layers)[success != 0]), collapse = ", ")
        )

    rhdf5::h5createGroup(file, "/col_attrs")
    exportLoom(colData(object), file, "col_attrs", colnames_attr)
    rhdf5::h5createGroup(file, "/row_attrs")
    if (is(object, "RangedSummarizedExperiment"))
        rowData <- rowRanges(object)
    else
        rowData <- rowData(object)
    exportLoom(rowData, file, "row_attrs", rownames_attr)

    if (length(colGraphs(object)) > 0)
        exportLoom(colGraphs(object), file, "col_graphs", colnames_attr)

    if (length(rowGraphs(object)) > 0)
        exportLoom(rowGraphs(object), file, "row_graphs", rownames_attr)

    invisible(file)
}

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
#' exportLoom(le, temp_file, rownames_attr="id", colnames_attr="id")
#' rhdf5::h5ls(temp_file)
#' rhdf5::h5dump(temp_file)
#' @export
#' @importFrom rhdf5 h5createGroup
setMethod("exportLoom", "LoomExperiment", .exportLoom.Experiment)

setMethod("exportLoom", "RangedLoomExperiment", .exportLoom.Experiment)

setMethod("exportLoom", "SingleCellLoomExperiment", .exportLoom.Experiment)

