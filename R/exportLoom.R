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
    #object <- as.data.frame(rowRanges(object))
    object <- as.data.frame(object)
    rhdf5::h5createGroup(con, name)
    .exportLoom(object, con, name, rowname_attr)
})    

setMethod(".exportLoom", "GenomicRangesList",
    function(object, con, name, rowname_attr)
{
#    warning(
#        "'.exportLoom()' does not support '", class(object),
#        "'; using mcols()",
#        call. = FALSE
#    )
#   object <- mcols(object, use.names=TRUE)
    granges_names <- paste0("granges", seq_along(object))
    name <- paste0(name, '/', granges_names)
    Map(.exportLoom, object, name = name, MoreArgs = list(con = con, rowname_attr=rowname_attr))
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
        rhdf5::h5writeAttribute("LoomExperiment", name="CreatedWith", h5obj=h5f)
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

    if (is(object, "SingleCellLoomExperiment")) {
        rhdf5::h5createGroup(con, "/singlecellexperiment")

        reducedDims_names <- "/singlecellexperiment/reducedDims"
        rhdf5::h5createGroup(con, reducedDims_names)
        reducedDims_names <- paste0(reducedDims_names, '/', names(reducedDims(object)))
        Map(.exportLoom, reducedDims(object), name = reducedDims_names, MoreArgs = list(con = con))

        rhdf5::h5createGroup(con, "/singlecellexperiment/int_colData")
        .exportLoom(object@int_colData, con, "singlecellexperiment/int_colData", NULL)

        rhdf5::h5createGroup(con, "/singlecellexperiment/int_elementMetadata")
        .exportLoom(object@int_elementMetadata, con, "singlecellexperiment/int_elementMetadata", NULL)
    }

    rhdf5::h5createGroup(con, "/col_attrs")
    .exportLoom(colData(object), con, "col_attrs", colnames_attr)
    rhdf5::h5createGroup(con, "/row_attrs")
    if (is(object, "RangedSummarizedExperiment"))
        rowData <- GRangesList(rowRanges(object))
    else
        rowData <- rowData(object)
    .exportLoom(rowData, con, "row_attrs", rownames_attr)

    if (length(colGraphs(object)) > 0)
        .exportLoom(colGraphs(object), con, "col_graphs")

    if (length(rowGraphs(object)) > 0)
        .exportLoom(rowGraphs(object), con, "row_graphs")

    invisible(con)
}

#' Method for exporting to .loom files
#'
#' @description
#'  A function for exporting a \code{LoomExperiment} object to a \code{.loom}
#'  con.
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
#' temp_con <- tempcon(conext=".h5")
#' .exportLoom(le, temp_con, rownames_attr="id", colnames_attr="id")
#' rhdf5::h5ls(temp_con)
#' rhdf5::h5dump(temp_con)
#' @export
#' @importFrom rhdf5 h5createGroup
#' @importFrom rtracklayer export
setMethod("export", signature=c("LoomExperiment", "character", "ANY"),
    .exportLoom.LoomExperiment)
