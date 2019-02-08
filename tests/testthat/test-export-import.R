
l1_f <- system.file(
    package="LoomExperiment", "extdata", "L1_DRG_20_example.loom"
)

sle_f <- tempfile(fileext=".loom")
rle_empty_f <- tempfile(fileext=".loom")
rle_some_f <- tempfile(fileext=".loom")
rle_full_f <- tempfile(fileext=".loom")
scle_f <- tempfile(fileext=".loom")

file_list <- list(sle_f, rle_empty_f, rle_some_f, rle_full_f, scle_f)

assay <- List(matrix=matrix(seq_len(120), 20))

lg <- LoomGraph(c(1, 2, 3), c(3, 2, 1), weight=c(4, 5, 6))
lgs <- LoomGraphs(PAM=lg, KAM=lg)

granges_one <- GRanges("chr1", IRanges(1, 1))
granges <- rep(granges_one, 20)
grangeslist_empty <- GRangesList(lapply(seq_len(20), function(x) GRanges()))
grangeslist_some <- grangeslist_empty
grangeslist_some[seq(from=4, to=20, by=4)] <- GRangesList(granges)
grangeslist_some[seq(from=3, to=20, by=4)] <- GRangesList(granges_one)
grangeslist_full <- GRangesList(lapply(seq_len(20), function(x) granges))

reducedDims_value <- matrix(seq_len(24), 6)
reducedDims <- List(KNN=reducedDims_value, MKNN=reducedDims_value)

sle <- LoomExperiment(assays=assay, colGraphs=lgs, rowGraphs=lgs)
rle_empty <- RangedLoomExperiment(assays=assay, rowRanges=grangeslist_empty, colGraphs=lgs, rowGraphs=lgs)
rle_some <- RangedLoomExperiment(assays=assay, rowRanges=grangeslist_some, colGraphs=lgs, rowGraphs=lgs)
rle_full <- RangedLoomExperiment(assays=assay, rowRanges=grangeslist_full, colGraphs=lgs, rowGraphs=lgs)
scle <- SingleCellLoomExperiment(assays=assay, rowRanges=granges, reducedDims=reducedDims, colGraphs=lgs, rowGraphs=lgs)

experiment_list <- list(sle, rle_empty, rle_some, rle_full, scle)

################################################################################
context("export method works")
################################################################################

.test_export <- function(experiment, file) {
    export(experiment, file)

    ls <- h5ls(file)

#    expect_true(nrow(ls_rle[ls_rle$group %in% "/row_attrs/GRanges",]) > 0)
#    expect_true(nrow(ls_scle[ls_scle$group %in% "/row_attrs/GRangesList",]) > 0)
}

test_that("export", {
    Map(.test_export, experiment_list, file_list)
})

################################################################################
context("import method works")
################################################################################

test_that("import", {
    rle2 <- import(rle_empty_f)
    expect_equal(dim(rle_empty), dim(rle2))
    rle_m <- as.matrix(assays(rle2)[[1]])
    colnames(rle_m) <- NULL
    rownames(rle_m) <- NULL
    expect_equal(rle_m, assays(rle_empty)[[1]])
    rle_mat0 <- rowData(rle_empty)
    rownames(rle_mat0) <- seq_len(20)
    expect_equal(rle_mat0, rowData(rle2))
    expect_equivalent(colData(rle_empty), colData(rle2))
## Possibly not important to return empty rowRanges?
#    expect_equal(rowRanges(rle_empty), rowRanges(rle2))

    rle2_some <- import(rle_some_f)
    expect_equal(rowRanges(rle_some), rowRanges(rle2_some))

    rle2_full <- import(rle_full_f)
    expect_equal(rowRanges(rle_full), rowRanges(rle2_full))

    sle2 <- import(sle_f)
    expect_equal(dim(sle), dim(sle2))
    sle_m <- as.matrix(assays(sle2)[[1]])
    #colnames(sle_m) <- NULL
    expect_equal(sle_m, as.matrix(assays(sle2)[[1]]))
    expect_equivalent(rowData(sle), rowData(sle2))
    expect_equivalent(colData(sle), colData(sle2))

    scle2 <- import(scle_f)
    expect_equal(dim(scle), dim(scle2))
    scle_m <- as.matrix(assays(sle2)[[1]])
    rownames(scle_m) <- NULL
    expect_equal(scle_m, as.matrix(assays(scle2)[[1]]))
#    expect_equivalent(rowData(scle), rowData(scle2))
    expect_equivalent(colData(scle), colData(scle2))
    expect_equal(rowRanges(scle), rowRanges(scle2))
    rd_scle <- lapply(reducedDims(scle), function(rd) {
        rownames(rd) <- seq_len(nrow(rd))
        rd
    })
    rd_scle <- List(rd_scle)
    expect_equal(rd_scle, reducedDims(scle2))
    expect_equal(scle@int_colData, scle2@int_colData)
    expect_equal(scle@int_elementMetadata, scle2@int_elementMetadata)
    expect_equal(scle@int_metadata, scle2@int_metadata)
    
    l1 <- import(l1_f, type="SingleCellLoomExperiment")
    expect_equal(nrow(l1), 20)
})

