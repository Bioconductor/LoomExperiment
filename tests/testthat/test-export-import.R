browser()

l1_f <- system.file(
    package="LoomExperiment", "extdata", "L1_DRG_20_example.loom"
)

rle_f <- system.file(
    package="LoomExperiment", "extdata", "rle_example.loom"
)

sle_f <- system.file(
    package="LoomExperiment", "extdata", "sle_example.loom"
)

scle_f <- system.file(
    package="LoomExperiment", "extdata", "scle_example.loom"
)

examples(SummarizedExperiment, echo=FALSE)
examples(SingleCellExperiment, echo=FALSE)

se <- as(rse, "SummarizedExperiment")
sce2 <- as(se, "SummarizedLoomExperiment")
rle2 <- as(rse, "RangedLoomExperiment")
scle2 <- as(scle, "SingleCellLoomExperiment")

################################################################################
context("import method works")
################################################################################

test_that("Import", {
    rle <- import(rle_f)
    expect_equal(dim(rle), dim(rle2))
    expect_equal(assays(rle)[[1]], assays(rle2)[[1]])
    expect_equal(rowData(rle)[[1]], rowData(rle2)[[1]])
    expect_equal(colData(rle), colData(rle2))
    expect_equal(rowRanges(rle), rowRanges(rle2))

    sle <- import(sle_f)
    expect_equal(dim(sle), dim(sle2))
    expect_equal(assays(sle)[[1]], assays(sle2)[[1]])
    expect_equal(rowData(sle)[[1]], rowData(sle2)[[1]])
    expect_equal(colData(sle), colData(sle2))

    scle <- import(scle_f)
    expect_equal(dim(scle), dim(scle2))
    expect_equal(assays(scle)[[1]], assays(scle2)[[1]])
    expect_equal(rowData(scle)[[1]], rowData(scle2)[[1]])
    expect_equal(colData(scle), colData(scle2))
    expect_equal(rowRanges(scle), rowRanges(scle2))
    expect_equal(reducedDims(scle), reducedDims(scle2))
    expect_equal(scle@int_colData, scle2@int_colData)
    expect_equal(scle@int_elementMetadata, scle2@int_elementMetadata)
    expect_equal(scle@int_metadata, scle2@int_metadata)
    
    l1 <- import(l1_f)
    expect_equal(nrow(l1), 282)
})

################################################################################
context("export method works")
################################################################################

test_that("export", {
    sle <- tempfile(fileext=".h5")
    rle <- tempfile(fileext=".h5")
    scle <- tempfile(fileext=".h5")

    export(sle2, sle)
    export(rle2, rle)
    export(scle2, scle)

    ls_sle <- h5ls(sle)
    ls_rle <- h5ls(rle)
    ls_scle <- h5ls(scle)

    expect_true(nrow(ls_rle[ls_rle$group %in% "/row_attrs/GRanges",]) > 0)
    expect_true(nrow(ls_scle[ls_scle$group %in% "/row_attrs/GRangesList",]) > 0)
})

