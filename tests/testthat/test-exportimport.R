file <- system.file(
    package="LoomExperiment", "extdata", "example.loom"
)

############################################################
context("export_loom and import_loom: Round-trip")
############################################################

test_that("Round-Trip", {
##  Likely not needed
#    temp_file <- tempfile(fileext=".h5")
#    le1 <- import_loom(file)
#    export_loom(le1, temp_file)
#    le2 <- import_loom(temp_file)
#    expect_equal(le1, le2)

#    temp_file2 <- tempfile(fileext=".h5")
#    le1 <- import_loom(file, rownames_attr="id", colnames_attr="id")
#    export_loom(le1, temp_file2)
#    le2 <- import_loom(temp_file2)
#    expect_equal(le1, le2)
})

############################################################
context("import_loom: Import")
############################################################

test_that("Import", {
#    le <- import_loom(file)

    
})

############################################################
context("export_loom: Export")
############################################################

test_that("Export", {
#    le <- import_loom(file)
#    temp <- tempfile(fileext=".h5")
#    export_loom(le, temp)
    
#    ls <- h5ls(file)
#    object <- subset(ls, group)
})
