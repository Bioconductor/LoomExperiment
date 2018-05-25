
###############################################################################
context("LoomExperiment: Constructors and Methods")
###############################################################################

ncells <- 100

experiments <- c(RangedLoomExperiment, SingleCellLoomExperiment, SummarizedLoomExperiment)
base_experiments <- c(RangedSummarizedExperiment, SingleCellExperiment, SummarizedExperiment)

v <- matrix(rnorm(20000), ncol=ncells)
u <- matrix(rpois(20000, 5), ncol=ncells)
w <- matrix(runif(20000), ncol=ncells)

rd <- DataFrame(stuff=runif(nrow(v)))
cd <- DataFrame(whee=runif(ncells))
sce <- SingleCellExperiment(u, rowData=rd, colData=cd)

rgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2, 3), b=c(2, 1, 3))),
                  KNN=LoomGraph(DataFrame(a=c(3, 2, 1), c=c(3, 1, 2))))
cgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2, 3), b=c(1, 2, 3))),
                  KNN=LoomGraph(DataFrame(a=c(3, 2, 1), c=c(3, 1, 2))))

new_rgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2), b=c(2, 1))),
                  KNN=LoomGraph(DataFrame(a=c(2, 1), c=c(1, 2))))
new_cgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2), b=c(1, 2))),
                  KNN=LoomGraph(DataFrame(a=c(2, 1), c=c(1, 2))))

bad_rgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2, 3), b=c(2, 1, 3))),
                      KNN=LoomGraph(DataFrame(a=c(3, 2, 100000), c=c(1, 3, 2))))
bad_cgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2, 3), b=c(2, 1, 3))),
                      KNN=LoomGraph(DataFrame(a=c(3, 2, -1), c=c(1, 3, 2))))

.test_constructors <- function(experiment) {
    le <- experiment(assay=list(counts=u, exprs=v))
    expect_equivalent(assay(le, "counts"), u)
    expect_equivalent(assay(le, "exprs"), v)

    assay(le, "exprs") <- w
    expect_equivalent(assay(le, "exprs"), w)

    le <- experiment(assay=v)
    expect_equivalent(assay(le), v)
}

test_that("creation through construction works", {
    browser()
    for (ex in experiments)
        .test_constructors(ex)
})

.test_coercion <- function(experiment, base_experiment) {
    le <- experiment(assay=v)
    be <- base_experiment(assay=v)

    le2 <- as(be, as.character(substitute(experiment)))
    be2 <- as(le, as.character(substitute(base_experiment)))

    expect_identical(le, le2)
    expect_identical(be, be2)
    
    le3 <- as(be, as.character(substitute(experiment)))
    be3 <- as(le2, as.character(substitute(base_experiment)))

    expect_identical(le2, le3)
    expect_identical(be2, be3)
}

test_that("creation through coercion works", {
    Map(.test_coercion, experiments, base_experiments)   
})

.test_LoomGraphs <- function(experiment) {
    le <- experiment(assay=v, colGraphs=cgs, rowGraphs=rgs)
    
    expect_identical(colGraphs(le), cgs)
    expect_identical(rowGraphs(le), rgs)

    colGraphs(le) <- rgs
    rowGraphs(le) <- cgs

    expect_identical(colGraphs(le), rgs)
    expect_identical(rowGraphs(le), cgs)

    colGraphs(le) <- cgs
    rowGraphs(le) <- rgs

    col_le <- le[,c(1,2)]
    expect_identical(colGraphs(col_le), cgs_new)
    expect_identical(dims(row_le), c(20000, 2))

    row_le <- le[c(1,2),]
    expect_identical(rowGraphs(row_le), rgs_new)
    expect_identical(dims(row_le), c(2, 100))
    
    both_le <- le[c(1,2), c(1,2)]
    expect_identical(colGraphs(both_le), cgs_new)
    expect_identical(rowGraphs(both_le), rgs_new)
    expect_identical(dims(row_le), c(2, 2))

    expect_error(colGraphs(le) <- bad_cgs)
    expect_error(rowGraphs(le) <- bad_rgs)
}

test_that("LoomGraphs work with LoomExperiments", {
    for (ex in experiments)
        .test_LoomGraphs(ex)
})
