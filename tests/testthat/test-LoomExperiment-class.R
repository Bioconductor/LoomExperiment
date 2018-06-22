library(SummarizedExperiment)
library(SingleCellExperiment)

###############################################################################
context("LoomExperiment: Constructors and Methods")
###############################################################################

ncells <- 100

experiments <- c(RangedLoomExperiment, SingleCellLoomExperiment, SummarizedLoomExperiment)
experiments_names <- c("RangedLoomExperiment", "SingleCellLoomExperiment", "SummarizedLoomExperiment")
base_experiments_con <- c(SummarizedExperiment, SingleCellExperiment, SummarizedExperiment)
base_experiments_class <- c("RangedSummarizedExperiment", "SingleCellExperiment", "SummarizedExperiment")

v <- matrix(rnorm(20000), ncol=ncells)
u <- matrix(rpois(20000, 5), ncol=ncells)
w <- matrix(runif(20000), ncol=ncells)

rd <- DataFrame(stuff=runif(nrow(v)))
cd <- DataFrame(whee=runif(ncells))
sce <- SingleCellExperiment(u, rowData=rd, colData=cd)

rgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2, 3), b=c(2, 1, 3))),
                  KNN=LoomGraph(DataFrame(a=c(3, 2, 1), b=c(3, 1, 2))))
cgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2, 3), b=c(1, 2, 3))),
                  KNN=LoomGraph(DataFrame(a=c(3, 2, 1), b=c(3, 1, 2))))

new_rgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2), b=c(2, 1))),
                  KNN=LoomGraph(DataFrame(a=c(2, 1), b=c(1, 2))))
new_cgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2), b=c(1, 2))),
                  KNN=LoomGraph(DataFrame(a=c(2, 1), b=c(1, 2))))

bad_rgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2, 3), b=c(2, 1, 3))),
                      KNN=LoomGraph(DataFrame(a=c(3, 2, 100000), b=c(1, 3, 2))))
bad_cgs <- LoomGraphs(PPN=LoomGraph(DataFrame(a=c(1, 2, 3), b=c(2, 1, 3))),
                      KNN=LoomGraph(DataFrame(a=c(3, 2, -1), b=c(1, 3, 2))))

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
    for (ex in experiments)
        .test_constructors(ex)
})

.test_coercion <- function(experiment, names, base_experiment, class) {
    le <- experiment(assay=v)
    be <- base_experiment(assay=v)

    le2 <- as(be, names)
    be2 <- as(le, class)

    expect_equivalent(le, le2)
    expect_equivalent(be, be2)
    
    le3 <- as(be, names)
    be3 <- as(le2, class)

    expect_equivalent(le2, le3)
    expect_equivalent(be2, be3)
}

test_that("creation through coercion works", {
    Map(.test_coercion, experiments, experiments_names, base_experiments_con, base_experiments_class)
})

.test_LoomGraphs <- function(experiment) {
    le <- experiment(assay=v, colGraphs=cgs, rowGraphs=rgs)
    
    expect_equivalent(colGraphs(le), cgs)
    expect_equivalent(rowGraphs(le), rgs)

    colGraphs(le) <- rgs
    rowGraphs(le) <- cgs

    expect_equivalent(colGraphs(le), rgs)
    expect_equivalent(rowGraphs(le), cgs)

    colGraphs(le) <- cgs
    rowGraphs(le) <- rgs

    col_le <- le[,c(1,2)]
    expect_equivalent(colGraphs(col_le), new_cgs)
    expect_equivalent(dim(col_le), c(200, 2))

    row_le <- le[c(1,2),]
    expect_equivalent(rowGraphs(row_le), new_rgs)
    expect_equivalent(dim(row_le), c(2, 100))
    
    both_le <- le[c(1,2), c(1,2)]
    expect_equivalent(colGraphs(both_le), new_cgs)
    expect_equivalent(rowGraphs(both_le), new_rgs)
    expect_equivalent(dim(both_le), c(2, 2))

    expect_error(colGraphs(le) <- bad_cgs)
    expect_error(rowGraphs(le) <- bad_rgs)
}

test_that("LoomGraphs work with LoomExperiments", {
    Map(.test_LoomGraphs, experiments)
})
