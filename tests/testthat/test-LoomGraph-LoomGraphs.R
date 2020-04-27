
################################################################################
context("LoomGraph: Constructor and Methods")
################################################################################

test_that("LoomGraph constructor works", {
    ## non-numeric arguments
    expect_error(LoomGraph(c('a', 'b', 'c'), c('a', 'b', 'c')))
    ## not between two and three columns
    expect_error(LoomGraph(c(1L, 2L, 3L), c(5L, 4L, 3L), weight=c(5, 4, 3), h=c(5, 4, 3)))

    df <- DataFrame(a=c(1, 2, 3), b=c(1, 2, 3), w=c(1, 2, 3))
    lg <- LoomGraph(c(1, 2, 3), c(1, 2, 3), weight=c(1, 2, 3))
    lg2 <- as(df, "LoomGraph")

    expect_equal(df, as(lg, "DataFrame"))
    expect_equal(df, as(lg2, "DataFrame"))
})

test_that("LoomGraph methods work", {
    lg <- LoomGraph(c(1, 2, 3, 4), c(4, 2, 1, 3), weight=c(4, 5, 6, 7))
    lg2 <- LoomGraph(c(7, 2, 8, 4), c(4, 2, 7, 8), weight=c(4, 5, 6, 7))
    lg_new <- LoomGraph(c(3), c(1), weight=c(6))
    
    expect_equivalent(lg[c(3)], lg_new)
    expect_equivalent(lg[-c(1, 2, 4)], lg_new)

    expect_equivalent(loomDropHits(lg, c(2, 4)), lg_new)

    loomDropHits(lg, c(1, 3)) <- c(7, 8)
    expect_equivalent(lg, lg2)
})

test_that("Test that subsetting issue is resolved", {
    counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
    sce <- SingleCellExperiment(assays = list(counts = counts))
    scle <- as(sce, "SingleCellLoomExperiment")

    a <- c(1, 2, 3)
    b <- c(3, 2, 1)
    w <- c(100, 10, 1)
    lg <- LoomGraph(a, b, weight=w)
    lgs <- LoomGraphs(lg, lg)
    names(lgs) <- c('lg1', 'lg2')

    colGraphs(scle) <- lgs
    rowGraphs(scle) <- lgs

    colData(scle)$foo = sample(c("A", "B"), ncol(scle), replace = TRUE)
    expect_error(scle[,scle$foo == "A"], NA)
})

################################################################################
context("LoomGraphs: Constructor and Methods")
################################################################################

test_that("LoomGraphs constructor works", {
    df <- DataFrame(a=c(1, 2, 3), b=c(1, 2, 3))
    lg <- as(df, "LoomGraph")

    ## error when entry is not a LoomGraph
    expect_error(LoomGraphs(lg, df))

    lgs <- LoomGraphs(lg, lg)
   expect_equal(lg, lgs[[1]])
})
