
################################################################################
context("LoomGraph: Constructor and Methods")
################################################################################

test_that("LoomGraph constructor work", {
    ## non-numeric arguments
    expect_error(LoomGraph(a=c('a', 'b', 'c'), b=c('a', 'b', 'c')))
    ## cols not named "a", "b", "w"
    expect_error(LoomGraph(a=c(1, 2, 3), c(3, 2, 1)))
    expect_error(LoomGraph(a=c(1, 2, 3), b=c(3, 2, 1), w1=c(4, 5, 6)))
    ## not between two and three columns
    expect_error(LoomGraph(a=c(1, 2, 3), b=c(5, 4, 3), w=c(5, 4, 3), h=c(5, 4, 3)))

    df <- DataFrame(a=c(1, 2, 3), b=(1, 2, 3), w=c(1, 2, 3))
    lg <- LoomGraph(a=c(1, 2, 3), b=(1, 2, 3), w=c(1, 2, 3))
    lg2 <- LoomGraph(df)

    expect_equal(df, as(lg, "DataFrame"))
    expect_equal(df, as(lg2, "DataFrame"))
})

test_that("LoomGraph methods work", {
    lg <- LoomGraph(a=c(1, 2, 3, 4), b=c(4, 2, 1, 3), w=c(4, 5, 6, 7))
    lg_new <- LoomGraph(a=c(3), b=c(1), w=c(7))
    
    expect_equal(lg[c(1, 3)], lg_new)
})

################################################################################
context("LoomGraphs: Constructor and Methods")
################################################################################

test_that("LoomGraphs constructor works", {
    df <- DataFrame(a=c(1, 2, 3), b=c(1, 2, 3))
    lg <- LoomGraph(df)

    ## error when entry is not a LoomGraph
    expect_error(LoomGraphs(lg, df))

    lgs <- LoomGraphs(lg, lg)
    expect_equal(lg, lgs[[1]])
})
