## ----options, include=FALSE, echo=FALSE----------------------------------
#library(BiocStyle)
#knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE)

## ----construct1----------------------------------------------------------
#library(LoomExperiment)
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(assays = list(counts = counts))
scle <- SingleCellLoomExperiment(sce)

## ----construct2----------------------------------------------------------
scle <- SingleCellLoomExperiment(assays = list(counts = counts))

## ----coerce--------------------------------------------------------------
scle <- as(sce, "SingleCellLoomExperiment")
scle

## ----load_l1-------------------------------------------------------------
l1_file <- system.file("extdata", "L1_DRG_20_example.loom", package = "LoomExperiment")
scle <- import(l1_file, type="SingleCellLoomExperiment")
scle

## ----construct_LoomGraph-------------------------------------------------
a <- c(1, 2, 3)
b <- c(3, 2, 1)
w <- c(100, 10, 1)
df <- DataFrame(a, b, w)
lg <- as(df, "LoomGraph")

## OR

lg <- LoomGraph(a, b, weight = w)
lg

## ----subset_LoomGraph----------------------------------------------------
lg[c(1, 2)]
lg[-c(2)]

## ----select_LoomGraph----------------------------------------------------
selectHits(lg, c(1, 3))

## ----drop_LoomGraph------------------------------------------------------
dropHits(lg, c(1, 3))

## ----dropreplace_LoomGraph-----------------------------------------------
dropHits(lg, 2) <- 5
lg

## ----construct_LoomGraphs------------------------------------------------
lgs <- LoomGraphs(lg, lg)
names(lgs) <- c('lg1', 'lg2')
lgs

## ----get_col_row_graphs--------------------------------------------------
colGraphs(scle)
rowGraphs(scle)

## ----replace_LoomGraphs--------------------------------------------------
colGraphs(scle) <- lgs
rowGraphs(scle) <- lgs

colGraphs(scle)
rowGraphs(scle)
colGraphs(scle)[[1]]
rowGraphs(scle)[[1]]

## ----subset_LoomExperiment-----------------------------------------------
scle2 <- scle[c(1, 3), 1:2]
colGraphs(scle2)[[1]]
rowGraphs(scle2)[[1]]

## ----selectdropHits_LoomExperiment---------------------------------------
selectHits(scle, c(1, 3))

dropHits(scle, c(1, 3))

dropHits(scle, c(1, 3)) <- c(7, 8)
colGraphs(scle2)[[1]]
rowGraphs(scle2)[[1]]

## ----export_LoomExperiment-----------------------------------------------
temp <- tempfile(fileext='.loom')
export(scle2, temp)

## ------------------------------------------------------------------------
sessionInfo()

