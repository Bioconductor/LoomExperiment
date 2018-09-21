.convert_subset_index <- function(subset, names) {
    if(is.character(subset)) {
        fmt <- paste0('index out of bounds: %s')
        subset <- SummarizedExperiment:::.SummarizedExperiment.charbound(subset, names, fmt)
    }
    return(as.vector(subset))
}
