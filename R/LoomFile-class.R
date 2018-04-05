setClass("loomFile", contains="RTLFile")

loomFile <- function(path)
{
    if (!isSingleString(path))
        stop("'filename' must be a single string, specifiying a path")
    new("loomFile", resource=path)
}

setMethod("import.loomFile", "ANY",
    function(con, ...)
{
    print("hi")
    #import(con, "LoomFile", ...)
})
