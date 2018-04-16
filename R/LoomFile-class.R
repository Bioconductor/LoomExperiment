setClass("LoomFile", contains="RTLFile")

LoomFile <- function(path)
{
    if (!isSingleString(path))
        stop("'filename' must be a single string, specifiying a path")
    new("LoomFile", resource=path)
}

#setMethod("import.LoomFile", "ANY",
#    function(con, ...)
#{
#    print("hi")
    #import(con, "LoomFile", ...)
#})
