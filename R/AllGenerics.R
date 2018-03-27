## colGraphs/rowGraphs getters and setters

setGeneric("colGraphs", function(x, ...) standardGeneric("colGraphs"))
setGeneric("colGraphs<-", function(x, ..., value) standardGeneric("colGraphs<-"))

setGeneric("rowGraphs", function(x, ...) standardGeneric("rowGraphs"))
setGeneric("rowGraphs<-", function(x, ..., value) standardGeneric("rowGraphs<-"))


## exportLoom

setGeneric(
    "exportLoom",
    function(object, file = tempfile(), ...) standardGeneric("exportLoom"),
    signature = "object"
)
