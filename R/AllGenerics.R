## colGraphs/rowGraphs getters and setters

setGeneric('colGraphs', function(x, ...) standardGeneric('colGraphs'))
setGeneric('colGraphs<-', function(x, ..., value) standardGeneric('colGraphs<-'))

setGeneric('rowGraphs', function(x, ...) standardGeneric('rowGraphs'))
setGeneric('rowGraphs<-', function(x, ..., value) standardGeneric('rowGraphs<-'))

setGeneric('selectHits', function(x, i, ...) standardGeneric('selectHits'))

setGeneric('dropHits', function(x, i, ...) standardGeneric('dropHits'))
setGeneric('dropHits<-', function(x, i, ..., value) standardGeneric('dropHits<-'))

setGeneric('import.LoomFile', function(con, ...) standardGeneric('import.LoomFile'))

## exportLoom

setGeneric('.exportLoom',
    function(object, con = tempfile(), name, ...) standardGeneric('.exportLoom'),
    signature = 'object'
)
