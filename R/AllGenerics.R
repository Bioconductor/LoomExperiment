## colGraphs/rowGraphs getters and setters

setGeneric('colGraphs', function(x, ...) standardGeneric('colGraphs'))
setGeneric('colGraphs<-', function(x, ..., value) standardGeneric('colGraphs<-'))

setGeneric('rowGraphs', function(x, ...) standardGeneric('rowGraphs'))
setGeneric('rowGraphs<-', function(x, ..., value) standardGeneric('rowGraphs<-'))

setGeneric(
    'loomSelectHits',
    function(x, i, ...) standardGeneric('loomSelectHits')
)

setGeneric('loomDropHits', function(x, i, ...) standardGeneric('loomDropHits'))
setGeneric(
    'loomDropHits<-',
    function(x, i, ..., value) standardGeneric('loomDropHits<-')
)

setGeneric(
    'import.LoomFile',
    function(con, ...) standardGeneric('import.LoomFile')
)

## exportLoom

setGeneric(
    '.exportLoom',
    function(object, con = tempfile(), name, ...) standardGeneric('.exportLoom'),
    signature = 'object'
)
