##
##  o c c u r s . R  Finding subsequences
##


occurs <- function(subseq, series){
    # Find all indices i such that series[i, ..., i+m-1] == subseq
    stopifnot(is.numeric(subseq), is.numeric(series))
    m <- length(subseq)
    n <- length(series)
    if (m > n)
        return(as.integer(c()))

    inds <- seq.int(length.out = n-m+1)
    for (i in seq.int(length.out = m))
        inds <- inds[ subseq[i] == series[inds + i - 1] ]
    inds
}
