##
##  Simple Evolutionary Algorithm
##


simpleEA <-
function(fun, lower, upper, N=100, 
                     eps=1e-6, scl=1/2, log=TRUE)
{
    # if ( length(lower) != length(upper) ) error()
    # if ( !scale < 1 ) error

    n <- length(lower)
    z <- (upper + lower)/2      # midpoint of rectangle
    h0 <- upper - lower         # side lengths of rectangle

    P <- 1                      # no. of parents
    Parents <- matrix(z, nrow=P)
    fvals   <- fun(z)

    h <- scl * h0
    while ( min(h) > eps ) {
        newParents <- rbind(    # TODO: keep new parents in rectangle !
            Parents,
            kronecker(Parents, rep(1, N)) + 
                matrix(runif(P*N*n, -1, 1), ncol=n) %*% diag(h) )
        newFvals <- apply(newParents, 1, fun)
        orderParents <- order(newFvals)
        Parents <- newParents[orderParents[1:N],]
        fvals <- apply(Parents, 1, fun)

        if (log) cat(min(fvals), "\n")

        # initialize new loop
        P <- N                  # TODO: avoid this
        h <- scl * h
    }

    return( list(par=Parents[1,], val=min(fvals), scl=max(h)) )
}

