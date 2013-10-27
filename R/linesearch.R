##
##  l i n e s e a r c h . R  Wolfe line search routines (weak and strong)
##


#-- ----------------------------------------------- weak Wolfe line search ---

# The weak Wolfe line search is far less complicated than the standard 
# strong Wolfe line search that is discussed in many texts. It appears
# to have no disadvantages compared to strong Wolfe when used with
# Newton or BFGS methods on smooth functions, and it is essential for 
# the application of BFGS or bundle to nonsmooth functions as done in HANSO.
# Weak Wolfe requires two conditions to be satisfied: sufficient decrease
# in the objective, and sufficient increase in the directional derivative
# (not reduction in its absolute value, as required by strong Wolfe).

linesearch_ww <- function( x0, d0, fn, gr = NULL, ...,
                           dir = c("central", "forward", "backward"),
                           c1 = 0, c2 = 0.5, fvalquit = -Inf, trace = 0 ) {
    stopifnot(is.numeric(x0), is.numeric(d0))
    if (c1 < 0 || c1 > c2 || c2 > 1)    # 0 <= c1 <= c2 <= 1/n
        stop("Arguments 'c1','c2' must satisfy: 0 <= c1 <= c2 <= 1/n.")

    dir <- match.arg(dir)
    fct <- match.fun(fn)
    fn  <- function(x) fct(x, ...)
    if (is.null(gr))
        gr <- function(x) ns.grad(fn, x, dir = dir)

    n <- length(x0)
    d <- d0
    fn0 <- fn(x0)
    gr0 <- as.matrix(gr(x0))

    #-- steplength parameters
    alpha  <- 0         # lower bound on steplength conditions
    xalpha <- x0
    falpha <- fn0
    galpha <- gr0       # need to pass grad0, not grad0'*d, in case line search fails
    beta   <- Inf       # upper bound on steplength satisfying weak Wolfe conditions
    gbeta  <- rep(NA, n)

    g0 <- sum(gr0 * d)
    if (g0 >= 0)
        warning("Linesearch: Argument 'd' is not a descent direction.")
    dnorm <- sqrt(sum(d * d))
    if (dnorm == 0)
        stop("Linesearch: Argument 'd' must have length greater zero.")

    t <- 1              # important to try step length one first
    nfeval   <- 0
    nbisect  <- 0
    nexpand  <- 0
    nbisectmax <- max(30, round(log2(1e5*dnorm)))   # allows more if ||d|| big
    nexpandmax <- max(10, round(log2(1e5/dnorm)))   # allows more if ||d|| small

    #-- main loop ----------------------
    fevalrec <- c()
    fail <- 0
    done <- FALSE
    while (!done) {
        x <- x0 + t*d
        fun <- fn(x)
        grd <- as.matrix(gr(x))
        nfeval <- nfeval + 1
        fevalrec <- c(fevalrec, fun)
        if (fun < fvalquit) {
            return( list(alpha = t, xalpha = x, falpha = fun, galpha = grd,
                        fail = fail, beta = beta, gbeta = gbeta, fevalrec = fevalrec))
        }

        gtd <- sum(grd * d)
        if (fun >= fn0 + c1*t*g0 || is.na(fun)) {# first condition violated
            beta  <- t
            gbeta <- grd
        } else if (gtd <= c2*g0 || is.na(gtd)) {# second condition violated   
            alpha  <- t
            xalpha <- x
            falpha <- fun
            galpha <- grd
        } else {# both conditions satisfied           
            return( list(alpha = t, xalpha = x, falpha = fun, galpha = grd,
                         fail = fail, beta = t, gbeta = grd, fevalrec = fevalrec))
        }

        # set up next function evaluation
        if (beta < Inf) {
            if (nbisect < nbisectmax) {
                nbisect <- nbisect + 1
                t <- (alpha + beta)/2           # bisection
            } else {
                done <- TRUE
            }
        } else {
            if (nexpand < nexpandmax) {
                nexpand <- nexpand + 1
                t <- 2*alpha                    # still in expansion mode
            } else {
                done <- TRUE
            }
        }
    } # end while

    # Wolfe conditions not satisfied; there are two cases:
    if (is.infinite(beta)) {# minimizer never bracketed
        fail <- -1
        if (trace > 1)
            cat("Linesearch: Failed to bracket point satisfying weak\n",
                " Wolfe conditions, function may be unbounded below.\n")
    } else {# point satisfying Wolfe conditions bracketed
        fail <- 1
        if (trace > 1)
            cat("Linesearch: Failed to satisfy weak Wolfe conditions\n",
                " although point satisfying conditions was bracketed.\n")
    }
    
    list(alpha = t, xalpha = x, falpha = fun, galpha = grd,
         fail = fail, beta = t, gbeta = grd, fevalrec = fevalrec)
}


#-- --------------------------------------------- strong Wolfe line search ---

# Based on Nocedal-Wright line search, see their book
# strong Wolfe line search with cubic interpolation
# not recommended for BFGS, especially if function may be nonsmooth: 
#   use the far simpler and preferable linesch_ww instead
# recommended for use with CG, where strong Wolfe condition needed for
# convergence analysis

linesearch_sw <- function( x0, d0, fn, gr = NULL, ...,
                           dir = c("central", "forward", "backward"),
                           c1 = 0, c2 = 0.5, fvalquit = -Inf, trace = 0)
{
    stopifnot(is.numeric(x0), is.numeric(d0))
    if (c1 < 0 || c1 > c2 || c2 > 1)    # 0 <= c1 <= c2 <= 1/n
        stop("Parameters 'c1','c2' are supposed to satisfy: 0 <= c1 <= c2 <= 1.")
    if (fvalquit != -Inf)
        stop("Linesearch: option 'fvalquit' has not yet been implemented.")

    dir <- match.arg(dir)
    fct <- match.fun(fn)
    fn  <- function(x) fct(x, ...)
    if (is.null(gr))
        gr <- function(x) ns.grad(fn, x, dir = dir)

    n <- length(x0)
    d <- d0
    f0 <- fn(x0)
    grad0 <- gr(x0)

    g0 <- sum(grad0*d)
    if (g0 >= 0)
        stop("Linesearch: 'grad * d' must be negative, not a descent direction.")

    dnorm <- sqrt(sum(d * d))
    if (dnorm == 0)
        stop("Linesearch: Argument 'd' must have length greater zero.")

    old <- 0
    fold <- f0
    gold <- g0
    new <- 1

    #-- main loop ----------------------
    nexpand <- max(50, round(log2(dnorm)))
    for(k in 1:nexpand) {
        xnew <- x0 + new*d
        fnew <- fn(xnew)
        gradnew <- gr(xnew)
        gnew <- sum(gradnew*d)

        if (fnew > f0 + c1*new*g0 || ((fnew >= fold) && k > 1)) {  # gone too far
            lsz <- lszoom(fn, gr, old, new, fold, fnew, gold, gnew, f0, g0, x0,
                           d, c1, c2, trace)
            return(list(alpha = lsz$alpha, xalpha = lsz$x, falpha = lsz$f, galpha = lsz$grd,
                        fail = lsz$fail, nsteps = lsz$nsteps))
        }
      
        if (abs(gnew) <= -c2*g0) {  # quit, conditions are satisfied first time
            return(list(alpha = new, xalpha = xnew, falpha = fnew, galpha = gradnew,
                        fail = 0, nsteps = k))
        }

        if (gnew >= 0) {  # gone past minimizer, but new point is better than old
            lsz <- lszoom(fn, gr, new, old, fnew, fold, gnew, gold, f0, g0, x0,
                           d, c1, c2, trace)
            return(list(alpha = lsz$alpha, xalpha = lsz$x, falpha = lsz$f, galpha = lsz$grd,
                        fail = lsz$fail, nsteps = lsz$nsteps))
        }

        # minimizer not bracketed yet
        old <- new
        fold <- fnew
        gold <- gnew
        new <- 2*new
    }  # end for

    if (trace > 1)
        cat("Linesearch: minimizer was not bracketed, function may be unbounded below.\n")
    alpha <- new
    x <- xnew
    f <- fnew
    grd <- gnew
    fail <- -1
    nsteps <- 0  # lszoom was never called
    
    list(alpha = alpha, xalpha = x, falpha = f, galpha = grd,
                fail = fail, nsteps = nsteps)
}


#-- ----------------------------------------------------- line search zoom ---

# Based on Nocedal Wright line search Zoom subroutine, see their book,
# does accurate line search efficiently if 2nd Wolfe parameter c2 is small
# strong Wolfe line search with cubic interpolation, does accurate line
# search efficiently if 2nd Wolfe parameter c2 is small
# intended only to be called by linesch_sw.m

lszoom <- function(fn, gr, lo, hi, flo, fhi, glo, ghi, f0, g0, x0,
                    d, c1, c2, trace)
{
    # Initialization
    fail   <- 0
    lo2    <- hi
    flo2   <- fhi
    glo2   <- ghi
    nsteps <- 0
  
    while ((lo != hi) && (nsteps < 50)) {
        nsteps <- nsteps + 1
        bisect <- (lo + hi)/2
        interp <- nso.cubic_interp(lo, lo2, flo, flo2, glo, glo2)

        if (nso.inside(interp, lo, bisect)) atry <- interp
        else atry <- bisect
      
        xtry <- x0 + atry * d
        ftry <- fn(xtry)
        gradtry <- gr(xtry)
        gtry <- sum(gradtry*d)
      
        if (ftry > f0 + c1*atry*g0 || ftry >= flo){
            hi <- atry
            lo2 <- hi
            flo2 <- ftry
            glo2 <- gtry
        } else {
            if(abs(gtry) <= -c2*g0) {  # strong Wolfe conditions are satisfied
                  alpha <- atry
                  x <- xtry
                  f <- ftry
                  grd <- gradtry
                  return(list(alpha = alpha, x = x, f = f, grd = grd,
                              fail = fail, nsteps = nsteps))
            }
            if (gtry*(hi - lo) >= 0) {
                  hi <- lo
            }
            lo2 <- lo
            flo2 <- flo
            glo2 <- glo
            lo <- atry
            flo <- ftry
            glo <- gtry
        }
    }  # end while
  
    if (trace > 1)
        cat("lszoom: failed to satisfy wolfe conditions; loop ran for 50 times.\n")

    alpha <- atry
    x <- xtry
    f <- ftry
    grd <- gradtry
    fail <- 1 
    list(alpha = alpha, x = x, f = f, grd = grd,
         fail = fail, nsteps = nsteps)
}


#-- ---------------------------------------------------- Utility functions ---

nso.cubic_interp <- function(x1, x2, f1, f2, g1, g2) {
  eta   <- g1 + g2 - 3*(f1 - f2)/(x1 - x2) 
  gamma <- sign(x2 - x1) * sqrt(eta^2 - g1 * g2)
  x2 - (x2 - x1) * (g2 + gamma - eta) / (g2 - g1 + 2*gamma)
}

nso.inside <- function(x, a, b) {
    in_p <- FALSE
    if (is.na(x)) return(in_p)
    if (a <= b) {
        if (x >= a & x <= b) in_p <- TRUE
    } else {
        if (x >= b & x <= a) in_p <- TRUE
    }
    in_p  
}

#--- EoF ---
