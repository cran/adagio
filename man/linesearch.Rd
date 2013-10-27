\name{linesearch}
\alias{linesearch_ww}
\alias{linesearch_sw}
\title{
  Wolfe Line Search
}
\description{
  Weak and strong Wolfe line search, based on the book by Nocedal and Wright.
}
\usage{
linesearch_ww(x0, d0, fn, gr = NULL, ...,
              dir = c("central", "forward", "backward"),
              c1 = 0, c2 = 0.5, fvalquit = -Inf, trace = 0 )

linesearch_sw(x0, d0, fn, gr = NULL, ...,
              dir = c("central", "forward", "backward"),
              c1 = 0, c2 = 0.5, fvalquit = -Inf, trace = 0)
}
\arguments{
  \item{x0}{initial point.}
  \item{d0}{search direction.}
  \item{fn}{function to be minimized.}
  \item{gr}{gradient function of \code{fn}; will be computed numerically
            if is \code{NULL}; should be passed as column vector.}
  \item{...}{additional parameters to be passed to the function.}
  \item{dir}{direction for computing derivatives; will not be used if 
             \code{gr} is not \code{NULL}.}
  \item{c1}{first Wolfe parameter for the sufficient decrease condition:\cr
            \code{f(x0 + u d) <  f0 + c1*u*t(grad0)\%*\%d}; default 0.}
  \item{c2}{second Wolfe parameter for the weak condition on directional derivative:\cr
            \code{(grad f) t(x0 + u d)\%*\%d > c2*t(grad0)\%*\%d}; default 0.5.}
  \item{fvalquit}{quit immediately if \code{f} drops below this value.}
  \item{trace}{printing level, default 0 for no printing at all.}
}
\details{
  Line search enforcing weak or strong Wolfe conditions, suitable for 
  minimizing both smooth and nonsmooth functions. Strong Wolfe line search 
  with cubic interpolation. For the Wolfe conditions, see the parameter 
  descriptions for \code{c1} and \code{c2} above.

  For functions that are known to be nonsmooth, setting the second Wolfe
  parameter to zero makes sense, especially for a bundle method, and for
  the Shor R-algorithm, for which it is essential.  However, it's not
  a good idea for BFGS, as for smooth functions this may prevent superlinear 
  convergence, and it can even make trouble for BFGS on, e.g., 
  \code{f(x) = x_1^2 + eps |x_2|}, when \code{eps} is small.
}
\value{
  Returns a list with the following components:
  \item{\code{alpha}}{step length satisfying weak Wolfe conditions if one 
        was found, otherwise left end point of interval bracketing such a 
        point (possibly 0)}
  \item{\code{xalpha}}{\code{x0 + alpha*d}.}
  \item{\code{falpha}}{\code{f(x0 + alpha d)}.}
  \item{\code{galpha}}{\code{(grad f)(x0 + alpha d)}.}
  \item{\code{fail}}{0 if both Wolfe conditions satisfied, or falpha < fvalquit\cr
                     1 if one or both Wolfe conditions not satisfied but an 
                     interval was found bracketing a point where both satisfied\cr
                     -1 if no such interval was found, function may be unbounded below.}
  \item{\code{beta} }{same as alpha if it satisfies weak Wolfe conditions,
                      otherwise right end point of interval bracketing such a point
                      (inf if no such finite interval found).}
  \item{\code{gbeta} }{\code{(grad f)(x0 + beta d)} (this is important for 
                          bundle methods), vector of nans if beta is Inf.}
  \item{\code{fvalrec}}{record of function evaluations.}
  \item{\code{nsteps}}{number of steps, strong Wolfe only.}
}
\note{
  The weak Wolfe line search is far less complicated that the standard 
  strong Wolfe line search that is discussed in many texts. It appears
  to have no disadvantages compared to strong Wolfe when used with
  Newton or BFGS methods on smooth functions, and it is essential for 
  the application of BFGS or bundle to nonsmooth functions as done in HANSO.
  However, it is NOT recommended for use with conjugate gradient methods,
  which require a strong Wolfe line search for convergence guarantees.
  Weak Wolfe requires two conditions to be satisfied: sufficient decrease
  in the objective, and sufficient increase in the directional derivative
  (not reduction in its absolute value, as required by strong Wolfe).
}
\references{
  Nocedal, Jorge, and Stephen  Wright (2006). Numerical Optimization.
  Second edition, Springer Science+Business Media, LLC.
}
\author{
  Matlab code by Michael Overton, 2010, under GPL license (HANSO toolbox);
  converted to R by Abhirup Mallik and Hans W. Borchers.
}
\examples{
  fnSphere <- function(x) sum(x * x)
  grSphere <- function(x) 2*x
  x0 <- c(1, 1)
  d0 <- -grSphere(x0)
  linesearch_ww(x0, d0, fnSphere, grSphere)
  # $alpha      0.5
  # $xalpha     0   0
  # $falpha     0
  # $galpha     0   0
  # $fail       0
  # $beta       0.5
  # $gbeta      0   0
  # $fevalrec   2   0

  fnExtSphere <- function(x) sum((1:length(x)) * x * x)
  grExtSphere <- function(x) 2 * (1:length(x)) * x
  x0 <- c(1, 1)
  d0 <- -grExtSphere(x0)
  linesearch_sw(x0, d0, fnExtSphere, grExtSphere)
  # $alpha      0.2777778
  # $xalpha     0.4444444   -0.1111111
  # $falpha     0.2222222
  # $galpha     0.8888889   -0.4444444
  # $fail       0
  # $steps      1
}

