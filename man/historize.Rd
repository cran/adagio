\name{Historize}
\alias{Historize}
\title{
  Historize function
}
\description{
  Provides storage for function calls.
}
\usage{
Historize(fun, len = 0, ...)
}
\arguments{
  \item{fun}{Function of one or several variables.}
  \item{len}{Allowed length of the input variables.}
  \item{...}{Additional parameters to be handed to the function.}
}
\details{
  Historize() adds storage to the function. If \code{len=0} then only
  function values will be stored in a vector, if \code{len>0} then
  variables and function values will be stored in a matrix, one line
  for every function call.

  If \code{Fn = Historize(...)} is the 'historized' function, then the
  storage can be retrieved with \code{Fn()} and can be cleared with
  \code{Fn(NULL)}.

  Filling the storage will take extra time and can slow down the
  function calls. Especially also storing the variables used in the
  call (with \code{len>0}) will make it considerably slower.

  Functions can have multivariate output; the user is asked to take
  care of handling the output vector or matrix correctly. The function
  may even require additional parameters.
}
\value{
  Returns a function identical to the original function, but the new
  function will now store its function calls in a vector or matrix.
}
\author{
  Hans W. Borchers
}
\note{
  Can also be applied to functions that output a vector (same length
  for every call).
}
\seealso{
  \code{trace}
}
\examples{
Fn <- Historize(fnRastrigin)
dm <- 10
lb <- rep(-5.2, dm); ub <- -lb

# Fn(NULL)
sol <- simpleDE(Fn, lower = lb, upper = ub)
fvalues <- Fn()
fvals <- cummin(fvalues)

\dontrun{
plot(fvalues, type = 'p', col = 7, pch = '.', cex = 2,
     main = "Simple DE optimization", xlab = '', ylab = '')
lines(fvals, col = 2, lwd = 1)
legend(length(fvalues), max(fvalues),
       c("Intermediate values", "cummulated min"),
       xjust = 1, col = c(7, 1), lwd = 2)
grid()}
}
\keyword{ manip }
