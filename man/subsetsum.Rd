\name{subsetsum}
\alias{subsetsum}
\title{
  Subset Sum Problem
}
\description{
  Subset sum routine for positive integers.
}
\usage{
subsetsum(S, t, method = "greedy")
}
\arguments{
  \item{S}{vector of positive integers.}
  \item{t}{target value.}
  \item{method}{can be ``greedy'' or ``dynamic'', where ``dynamic'' stands
                for the dynamic programming approach.}
}
\details{
  Searching for a set of elements in \code{S} that sum up to \code{t} by
  continuously adding more elements of \code{S}.

  The first components will be preferred, i.e., if \code{S} is decreasing,
  the sum with larger elements will be found, if increasing, the sum with
  smaller elements.

  The dynamic method may be faster for large sets, but will also require
  much more memory if the target value is large.
}
\value{
  List with the target value, if reached, and vector of indices of elements
  in \code{S} that sum up to \code{t}.

  If no solution is found, the dynamic method will return indices for the
  largest value below the target, the greedy method witll return NULL.
}
\note{
  Will be replaced by a compiled version.
}
\author{
  HwB  email: <hwborchers@googlemail.com>
}
\references{
  Horowitz, E., and S. Sahni (1978). Fundamentals of Computer Algorithms.
  Computer Science Press, Rockville, ML.
}
\seealso{
  \code{\link{maxsub}}
}
\examples{
\dontrun{
amount <- 4748652
products <- 
c(30500,30500,30500,30500,42000,42000,42000,42000,
  42000,42000,42000,42000,42000,42000,71040,90900,
  76950,35100,71190,53730,456000,70740,70740,533600,
  83800,59500,27465,28000,28000,28000,28000,28000,
  26140,49600,77000,123289,27000,27000,27000,27000,
  27000,27000,80000,33000,33000,55000,77382,48048,
  51186,40000,35000,21716,63051,15025,15025,15025,
  15025,800000,1110000,59700,25908,829350,1198000,1031655)

# prepare set
prods <- products[products <= amount]  # no elements > amount
prods <- sort(prods, decreasing=TRUE)  # decreasing order

# now find one solution
system.time(is <- subsetsum(prods, amount))
#  user  system elapsed 
# 0.320   0.032   0.359 

prods[is]
#  [1]   70740   70740   71190   76950   77382   80000   83800
#  [8]   90900  456000  533600  829350 1110000 1198000

sum(prods[is]) == amount
# [1] TRUE}
}
\keyword{ optimization }
