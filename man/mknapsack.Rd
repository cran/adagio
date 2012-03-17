\name{mknapsack}
\alias{mknapsack}
\title{
  Multiple 0-1 Knapsack Problem
}
\description{
  Solves the 0-1 (binary) multiple knapsack problem.
}
\usage{
mknapsack(p, w, k, bck = -1)
}
\arguments{
  \item{p}{integer vector of profits.}
  \item{w}{integer vector of weights.}
  \item{k}{integer vector of capacities of the knapsacks.}
  \item{bck}{maximal number of backtrackings allowed; default: -1.}
}
\details{
  Solves the 0-1 multiple knapsack problem for integer profits and weights

  A multiple 0-1 knapsack problem can be formulated as:

  \code{
  maximize  vstar = p(1)*(x(1,1) + ... + x(m,1)) + ...
                        ... + p(n)*(x(1,n) + ... + x(m,n))
  subject to
      w(1)*x(i,1) + ... + w(n)*x(i,n) <= k(i)   for  i=1,...,m
      x(1,j) + ... + x(m,j) <= 1   for  j=1,...,n
      x(i,j) = 0 or 1   for  i=1,...,m ,  j=1,...,n ,
  }

  The input problem description must satisfy the following conditions:
  \itemize{
  \item \code{vs=-1}  if  n < 2  or  m < 1
  \item \code{vs=-2}  if some  p(j) ,  w(j)  or  k(i) are not positive
  \item \code{vs=-3}  if a knapsack cannot contain any item
  \item \code{vs=-4}  if an item cannot fit into any knapsack
  \item \code{vs=-5}  if knapsack  m  contains all the items
  \item \code{vs=-7}  if array  k is not correctly sorted
  \item \code{vs=-8}  [should not happen]
  }
}
\value{
  A list with compomnents, \code{ksack} the knapsack numbers the items are
  assigned to, \code{value} the total value/profit of the solution found, and 
  \code{bs} the number of backtracks used.
}
\note{
  With some care, this function can be used for the bounded and unbounded
  single knapsack problem as well.
}
\author{
  The Fortran source code is adapted from the free NSCW Library of 
  Mathematical Subroutines.

  The wrapping code has been written by yours package maintainer,\cr
  HwB  email: <hwborchers@googlemail.com>
}
\references{
  Kellerer, H., U. Pferschy, and D. Pisinger (2004). Knapsack Problems.
  Springer-Verlag, Berlin Heidelberg.

  Martello, S., and P. Toth (1990). Knapsack Problems: Algorithms and
  Computer Implementations. John Wiley & Sons, Ltd.
}
\seealso{
  Other packages implementing knapsack routines.
}
\examples{
## Example 1: single knapsack
p <- c(15, 100, 90, 60, 40, 15, 10,  1)
w <- c( 2,  20, 20, 30, 40, 30, 60, 10)
cap <- 102
(is <- mknapsack(p, w, cap))
which(is$ksack == 1)
# [1] 1 2 3 4 6 , capacity 102 and total profit 280

## Example 2: multiple knapsack
p <- c(110, 150, 70, 80, 30, 5)
w <- c( 40,  60, 30, 40, 20, 5)
k <- c(65, 85)
is <- mknapsack(p, w, k)
# kps 1: 2,6; kps 2: 1,4; value: 345; backtracks: 14

## Example 3: multiple knapsack
p <- c(78, 35, 89, 36, 94, 75, 74, 79, 80, 16)
w <- c(18,  9, 23, 20, 59, 61, 70, 75, 76, 30)
k <- c(103, 156)
is <- mknapsack(p, w, k)
# kps 1: 1,3,6; kps 2: 4,5,9; value: 452; backtracks: 4 

##  Example 4: subset sum
p <- seq(2, 44, by = 2)^2
w <- p
is <- mknapsack(p, w, 2012)
sum((2 * which(is$ksack == 1))^2)

##  Example 5: maximize number of items
w <- seq(2, 44, by = 2)^2
p <- numeric(22) + 1
is <- mknapsack(p, w, 2012)
}
