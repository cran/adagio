\name{occurs}
\alias{occurs}
\title{
  Finding Subsequences
}
\description{
  Find subsequences of (integer) sequences.
}
\usage{
  occurs(subseq, series)
}
\arguments{
  \item{subseq}{vector of integers.}
  \item{series}{vector of integers.}
}
\details{
  If \code{m} and \code{n} are the lengths of \code{s} and \code{S} resp.,
  \code{occurs(s, S)} determines all positions \code{i} such that
  \code{s == S[i, ..., i+m-1]}.

  The code is vectorized and relatively fast. It is intended to complement
  this with an implementation of Rabin-Karp, and possibly Knuth-Morris-Pratt
  and Boyer-Moore algorithms.
}
\value{
  Returns a vector of indices.
}
\examples{
##  Examples
patrn <- c(1,2,3,4)
exmpl <- c(3,3,4,2,3,1,2,3,4,8,8,23,1,2,3,4,4,34,4,3,2,1,1,2,3,4)
occurs(patrn, exmpl)
## [1]  6 13 23

\dontrun{
set.seed(2437)
p = sample(1:20, 1000000, replace=TRUE)
system.time(i <- occurs(c(1,2,3,4,5), p))  #=>  [1] 799536
##  user  system elapsed 
## 0.017   0.000   0.017 [sec]
}
}
