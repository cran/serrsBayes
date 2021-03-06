% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sumDexp}
\alias{sumDexp}
\title{Sum log-likelihoods of i.i.d. exponential.}
\usage{
sumDexp(x, rate)
}
\arguments{
\item{x}{Vector of i.i.d. exponential random varibles}

\item{rate}{parameter of the exponential distribution}
}
\value{
log-likelihood of x
}
\description{
This is an internal function that is only exposed on the public API for unit testing purposes.
}
\details{
The sum of the log-likelihoods (log of the product of the likelihoods)
for independent, identically-distributed, exponential random variables. 
Note: this Rcpp function is thread-safe, unlike the equivalent alternatives.
}
\seealso{
\code{sum(dexp(x, rate, log=TRUE))}
}
