% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen.data.R
\name{mixtwo}
\alias{mixtwo}
\title{Generate data with correlated and mixed-type predictors}
\usage{
mixtwo(n, sigma, binary)
}
\arguments{
\item{n}{The number of observations.}

\item{sigma}{The error variance.}

\item{binary}{A boolean argument: \code{binary = TRUE} indicates that binary responses are generated and \code{binary = FALSE} 
indicates that continuous responses are generated.}
}
\value{
Return a list with the following components.
\item{X}{An n by p data frame representing predictors values, with each row corresponding an observation.}
\item{Y}{A vector of length n representing response values.}
\item{f0}{A vector of length n representing the values of \eqn{f0(x)}.}
\item{sigma}{The error variance which is only returned when \code{binary = FALSE}.}
\item{prob}{A vector of length n representing the values of \eqn{\Phi(f0(x))}, which is only returned when \code{binary = TRUE}.}
}
\description{
Generate data including responses and predictors values, of which predictors are correlated and of mixed types.
}
\details{
Sample the predictors \eqn{x_1, ..., x_{20}} from Bernoulli(0.2) independently,
\eqn{x_{21}, ..., x_{40}} from Bernoulli(0.5) independently,
and \eqn{x_{41}, ..., x_{84}} from a multivariate normal distribution with mean 0, variance 1 and correlation 0.3.
If \code{binary = FALSE}, sample the continuous response \eqn{y} from Normal(\eqn{f0(x), \sigma^2}), where 
\deqn{f0(x) = -4 + x_1 + sin(\pi x_1*x_{44}) - x_{21} + 0.6x_{41}*x_{42} - exp[-2(x_{42}+1)^2] - x_{43}^2 + 0.5x_{44}.}
If \code{binary = TRUE}, sample the binary response \eqn{y} from Bernoulli(\eqn{\Phi(f0(x))}) where \eqn{f0} is defined above and 
\eqn{\Phi} is the cumulative density function of the standard normal distribution.
}
\examples{
data = mixtwo(100, 1, FALSE)
}
\references{
Luo, C. and Daniels, M. J. (2021)
  "Variable Selection Using Bayesian Additive Regression Trees."
  \emph{arXiv preprint arXiv:2112.13998}.
}
\author{
Chuji Luo: \email{cjluo@ufl.edu} and Michael J. Daniels: \email{daniels@ufl.edu}.
}
