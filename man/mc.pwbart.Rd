% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc.pwbart.R
\name{mc.pwbart}
\alias{mc.pwbart}
\title{Predicting new observations based on a previously fitted BART model with parallel computation}
\usage{
mc.pwbart(
  x.test,
  treedraws,
  rm.const,
  mu = 0,
  mc.cores = 2L,
  transposed = FALSE,
  dodraws = TRUE,
  nice = 19L
)
}
\arguments{
\item{x.test}{A matrix or a data frame of predictors values for prediction with each row corresponding to an observation 
and each column corresponding to a predictor.}

\item{treedraws}{A list which is the \code{$treedraws} returned from the function \code{wbart()} or \code{pbart()}.}

\item{rm.const}{A vector which is the \code{$rm.const} returned from the function \code{wbart()} or \code{pbart()}.}

\item{mu}{Mean to add on to \code{y} prediction.}

\item{mc.cores}{The number of threads to utilize.}

\item{transposed}{A Boolean argument indicating whether the matrix \code{x.test} is transposed. When
running \code{pwbart()} or \code{mc.pwbart()} in parallel, it is more memory-efficient to transpose \code{x.test} prior to
calling the internal versions of these functions.}

\item{dodraws}{A Boolean argument indicating whether to return the draws themselves (the default), or whether to return the
mean of the draws as specified by \code{dodraws=FALSE}.}

\item{nice}{Set the job niceness. The default niceness is \eqn{19} and niceness goes from \eqn{0} (highest) to \eqn{19} 
(lowest).}
}
\value{
Returns the predictions for \code{x.test}. If \code{dodraws=TRUE}, return a matrix of prediction with each row 
corresponding to a draw and each column corresponding to a new observation; if \code{dodraws=FALSE}, return a vector of 
predictions which are the mean of the draws.
}
\description{
BART is a Bayesian approach to nonparametric function estimation and inference using a sum of trees.\cr
For a continuous response \eqn{y} and a \eqn{p-}dimensional vector of predictors \eqn{x = (x_1, ..., x_p)'}, 
BART models \eqn{y} and \eqn{x} using \deqn{y = f(x) + \epsilon,}
where \eqn{f} is a sum of Bayesian regression trees function and \eqn{\epsilon ~ N(0, \sigma^2)}.\cr
For a binary response \eqn{y}, probit BART models \eqn{y} and \eqn{x} using \deqn{P(Y=1|x)=\Phi[f(x)],}
where \eqn{\Phi} is the CDF of the standard normal distribution and \eqn{f} is a sum of Bayesian regression 
trees function.\cr
The function \code{mc.pwbart()} is inherited from the CRAN R package 'BART'.
}
\examples{
 
## simulate data (Scenario C.M.1. in Luo and Daniels (2021))
set.seed(123)
data = mixone(100, 10, 1, FALSE)
## run wbart() function
res = wbart(data$X, data$Y, ntree=10, nskip=100, ndpost=100)
## parallel::mcparallel/mccollect do not exist on windows
if(.Platform$OS.type=='unix') {
## test pwbart() function
  x.test = mixone(5, 10, 1, FALSE)$X
  pred = mc.pwbart(x.test, res$treedraws, res$rm.const, mu=mean(data$Y), mc.cores=2)
}
}
\references{
Chipman, H. A., George, E. I. and McCulloch, R. E. (2010). 
  "BART: Bayesian additive regression trees."
   \emph{Ann. Appl. Stat.} \strong{4} 266--298.
   
Linero, A. R. (2018). 
  "Bayesian regression trees for high-dimensional prediction and variable selection." 
  \emph{J. Amer. Statist. Assoc.} \strong{113} 626--636.

Luo, C. and Daniels, M. J. (2021)
  "Variable Selection Using Bayesian Additive Regression Trees."
  \emph{arXiv preprint arXiv:2112.13998}.
  
Rockova V, Saha E (2019). 
  “On theory for BART.” 
  \emph{In The 22nd International Conference on Artificial Intelligence and Statistics} (pp. 2839–2848). PMLR.
  
Sparapani, R., Spanbauer, C. and McCulloch, R. (2021).
  "Nonparametric machine learning and efficient computation with bayesian additive regression trees: the BART R package."
  \emph{J. Stat. Softw.} \strong{97} 1--66.
}
\seealso{
\code{\link{wbart}}, \code{\link{pbart}} and \code{\link{mc.pwbart}}.
}
\author{
Chuji Luo: \email{cjluo@ufl.edu} and Michael J. Daniels: \email{daniels@ufl.edu}.
}
