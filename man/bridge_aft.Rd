\name{bridge_aft}

\alias{bridge_aft}

\title{
Bridge-Penalized AFT Regression via Iteratively Reweighted LASSO
}

\description{
\code{bridge_aft} fits an accelerated failure time (AFT) model using an iterative reweighted LASSO scheme to approximate a bridge (\eqn{L_\gamma}) penalty on the regression coefficients.
}

\usage{
bridge_aft(y, X, gamma = 0.5, alpha = 1, max_iter = 100, tol = 1e-05)
}

\arguments{
  \item{y}{Response; a list of two elements \code{time} and \code{status}.}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{gamma}{Bridge penalty exponent \eqn{\gamma} (default \code{0.5}). Values in \eqn{(0,1]} approximate nonconvex penalties; \eqn{\gamma=1} reduces to LASSO-like weighting.}
  %
  \item{alpha}{Elastic-net mixing parameter passed to \code{glmnet::cv.glmnet} (default \code{1} for LASSO; \code{0} is ridge; \code{(0,1)} is elastic net).}
  %
  \item{max_iter}{Maximum number of outer reweighting iterations (default \code{100}).}
  %
  \item{tol}{Convergence tolerance on the \eqn{L_1} change in coefficients between iterations (default \code{1e-5}).}
}

%%\details{
%%Details to come...
%%}

\value{
A list with components:
  \item{beta}{Numeric vector of estimated coefficients of length \eqn{p}.}
  \item{gamma}{The bridge exponent used in the fit.}
  \item{iterations}{Number of outer reweighting iterations performed.}
}

\references{
Jian Huang and Shuangge Ma. Variable selection in the accelerated failure time model via the bridge method. Lifetime Data Analysis, 16(2):176-195, 2010.
}

\author{
Padmore Prempeh <pprempeh@albany.edu>, Nilotpal Sanyal <nsanyal@utep.edu>
}

%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link[glmnet]{glmnet}}, \code{\link[glmnet]{cv.glmnet}}, \code{\link[survival]{Surv}}, \code{\link[survival]{survfit}}
}

\examples{
set.seed(1)
n  <- 50
p  <- 10
X  <- matrix(rnorm(n * p), n, p)
beta_true <- c(runif(10, -1.5, 1.5), rep(0, p - 10))
linpred   <- as.vector(X \%*\% beta_true)

## Generate log-normal AFT survival times (no censoring in this simple example)
sigma <- 0.6
logT  <- linpred + rnorm(n, sd = sigma)
time  <- exp(logT)
delta <- rep(1, n)  # all events (censoring ignored by current implementation)

y_surv <- list(time = time, status = delta)
fit <- bridge_aft(y_surv, X, gamma = 0.5, alpha = 1, max_iter = 50, tol = 1e-5)
str(fit)
fit$beta[1:10]
}

