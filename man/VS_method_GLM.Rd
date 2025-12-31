\name{VS_method_GLM}

\alias{VS_method_GLM}

\title{
Variable Selection in Leading Sets for Generalized Linear Models under the S3VS Framework
}

\description{
\code{VS_method} applies the chosen variable-selection algorithm for generalized linear models to each leading set produced by S3VS at every iteration.
}

\usage{
VS_method_GLM(y, X, vsel_method, alpha = 0.5, verbose = FALSE)
}

\arguments{
  \item{y}{Response. A numeric/integer/logical vector with values in \{0,1\}.}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{vsel_method}{Character string indicating the variable-selection engine used at each iteration. Available options are \code{"NLP"} and \code{"LASSO"}. See Details for more information.}
  %
  \item{alpha}{Only used when \code{vsel_method == "ENET"}. Elastic net mixing parameter, with \eqn{\alpha \in (0,1)}.}
  %
  \item{verbose}{If \code{TRUE}, some information is printed. \code{FALSE}, by default.}
}

%%\details{
%%Details to come...
%%}

\value{
A list containing:
\item{sel}{Character vector with names of the selected predictors.}
\item{nosel}{Character vector with names of the predictors not selected.}
}

%%\references{
%%Reference to come...
%%}

\author{
Nilotpal Sanyal <nsanyal@utep.edu>, Padmore N. Prempeh <pprempeh@albany.edu>
}

%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link[mombf]{modelSelection}}, \link[glmnet]{cv.glmnet}, \link[ncvreg]{cv.ncvreg}
}

\examples{
# Simulate binary data
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
eta <- X[,1] + 0.5 * X[,2]
prob <- 1 / (1 + exp(-eta))
y <- rbinom(n, size = 1, prob = prob)
# Run VS_method
VS_method_GLM(y, X, vsel_method = "LASSO", verbose = FALSE)
}

