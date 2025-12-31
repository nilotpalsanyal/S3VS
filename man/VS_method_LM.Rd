\name{VS_method_LM}

\alias{VS_method_LM}

\title{
Variable Selection in Leading Sets for Linear Models under the S3VS Framework
}

\description{
\code{VS_method} applies the chosen variable-selection algorithm for linear models to each leading set produced by S3VS at every iteration.
}

\usage{
VS_method_LM(y, X, vsel_method, alpha = 0.5, verbose = FALSE)
}

\arguments{
  \item{y}{Response. A numeric vector.}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{vsel_method}{Character string indicating the variable-selection engine used at each iteration. Available options are \code{"NLP"}, \code{"LASSO"}, \code{"SCAD"}, \code{"MCP"}. See Details for more information.}
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
# Simulate continuous data
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
y <- X[,1] + 0.5 * X[,2] + rnorm(n)
# Run VS_method
VS_method_LM(y, X, vsel_method = "NLP", verbose = FALSE)
}

