\name{get_leadvars_GLM}

\alias{get_leadvars_GLM}

\title{
Screening Predictors As 'Leading Variables' By Evaluating Predictor-Response Associations In Generalized Linear Models
}

\description{
\code{get_leadvars_GLM} screens some predictors as "leading variables" based on predictor-response associations in generalized linear models.}

\usage{
get_leadvars_GLM(y, X, method = c("topk", "fixedetasqthresh", "percetasqthresh"), param)
}

\arguments{
  \item{y}{Response. A numeric/integer/logical vector with values in \{0,1\}.}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{method}{Screening rule, one of \code{c("topk", "fixedthresh", "percthresh")}. The association measure is eta-squared. \code{"topk"} keeps the predictors with the largest \eqn{k} association values; \code{"fixedthresh"} keeps predictors whose association is greater than or equal to a specified threshold; \code{"percthresh"} keeps predictors whose association is within a given percentage of the best.}
  %
  \item{param}{Tuning parameter for \code{method}. If \code{"topk"}, supply an integer \eqn{k} (keep the top \eqn{k}). If \code{"fixedthresh"}, supply a numeric threshold (keep predictors with association \eqn{\ge} threshold). If \code{"percthresh"}, supply a percentage in \eqn{(0,100]} (keep predictors with association \eqn{\ge} that percent of the highest association).}
}

%%\details{
%%Details to come...
%%}

\value{A character vector containing the names of the leading varibales.}

%%\references{
%%Reference to come...
%%}

\author{
Nilotpal Sanyal <nsanyal@utep.edu>, Padmore N. Prempeh <pprempeh@albany.edu>
}

%%\note{
%%  ~~further notes~~
%%}

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
# Select leading variables
leadvars <- get_leadvars_GLM(y = y, X = X, method = "topk", param = list(k=2))
leadvars
}

