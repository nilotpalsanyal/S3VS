\name{get_leadvars_LM}

\alias{get_leadvars_LM}

\title{
Screening Predictors As 'Leading Variables' By Evaluating Predictor-Response Associations In Linear Models 
}

\description{
\code{get_leadvars_LM} screens some predictors as "leading variables" based on predictor-response associations in linear models.}

\usage{
get_leadvars_LM(y, X, method = c("topk", "fixedcorthresh", "perccorthresh"), param)
}

\arguments{
  \item{y}{Response. A numeric vector.}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{method}{Screening rule, one of \code{c("topk", "fixedthresh", "percthresh")}. The association measure is correlation. \code{"topk"} keeps the predictors with the largest \eqn{k} association values; \code{"fixedthresh"} keeps predictors whose association is greater than or equal to a specified threshold; \code{"percthresh"} keeps predictors whose association is within a given percentage of the best.}
  %
  \item{param}{Tuning parameter for \code{method}. If \code{"topk"}, supply an integer \eqn{k} (keep the top \eqn{k}). If \code{"fixedthresh"}, supply a numeric threshold (keep predictors with association \eqn{\ge} threshold). If \code{"percthresh"}, supply a percentage in \eqn{(0,100]} (keep predictors with association \eqn{\ge} that percent of the highest association).).}
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
# Simulate continuous data
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
y <- X[,1] + 0.5 * X[,2] + rnorm(n)
# Select leading variables
leadvars <- get_leadvars_LM(y = y, X = X, method = "topk", param = list(k=2))
leadvars
}

