\name{get_leadsets}

\alias{get_leadsets}

\title{
Identify Leading Sets of Covariates via Inter-Predictor Associations
}

\description{
\code{get_leadsets} identifies, for a specified \emph{leading variable}, a set of associated predictors, the leading set, based on inter-predictor associations (absolute value of the correlation coefficient).
}

\usage{
get_leadsets(x_lead, X, method = c("topk", "fixedthresh", "percthresh"), param)
}

\arguments{
  \item{x_lead}{Vector with values of the \emph{leading variable}}
  %
  \item{X}{Predictor matrix. Must contain the \emph{leading variable}. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{method}{Rule for constructing, for each \emph{leading variable}, the set of associated predictors (the "leading set") using inter-predictor association (absolute value of the correlation coefficient); one of \code{c("topk", "fixedthresh", "percthresh")}. \code{"topk"} keeps the predictors with the largest \eqn{k} association values; \code{"fixedthresh"} keeps predictors whose association is greater than or equal to a specified threshold; \code{"percthresh"} keeps predictors whose association is within a given percentage of the best.}
  %
  \item{param}{Tuning parameter for \code{method}; If \code{"topk"}, supply an integer \eqn{k} (keep the top \eqn{k}). If \code{"fixedthresh"}, supply a numeric threshold (keep predictors with association \eqn{\ge} threshold). If \code{"percthresh"}, supply a percentage in \eqn{(0,100]} (keep predictors with association \eqn{\ge} that percent of the highest association).}
}

%%\details{
%%Details to come...
%%}

\value{
A character vector containing the names of the predictors.
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

\examples{
# Simulate continuous data
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
y <- X[,1] + 0.5 * X[,2] + rnorm(n)
leadvars <- get_leadvars_LM(y = y, X = X, method = "topk", param = list(k=2))
get_leadsets(X[,leadvars[1]], X, method = "percthresh", param = list(thresh = 0.2))
}

