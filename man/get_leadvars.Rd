\name{get_leadvars}

\alias{get_leadvars}

\title{
Screening Predictors As 'Leading Variables' By Evaluating Predictor-Response Associations 
}

\description{
\code{get_leadvars} screens some predictors as "leading variables" based on predictor-response associations in linear, generalized linear, and survival models.
}

\usage{
get_leadvars(y, X, family = c("normal","binomial","survival"), 
  surv_model = c("AFT", "COX"), 
  method = c("topk", "fixedthresh", "percthresh"), param, 
  varsselected = NULL, varsleft = colnames(X), parallel = FALSE)
}

\arguments{
  \item{y}{Response. If \code{family = "normal"}, a numeric vector. If \code{family = "binomial"}, a numeric/integer/logical vector with values in \{0,1\}. If \code{family = "survival"}, a list with components \code{time} and \code{status} (1 = event, 0 = censored).}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{family}{Model family; one of \code{c("normal","binomial","survival")}. Determines which engine is called (\code{get_leadvars_LM}, \code{get_leadvars_GLM}, or \code{get_leadvars_SURV}).}
  %
  \item{surv_model}{Character string specifying the survival model (\code{family="survival"} only). Must be explicitly provided; there is no default. Values are \code{"Cox"} for proportional hazards models, \code{"AFT"} for accelerated failure time models.}
  %
  \item{method}{Screening rule, one of \code{c("topk", "fixedthresh", "percthresh")}. The association measure depends on \code{family} (e.g., correlation for \code{"normal"}, eta-squared for \code{"binomial"}, or marginal utility for \code{"survival"}). \code{"topk"} keeps the predictors with the largest \eqn{k} association values; \code{"fixedthresh"} keeps predictors whose association is greater than or equal to a specified threshold; \code{"percthresh"} keeps predictors whose association is within a given percentage of the best.}
  %
  \item{param}{Tuning parameter for \code{method}. If \code{"topk"}, supply an integer \eqn{k} (keep the top \eqn{k}). If \code{"fixedthresh"}, supply a numeric threshold (keep predictors with association \eqn{\ge} threshold). If \code{"percthresh"}, supply a percentage in \eqn{(0,100]} (keep predictors with association \eqn{\ge} that percent of the highest association).}
  %
  \item{varsselected}{Used only when \code{family=survival}. A character vector containing the predictors that are already selected in previous iterations. The association measure, conditional utility, is computed controling for these predictors. \code{NULL}, by default.}
  %
  \item{varsleft}{Used only when \code{family=survival}. A character vector containing the predictors that are neither selected, nor removed from consideration in previous iterations. Leading predictors are chosen from these predictors. \code{colnames(X)}, by default.}
  %
  \item{parallel}{Logical. If \code{TRUE}, attempts to perform some computations in parallel mode in \code{binomial} and \code{survival} families, which is strongly recommended for faster execution. Defaults to \code{FALSE}.}
}

%%\details{
%%Details to come...
%%}

\value{A character vector containing the names of the \emph{leading variables}.}

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
\code{\link{get_leadvars_LM}}, \code{\link{get_leadvars_GLM}}, \code{\link{get_leadvars_SURV}}, 
}

\examples{
# Simulate continuous data
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
y <- X[,1] + 0.5 * X[,2] + rnorm(n)
# Select leading variables
leadvars <- get_leadvars(y = y, X = X, family = "normal", 
                         method = "topk", param = list(k=2))
leadvars
}

