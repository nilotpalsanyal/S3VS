\name{get_leadvars_SURV}

\alias{get_leadvars_SURV}

\title{
Screening Predictors As "Leading Variables" By Evaluating Predictor-Response Associations In Survival Models 
}

\description{
\code{get_leadvars_SURV} screens some predictors as "leading variables" based on predictor-response associations in survival models.}

\usage{
get_leadvars_SURV(y, X, surv_model = c("AFT", "COX"), 
  method = c("topk", "fixedmuthresh", "percmuthresh"), param, 
  varsselected = NULL, varsleft = colnames(X), parallel = FALSE)
}

\arguments{
  \item{y}{Response. A list with components \code{time} and \code{status} (1 = event, 0 = censored).}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{surv_model}{Character string specifying the survival model. Must be explicitly provided; there is no default. Values are \code{"Cox"} for proportional hazards models, \code{"AFT"} for accelerated failure time models.}
  %
  \item{method}{Screening rule, one of \code{c("topk", "fixedthresh", "percthresh")}. The association measure is marginal utility. \code{"topk"} keeps the predictors with the largest \eqn{k} association values; \code{"fixedthresh"} keeps predictors whose association is greater than or equal to a specified threshold; \code{"percthresh"} keeps predictors whose association is within a given percentage of the best.}
  %
  \item{param}{Tuning parameter for \code{method}. If \code{"topk"}, supply an integer \eqn{k} (keep the top \eqn{k}). If \code{"fixedthresh"}, supply a numeric threshold (keep predictors with association \eqn{\ge} threshold). If \code{"percthresh"}, supply a percentage in \eqn{(0,100]} (keep predictors with association \eqn{\ge} that percent of the highest association).}
  %
  \item{varsselected}{A character vector containing the predictors that are already selected in previous iterations. The association measure, conditional utility, is computed controling for these predictors. \code{NULL}, by default.}
  %
  \item{varsleft}{A character vector containing the predictors that are neither selected, nor removed from consideration in previous iterations. Leading predictors are chosen from these predictors.}
  %
  \item{parallel}{Logical. If \code{TRUE}, attempts to perform some computations in parallel mode in \code{binomial} and \code{survival} families, which is strongly recommended for faster execution. Defaults to \code{colnames(X)}.}
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
# Simulate survival data (Cox)
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
eta <- X[,1] + 0.5 * X[,2]
base_rate <- 0.05
T_event <- rexp(n, rate = base_rate * exp(eta))
C <- rexp(n, rate = 0.03)
time <- pmin(T_event, C)
status <- as.integer(T_event <= C)
y_surv <- list(time = time, status = status)
# Select leading variables
leadvars <- get_leadvars_SURV(y = y_surv, X = X, surv_model = "COX", 
                              method = "topk", param = list(k=2), 
                              varsselected = NULL, varsleft = colnames(X))
leadvars
}

