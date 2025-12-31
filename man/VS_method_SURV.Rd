\name{VS_method_SURV}

\alias{VS_method_SURV}

\title{
Variable Selection in Leading Sets for Survival Models under the S3VS Framework
}

\description{
\code{VS_method} applies the chosen variable-selection algorithm for survival models to each leading set produced by S3VS at every iteration.
}

\usage{
VS_method_SURV(y, X, surv_model, vsel_method, alpha = 0.5,
  p_thresh = 0.1, gamma = 0.9, verbose = FALSE, ...)
}

\arguments{
  \item{y}{Response. A list with components \code{time} and \code{status} (1 = event, 0 = censored).}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{surv_model}{Character string specifying the survival model. Must be explicitly provided; there is no default. Values are \code{"Cox"} for proportional hazards models, \code{"AFT"} for accelerated failure time models.}
  %
  \item{vsel_method}{Character string indicating the variable-selection engine used at each iteration. Available options are \code{"COXGLMNET"} for \code{surv_model = "COX"} and \code{"AFTREG"}, \code{"AFTGEE"}, \code{"BRIDGE"}, and \code{"PVAFT"} for \code{surv_model = "AFT"}. See Details for more information.}
  %
  \item{alpha}{Only used when \code{vsel_method == "ENET"}. Elastic net mixing parameter, with \eqn{\alpha \in (0,1)}.}
  %
  \item{p_thresh}{Only used with \code{vsel_method = "AFTGEE"}. p-value threshold for variable selection.}
  %
  \item{gamma}{Only used with \code{vsel_method = "BRIDGE"}. Numeric scalar (default \code{0.5}) giving the exponent in the bridge penalty. Must be \eqn{> 0}. See details.}
  %
  \item{verbose}{If \code{TRUE}, some information is printed. \code{FALSE}, by default.}
  %
  \item{...}{Other arguments to be passed inside \code{eha::aftreg()} when \code{surv_model = "AFT"} and \code{vsel_method = "AFTREG"}}
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
\code{\link[glmnet]{cv.glmnet}}, \code{\link[eha]{aftreg}}, \code{\link[aftgee]{aftgee}}, \link{bridge_aft}, \link[afthd]{pvaft}
}

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
# Run VS_method
VS_method_SURV(y_surv, X, surv_model = "COX", vsel_method = "COXGLMNET", verbose = FALSE)
}

