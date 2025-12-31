\name{VS_method}

\alias{VS_method}

\title{
Variable Selection in Leading Sets under the S3VS Framework
}

\description{
\code{VS_method} applies the chosen variable-selection algorithm to each leading set produced by S3VS at every iteration.
}

\usage{
VS_method(y, X, family, surv_model = NULL, vsel_method, alpha = 0.5,
          p_thresh = 0.1, gamma = 0.9, verbose = FALSE)
}

\arguments{
  \item{y}{Response. If \code{family = "normal"}, a numeric vector. If \code{family = "binomial"}, a numeric/integer/logical vector with values in \{0,1\}. If \code{family = "survival"}, a list with components \code{time} and \code{status} (1 = event, 0 = censored).}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{family}{Model family; one of \code{c("normal","binomial","survival")}. Determines which engine is called (\code{VS_method_LM}, \code{VS_method_GLM}, or \code{VS_method_SURV}).}
  %
  \item{surv_model}{Character string specifying the survival model (\code{family="survival"} only). Must be explicitly provided; there is no default. Values are \code{"Cox"} for proportional hazards models, \code{"AFT"} for accelerated failure time models.}
  %
  \item{vsel_method}{Character string indicating the variable-selection engine used inside \code{S3VS} at each iteration. Allowed values depend on \code{family}: for \code{family = "normal"}(\code{S3VS_LM}) and \code{"binomial"} (\code{S3VS_GLM}), available options are \code{"NLP"}, \code{"ENET"}, \code{"LASSO"}, \code{"SCAD"}, \code{"MCP"}; for \code{family = "survival"}, available options are \code{"LASSO"}, \code{"ENET"} for \code{surv_model = "COX"} and \code{"AFTGEE"}, \code{"BRIDGE"}, \code{"PVAFT"} for \code{surv_model = "AFT"}. See Details for more information.}
  %
  \item{alpha}{Only used when \code{vsel_method == "ENET"}. Elastic net mixing parameter, with \eqn{\alpha \in (0,1)}.}
  %
  \item{p_thresh}{Only used for \code{surv_model = "AFT"} with \code{vsel_method = "AFTGEE"}. p-value threshold for variable selection.}
  %
  \item{gamma}{Only used for \code{surv_model = "AFT"} with \code{vsel_method = "BRIDGE"}. Numeric scalar (default \code{0.5}) giving the exponent in the bridge penalty. Must be \eqn{> 0}. See details.}
  %
  \item{verbose}{If \code{TRUE}, some information is printed. \code{FALSE}, by default.}
}

\details{
Details to come...
%% Numeric scalar (default \code{0.5}) giving the exponent in the bridge penalty \eqn{\sum_{j=1}^p |\beta_j|^{\gamma}}. It controls the penalty shape and sparsity: \code{gamma = 1} reduces to the lasso, \code{gamma = 2} to ridge, and \code{0 < gamma < 1} gives a concave "bridge" penalty that promotes stronger sparsity with less shrinkage of large effects (but increases non-convexity).
}

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
\code{\link{VS_method_LM}}, \code{\link{VS_method_GLM}}, \code{\link{VS_method_SURV}}
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
VS_method(y, X, family = "normal", vsel_method = "NLP", verbose = FALSE)
}

