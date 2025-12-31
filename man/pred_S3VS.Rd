\name{pred_S3VS}

\alias{pred_S3VS}

\title{
Prediction Using S3VS-Selected Predictors
}

\description{
\code{pred_S3VS} performs prediction using predictors selected by S3VS in linear, generalized linear, and survival models.
}

\usage{
pred_S3VS(y, X, family, surv_model = NULL, method)
}

\arguments{
  \item{y}{Response. If \code{family = "normal"}, a numeric vector. If \code{family = "binomial"}, a numeric/integer/logical vector with values in \{0,1\}. If \code{family = "survival"}, a list with components \code{time} and \code{status} (1 = event, 0 = censored).}
  %
  \item{X}{Predictor matrix. This should include predictors selected by S3VS. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{family}{Model family; one of \code{c("normal","binomial","survival")}. Determines which engine is called (\code{pred_LM}, \code{pred_GLM}, or \code{pred_SURV}).}
  %
  \item{surv_model}{Character string specifying the survival model (\code{family="survival"} only). Must be explicitly provided; there is no default. Values are \code{"Cox"} for proportional hazards models, \code{"AFT"} for accelerated failure time models.}
  %
  \item{method}{Character string indicating the prediction method used. Allowed values depend on \code{family}: for \code{family = "normal"}(functions \code{pred_S3VS_LM}), available options are \code{"NLP"}, \code{"LASSO"}, \code{"SCAD"}, \code{"MCP"}; for \code{"binomial"} (\code{S3VS_GLM}), available options are \code{"NLP"}, \code{"LASSO"}; for \code{family = "survival"} (\code{S3VS_SURV}), available options are \code{"COXGLMNET"} for \code{surv_model = "COX"} and .... for \code{surv_model = "AFT"}. See Details for more information.}
}

%%\details{
%%Details to come...
%%}

\value{
A list containing:
\item{y.pred}{Predicted response}
\item{coef}{Coefficient estimates of the predictors used for prediction}
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
\code{\link{pred_S3VS_LM}}, \code{\link{pred_S3VS_GLM}}, \code{\link{pred_S3VS_SURV}}
}

\examples{
# Simulate continuous data
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
y <- X[,1] + 0.5 * X[,2] + rnorm(n)
# Run S3VS for LM
res_lm <- S3VS(y = y, X = X, family = "normal",
               method_xy = "topk", param_xy = list(k=1),
               method_xx = "topk", param_xx = list(k=3),
               vsel_method = "LASSO", method_sel = "conservative", 
               method_rem = "conservative_begin", rem_regout = FALSE, 
               m = 100, nskip = 3, verbose = TRUE, seed = 123)
pred_lm <- pred_S3VS(y = y, X = X[,res_lm$selected], family = "normal", method = "LASSO")
}

