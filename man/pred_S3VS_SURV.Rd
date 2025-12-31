\name{pred_S3VS_SURV}

\alias{pred_S3VS_SURV}

\title{
Predicted Survival Probabilities Using S3VS-Selected Predictors in Generalized Linear Models
}

\description{
\code{pred_S3VS} returns predicted survival probabilities using predictors selected by S3VS in generalized linear models.
}

\usage{
pred_S3VS_SURV(y, X, surv_model = c("AFT", "COX"), method = c("AFTREG", "AFTGEE"), times)
}

\arguments{
  \item{y}{Response. A list with components \code{time} and \code{status} (1 = event, 0 = censored).}
  %
  \item{X}{Predictor matrix. This should include predictors selected by S3VS. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{surv_model}{Character string specifying the survival model. Must be explicitly provided; there is no default. Values are \code{"Cox"} for proportional hazards models, \code{"AFT"} for accelerated failure time models.}
  %
  \item{method}{Character string indicating the prediction method used. Available options are \code{"COXGLMNET"} for \code{surv_model = "COX"} and \code{"AFTREG"} and \code{"AFTGEE"} for \code{surv_model = "AFT"}. See Details for more information.}
  %
  \item{times}{Vector of time points where predicted survival probabilities will be computed.}
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
\code{\link[glmnet]{cv.glmnet}}, \code{\link[survival]{coxph}}, \code{\link[eha]{aftreg}}, \code{\link[aftgee]{aftgee}}
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
# Run S3VS for linear models
res_surv <- S3VS(y = y_surv, X = X, family = "survival", 
                 surv_model = "COX", vsel_method = "COXGLMNET",
                 method_xy = "topk", param_xy = list(k = 1),
                 method_xx = "topk", param_xx = list(k = 3),
                 method_sel = "conservative", method_rem = "conservative_begin",
                 sel_regout = FALSE, rem_regout = FALSE, 
                 m = 100, nskip = 3, verbose = TRUE, seed = 123)
pred_surv <- pred_S3VS_SURV(y = y_surv, X = X[,res_surv$selected], 
                            surv_model = "COX", method = "COXGLMNET")
pred_surv
}

