\name{pred_S3VS_LM}

\alias{pred_S3VS_LM}

\title{
Prediction Using S3VS-Selected Predictors in Linear Models
}

\description{
\code{pred_S3VS} performs prediction using predictors selected by S3VS in linear models.
}

\usage{
pred_S3VS_LM(y, X, method)
}

\arguments{
  \item{y}{Response. A numeric vector.}
  \item{X}{Predictor matrix. This should include predictors selected by S3VS. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  \item{method}{Character string indicating the prediction method used. Available options are \code{"NLP"}, \code{"LASSO"}, \code{"SCAD"}, \code{"MCP"}}
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
\code{\link[mombf]{rnlp}}, \link[glmnet]{cv.glmnet}, \link[ncvreg]{cv.ncvreg}
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
pred_lm <- pred_S3VS_LM(y = y, X = X[,res_lm$selected], method = "LASSO")
pred_lm
}

