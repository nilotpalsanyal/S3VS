\name{pred_S3VS_GLM}

\alias{pred_S3VS_GLM}

\title{
Prediction Using S3VS-Selected Predictors in Survival Models
}

\description{
\code{pred_S3VS} performs prediction using predictors selected by S3VS in survival models.
}

\usage{
pred_S3VS_GLM(y, X, method = c("NLP", "LASSO"))
}

\arguments{
  \item{y}{Response. A numeric/integer/logical vector with values in \{0,1\}.}
  \item{X}{Predictor matrix. This should include predictors selected by S3VS. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  \item{method}{Character string indicating the prediction method used. Available options are \code{"NLP"}, \code{"LASSO"}.}
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
\code{\link[mombf]{rnlp}}, \link[glmnet]{cv.glmnet}
}

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
# Predict
pred_glm <- pred_S3VS_GLM(y = y, X = X[,1:3], method = "LASSO")
pred_glm
}

