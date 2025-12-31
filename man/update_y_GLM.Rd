\name{update_y_GLM}

\alias{update_y_GLM}

\title{
Update Response Accounting for Selected Predictors in Generalized Linear Models
}

\description{
\code{update_y_LM} updates the response accounting for the selected predictors in generalized linear models.
}

\usage{
update_y_GLM(y, X, vars, update_y_thresh)
}

\arguments{
  \item{y}{Response. A numeric/integer/logical vector with values in \{0,1\}.}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{vars}{Character vector containing the names of predictors that need to be accounted for. They must appear in \code{X}.}
  %
  \item{update_y_thresh}{Numeric scalar threshold. When \eqn{|y - fitted\_y| > update\_y\_thresh}, \code{y} is kept, else \code{y} replaced by the rounded value of fitted_y, which is the fitted probability from the logistic model. The default value is 0.5.}
}

%%\details{
%%Details to come...
%%}

\value{
Returns the updated (binary) response vector.
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
# Simulate binary data
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
eta <- X[,1] + 0.5 * X[,2]
prob <- 1 / (1 + exp(-eta))
y <- rbinom(n, size = 1, prob = prob)
update_y(family = "binomial", y = y, X = X, vars = c("V1","V4"), update_y_thresh = 0.8)
}

