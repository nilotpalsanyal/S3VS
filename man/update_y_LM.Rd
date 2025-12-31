\name{update_y_LM}

\alias{update_y_LM}

\title{
Update Response Accounting for Selected Predictors in Linear Models
}

\description{
\code{update_y_LM} updates the response accounting for the selected predictors in linear models.
}

\usage{
update_y_LM(y, X, vars)
}

\arguments{
  \item{y}{Response. A numeric vector of length \eqn{n}.}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{vars}{Character vector containing the names of predictors that need to be accounted for. They must appear in \code{X}.}
}

%%\details{
%%Details to come...
%%}

\value{
Returns the updated response vector.
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
update_y(family = "normal", y = y, X = X, vars = c("V1","V4"))
}

