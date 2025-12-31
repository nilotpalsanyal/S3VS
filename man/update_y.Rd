\name{update_y}

\alias{update_y}

\title{
Update Response Accounting for Selected Predictors
}

\description{
\code{update_y} updates the response accounting for the selected predictors in linear models, and selected or removed predictors in generalized linear models.
}

\usage{
update_y(y, X, family, vars, update_y_thresh = NULL)
}

\arguments{
  \item{y}{Response. If \code{family = "normal"}, a numeric vector. If \code{family = "binomial"}, a numeric/integer/logical vector with values in \{0,1\}.}
  %
  \item{family}{Model family; one of \code{c("normal","binomial")}. Determines which engine is called (\code{update_y_LM}, \code{update_y_GLM}).}
  %
  \item{X}{Predictor matrix. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{vars}{Character vector containing the names of predictors that need to be accounted for. They must appear in \code{X}.}
  %
  \item{update_y_thresh}{Numeric scalar threshold used only \code{family = "binomial"}. When \eqn{|y - fitted\_y| > update\_y\_thresh}, \code{y} is kept, else \code{y} replaced by the rounded value of \eqn{fitted_y}, which is the fitted probability from the logistic model. The default value is 0.5.}
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

\seealso{
\code{\link{update_y_LM}}, \code{\link{update_y_GLM}}
}

\examples{
# Simulate continuous data
set.seed(123)
n <- 100
p <- 150
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", 1:p)
y <- X[,1] + 0.5 * X[,2] + rnorm(n)
update_y(y = y, X = X, family = "normal", vars = c("V1","V4"))
}

