\name{looprun}

\alias{looprun}

\title{
Decide Whether to Run Another S3VS Iteration
}

\description{
\code{looprun} evaluates simple stopping criteria for the S3VS procedure and returns an indicator of whether \emph{one more iteration} should be executed.
}

\usage{
looprun(varsselected, varsleft, max_nocollect, m, nskip)
}

\arguments{
  \item{varsselected}{Character vector with names of predictors selected so far. Only its length is used; \code{NULL} is treated as length 0.}
  %
  \item{varsleft}{Character vector with names of candidate predictors that remain available for selection in future iterations. Only its length is used; \code{NULL} is treated as length 0.}
  %
  \item{max_nocollect}{Integer count of iterations \emph{up to now} in which no new predictors were selected.}
  %
  \item{m}{Maximum allowed number of selected predictors (target cap for \code{length(varsselected)}).}
  %
  \item{nskip}{Maximum allowed number of "no-collection" iterations before stopping.}
}

\details{
An additional S3VS iteration is recommended \emph{iff} all three conditions hold:
\deqn{|\texttt{varsselected}| < m,}
\deqn{|\texttt{varsleft}| > 0,}
\deqn{\texttt{max\_nocollect} < \texttt{nskip}.}
}

\value{
\code{1} if another iteration should run, \code{0} otherwise.
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
looprun(varsselected = c("x1","x2","x3"),
        varsleft     = paste0("x", 4:23),
        max_nocollect = 0,
        m = 10,
        nskip = 2)
}

