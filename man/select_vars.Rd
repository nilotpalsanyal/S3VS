\name{select_vars}

\alias{select_vars}

\title{
Aggregate Selected Predictors Across Multiple Leading Sets
}

\description{
\code{select_vars} combines variable selections obtained from multiple leading sets into a single set, using either a liberal (union) or conservative (progressive intersection) rule.
}

\usage{
select_vars(listselect, method = c("conservative", "liberal"))
}

\arguments{
  \item{listselect}{A \code{list} of vectors, each containing names of the predictors selected in the corresponding leading set.}
  %
  \item{method}{Aggregation rule. One of \code{"conservative"} or \code{"liberal"} (partial matching allowed). Referring to the sets (vectors included in \code{listselect}) of selected predictors:
    \describe{
      \item{\code{"liberal"}}{Returns the \emph{union} of all sets: \code{unique(unlist(listselect))}.}
      \item{\code{"conservative"}}{Returns the \emph{last non-empty intersection} when intersecting the first \eqn{i=1,2,\dots} sets in order. The procedure stops once the running intersection becomes empty and returns the previous (last non-empty) intersection. If the first set is empty, returns \code{character(0)}. Order of vectors in \code{listselect} matters for this method.}
    }}
}

\details{
The liberal rule favors inclusiveness, while the conservative rule favors stability.
}

\value{
Vector with names of the retained predictors (considered selected in the current iteration of S3VS); if no predictors are retained, \code{character(0)}).
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
listselect <- list(
  c("V1","V2","V23"),
  c("V4","V2","V23"),
  c("V4","V5","V23")
)
select_vars(listselect, method="conservative")
}

