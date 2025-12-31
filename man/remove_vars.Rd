\name{remove_vars}

\alias{remove_vars}

\title{
Aggregate Not-Selected Predictors for Removal Across Multiple Leading Sets
}

\description{
\code{remove_vars} combines lists of predictors that were \emph{not selected} from multiple leading sets into a single set to remove, using either a liberal (union) rule or a conservative (progressive intersection) rule.
}

\usage{
remove_vars(listnotselect, 
  method = c("conservative_begin", "conservative_end", "liberal"))
}

\arguments{
  \item{listnotselect}{A \code{list} of vectors, each containing names of the predictors not selected in the corresponding leading set.}
  %
  \item{method}{Aggregation rule; one of \code{"conservative_begin"}, \code{"conservative_end"}, or \code{"liberal"}. Referring to the sets (vectors included in \code{listselect}) of not-selected predictors:
    \describe{
      \item{\code{"liberal"}}{Returns the \emph{union} of all sets: \code{unique(unlist(listnotselect))}.}
      \item{\code{"conservative_begin"}}{Returns the \emph{last non-empty intersection} when intersecting the \strong{first} \eqn{i=1,2,\dots} sets in order (note that, the first set is assumed to be non-empty, because that will automatically be true if \code{remove_vars} is being called by \code{S3VS} or its family-specific engines). The procedure stops once the running intersection becomes empty and returns the previous (last non-empty) intersection. Order of \code{listnotselect} matters.}
      \item{\code{"conservative_end"}}{Returns the \emph{last non-empty intersection} when intersecting the \strong{last} \eqn{i=1,2,\dots} sets in order (if the last set is empty, the function finds the first non-empty set from the end, and then, starts the intersection process from that set). The procedure stops once the running intersection becomes empty and returns the previous (last non-empty) intersection. Order of \code{listnotselect} matters.}
    }}
}

\details{
The liberal rule favors inclusiveness (drop all predictors that were not selected in an iteration), whereas the conservative rule favors stability across earlier/latter leading sets (drop only predictors consistently absent in earlier/latter leading sets).
}

\value{
Vector with names of the predictors that are not selected till the current S3VS iteration and to be removed from all future iterations.
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
remove_vars(listselect, method="liberal")
}

