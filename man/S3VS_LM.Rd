\name{S3VS_LM}

\alias{S3VS_LM}

\title{
Structured Screen-and-Select Variable Selection in Linear Models
}
\description{
\code{S3VS_LM} performs variable selection based on the structured screen-and-select framework in linear models.
}

\usage{
S3VS_LM(y, X, cor_xy = NULL, 
  method_xy = c("topk", "fixedcorthresh", "perccorthresh"), param_xy, 
  method_xx = c("topk", "fixedcorthresh", "perccorthresh"), param_xx, 
  vsel_method = c("NLP", "LASSO", "ENET", "SCAD", "MCP"), 
  alpha = 0.5,
  method_sel = c("conservative", "liberal"), 
  method_rem = c("conservative_begin", "conservative_end", "liberal"), 
  rem_regout = FALSE, 
  m = 100, nskip = 3, verbose = FALSE, seed = NULL)
}

\arguments{
  \item{y}{Response. A numeric vector.}
  %
  \item{X}{Design matrix of predictors. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{cor_xy}{Optional numeric vector of precomputed marginal correlations between \code{y} and each column of \code{X}. Used to speed up or reproduce screening by \eqn{|cor(y, X_j)|}. If \code{NULL}, correlations are computed internally.}
  %
  \item{method_xy}{Rule for screening some predictors as 'leading variables' based on their association with the response; one of \code{c("topk", "fixedthresh", "percthresh")}. The association measure is correlation. 
  
    \code{"topk"} keeps the predictors with the largest \eqn{k} association values; \code{"fixedthresh"} keeps predictors whose association is greater than or equal to a specified threshold; \code{"percthresh"} keeps predictors whose association is within a given percentage of the best.}
  %
  \item{param_xy}{Tuning parameter for \code{method_xy}. If \code{"topk"}, supply a list with an integer \code{k} (keep the top \eqn{k}). If \code{"fixedthresh"}, supply a list with a numeric threshold \code{thresh} (keep predictors with association \eqn{\ge} threshold). If \code{"percthresh"}, supply a list with a numeric percentage \code{thresh} in \eqn{(0,100]} (keep predictors with association \eqn{\ge} that percent of the highest association).}
  %
  \item{method_xx}{Rule for constructing, for each \emph{leading variable}, the set of associated predictors (the "leading set") using inter-predictor association (absolute value of the correlation coefficient); one of \code{c("topk", "fixedthresh", "percthresh")} with same interpretation as \code{method_xy}.}
  %
  \item{param_xx}{Tuning parameter for \code{method_xx}; same interpretation as \code{param_xy} but applied to inter-predictor association (absolute value of the correlation coefficient).}
  %
  \item{vsel_method}{Character string specifying the variable selection method to be used within each \emph{leading set}. Available options are \code{"NLP"}, \code{"LASSO"}, \code{"ENET"}, \code{"SCAD"}, \code{"MCP"}.}
  %
  \item{alpha}{Only used when \code{vsel_method == "ENET"}. Elastic net mixing parameter, with \eqn{\alpha \in (0,1)}.}
  %
  \item{method_sel}{Policy for aggregating predictors selected across leading sets in an iteration; one of \code{c("conservative","liberal")}. \code{"conservative"} selects the smallest admissible set of predictors by intersecting the selected sets of predictors across leading sets, beginning with all and gradually reducing from the end until a non-empty intersection is found; this ensures only predictors consistently selected across leading sets are retained. \code{"liberal"} selects the largest admissible set of predictors by taking the union of all selected sets of predictors, so any predictor chosen in at least one leading set is included. If no predictor is selected from the first leading set, the iteration does not contribute to final selection and exclusion rules (\code{method_rem}) are applied instead.}
  %
  \item{method_rem}{Policy for excluding predictors when no selections are made in an iteration; one of \code{c("conservative_begin","conservative_end","liberal")}. \code{"conservative_begin"} excludes the smallest admissible set of predictors by intersecting the non-selected sets of predictors starting from the first leading set; \code{"conservative_end"} does the same but begins from the last leading set and moves backward; \code{"liberal"} excludes the largest admissible set of predictors by taking the union of all non-selected sets of predictor. Predictors excluded under this rule are removed from subsequent iterations.}
  %
  \item{rem_regout}{Logical. If \code{TRUE}, when no predictors are selected in an iteration and some are removed instead, the working response \code{y} is updated using the removed predictors via \code{update_y_LM}.}
  %
  \item{m}{Integer. Maximum number of S3VS iterations to perform. Defaults to \code{100}.}
  %
  \item{nskip}{Integer. Maximum number of iterations in which no new predictors are selected before the algorithm stops. Defaults to \code{3}.}
  %
  \item{verbose}{Logical. If \code{TRUE}, prints detailed progress information at each iteration (e.g., iteration number, predictors selected or removed). Defaults to \code{FALSE}.}
  %
  \item{seed}{If supplied, sets the random seed via \code{set.seed()} to ensure reproducibility of stochastic components. If \code{NULL}, no seed is set.}
}

\details{
  For a continuous response, S3VS considers the linear model (LM)
  \deqn{
    \boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}
  }

  For the S3VS algorithm, see the manual of the top-level function \code{S3VS}.
}

\value{
A list with the following components:
\item{selected}{A character vector of predictor names that were selected across all iterations.}
\item{selected_iterwise}{A list recording the predictors selected at each iteration, in the order they were considered.}
\item{runtime}{Runtime in seconds.}
}

%%\references{
%%Reference to come...
%%}

\seealso{
  \code{\link{get_leadvars_LM}}, \code{\link{get_leadsets}}, \code{\link{VS_method_LM}}, \code{\link{select_vars}}, \code{\link{remove_vars}}, \code{\link{update_y_LM}}
}

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
# Run S3VS for LM
res_lm <- S3VS_LM(y = y, X = X,
               method_xy = "topk", param_xy = list(k=1),
               method_xx = "topk", param_xx = list(k=3),
               vsel_method = "LASSO", method_sel = "conservative", 
               method_rem = "conservative_begin", rem_regout = FALSE, 
               m = 100, nskip = 3, verbose = TRUE, seed = 123)
# View selected predictor
res_lm$selected
}

%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory (show via RShowDoc("KEYWORDS")):
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }
%% Use only one keyword per line.
%% For non-standard keywords, use \concept instead of \keyword:
%% \concept{ ~cpt1 }
%% \concept{ ~cpt2 }
%% Use only one concept per line.
