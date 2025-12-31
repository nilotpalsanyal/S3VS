\name{S3VS_GLM}

\alias{S3VS_GLM}

\title{
Structured Screen-and-Select Variable Selection in Generalized Linear Models
}

\description{
\code{S3VS_GLM} performs variable selection based on the structured screen-and-select framework in generalized linear models.
}

\usage{
S3VS_GLM(y, X, 
  method_xy = c("topk", "fixedetasqthresh", "percetasqthresh"), param_xy, 
  method_xx = c("topk", "fixedcorthresh", "perccorthresh"), param_xx, 
  vsel_method = c("NLP", "LASSO", "ENET", "SCAD", "MCP"), 
  alpha = 0.5,
  method_sel = c("conservative", "liberal"), 
  method_rem = c("conservative_begin", "conservative_end", "liberal"), 
  sel_regout = FALSE, rem_regout = FALSE, update_y_thresh = NULL, 
  m = 100, nskip = 3, verbose = FALSE, seed = NULL, parallel = FALSE)
}

\arguments{
  \item{y}{Response. A numeric/integer/logical vector with values in \{0,1\}.}
  %
  \item{X}{Design matrix of predictors. Can be a base matrix or something \code{as.matrix()} can coerce. No missing values are allowed.}
  %
  \item{method_xy}{Rule for screening some predictors as "leading variables" based on their association with the response; one of \code{c("topk", "fixedthresh", "percthresh")}. The association measure is eta-squared. 
  
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
  \item{sel_regout}{Logical. If \code{TRUE}, when predictors are selected in an iteration, the working response \code{y} is updated using the selected predictors via \code{update_y_GLM}.}
  %
  \item{rem_regout}{Logical. If \code{TRUE}, when no predictors are selected in an iteration and some are removed instead, the working response \code{y} is updated using the removed predictors via \code{update_y_GLM}.}
  %
  \item{update_y_thresh}{Numeric scalar threshold controlling how the working response \code{y} is updated when \code{sel_regout=TRUE} or \code{rem_regout=TRUE}. When \eqn{|y - fitted\_y| > update\_y\_thresh}, \code{y} is kept, else \code{y} replaced by the rounded value of fitted_y, where fitted_y is the fitted probability from the logistic model. The default value is 0.5.}
  %
  \item{m}{Integer. Maximum number of S3VS iterations to perform. Defaults to \code{100}.}
  %
  \item{nskip}{Integer. Maximum number of iterations in which no new predictors are selected before the algorithm stops. Defaults to \code{3}.}
  %
  \item{verbose}{Logical. If \code{TRUE}, prints detailed progress information at each iteration (e.g., iteration number, predictors selected or removed). Defaults to \code{FALSE}.}
  %
  \item{seed}{If supplied, sets the random seed via \code{set.seed()} to ensure reproducibility of stochastic components. If \code{NULL}, no seed is set.}
  %
  \item{parallel}{Logical. If \code{TRUE}, attempts to perform some computations in parallel mode, which is strongly recommended for faster execution. Defaults to \code{FALSE}.}
}

\details{
  For a binary response, S3VS considers the generalized linear model (GLM)
  \deqn{
    g\!\left( E\!\left( \boldsymbol{y} \mid \boldsymbol{X} \right) \right)
    = \boldsymbol{X}\boldsymbol{\beta}
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
  \code{\link{get_leadvars_GLM}}, \code{\link{get_leadsets}}, \code{\link{VS_method_GLM}}, \code{\link{select_vars}}, \code{\link{remove_vars}}, \code{\link{update_y_GLM}}
}

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
# Run S3VS for for GLM (logistic)
res_glm <- S3VS_GLM(y = y, X = X,
                method_xy = "topk", param_xy = list(k = 1),
                method_xx = "topk", param_xx = list(k = 3),
                vsel_method = "LASSO", 
                method_sel = "conservative", method_rem = "conservative_begin", 
                sel_regout = FALSE, rem_regout = FALSE,
                m = 100, nskip = 3, verbose = TRUE, seed = 123)
# View selected predictors
res_glm$selected
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
