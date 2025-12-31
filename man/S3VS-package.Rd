\name{S3VS-package}

\alias{S3VS-package}

\docType{package}

\title{
  Structured Screen-and-Select Variable Selection in Linear, Generalized Linear, and Survival Models
}

\description{
  Performs variable selection using the structured screen-and-select (S3VS) framework in linear models, generalized linear models with binary data, and survival models such as the Cox model and accelerated failure time (AFT) model.
}

\details{
The \pkg{S3VS} package implements the Structured Screen-and-Select Variable Selection (S3VS) framework for linear models, generalized linear models with binary responses, and survival models (Cox proportional hazards and accelerated failure time models).

The central entry point is \code{\link{S3VS}}, which dispatches to a family-specific routine via the argument \code{family}:
\itemize{
  \item \code{\link{S3VS_LM}} for linear models,
  \item \code{\link{S3VS_GLM}} for generalized linear models with binary outcomes,
  \item \code{\link{S3VS_SURV}} for survival models.
}

The S3VS workflow proceeds through the following steps, each handled by helper functions:
\describe{
  \item{Stopping rule check}{\code{\link{looprun}} determines whether the iterative screen-and-select process should continue.}
  
  \item{Leading variable identification}{\code{\link{get_leadvars}} identifies leading variables; family-specific versions are \code{\link{get_leadvars_LM}}, \code{\link{get_leadvars_GLM}}, and \code{\link{get_leadvars_SURV}}.}
  
  \item{Leading set identification}{\code{\link{get_leadsets}} identifies the leading set for each leading variable.}
  
  \item{Selection within leading sets}{\code{\link{VS_method}} performs selection within leading sets; family-specific methods include \code{\link{VS_method_LM}}, \code{\link{VS_method_GLM}}, \code{\link{VS_method_SURV}}, and \code{\link{bridge_aft}} implements BRIDGE specifically for AFT models.}
  
  \item{Aggregation of selected variables}{\code{\link{select_vars}} retains promising variables as selected from an iteration.}
  
  \item{Aggregation of non-selected variables (optional)}{\code{\link{remove_vars}} removes variables deemed uninformative from future iterations (if no variable is selected in the current iteration by \code{\link{select_vars}}).}
  
  \item{Response update (optional)}{\code{\link{update_y}} enables iterative response updates; family-specific variants include \code{\link{update_y_LM}} and \code{\link{update_y_GLM}}.}
}

Together, these functions form a structured, iterative pipeline for efficient variable screening and selection in high-dimensional regression and survival analysis.

\describe{
  \item{Prediction}{\code{\link{pred_S3VS}} produces predictions using variables selected by S3VS, calling \code{\link{pred_S3VS_LM}}, \code{\link{pred_S3VS_GLM}}, or \code{\link{pred_S3VS_SURV}} as appropriate.}
}
}


\author{
Nilotpal Sanyal <nsanyal@utep.edu>, Padmore N. Prempeh <pprempeh@albany.edu>

Maintainer: Nilotpal Sanyal <nsanyal@utep.edu>
}

%%\references{
%%Reference to come...
%%}

%\keyword{package}

%% Uncomment below to imitate parts of library(help = S3VS)
%\section{The \file{DESCRIPTION} File}{\packageDESCRIPTION{S3VS}}
%\section{Documentation Index}{\packageIndices{S3VS}}

%\seealso{
%%  ~~ Optional links to other man pages, e.g. ~~
%%  ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%%  ~~ Optional simple examples of the most important functions ~~
%}
