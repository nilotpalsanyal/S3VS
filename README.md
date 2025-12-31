
<!-- README.md is generated from README.Rmd. Please edit that file -->

# S3VS

<!-- badges: start -->

<!-- badges: end -->

The goal of S3VS is to perform variable selection using the structured
screen-and-select (S3VS) framework in linear models, generalized linear
models with binary data, and survival models such as the Cox model and
accelerated failure time (AFT) model.

## Installation

You can install the development version of S3VS like so:

``` r
devtools::install_github("nilotpalsanyal/S3VS")
```

## Description

The central entry point is `S3VS()`, which dispatches to a
family-specific routine via the argument `family`:

- `S3VS_LM()` for linear models,
- `S3VS_GLM()` for generalized linear models with binary outcomes,
- `S3VS_SURV()` for survival models.

The S3VS workflow proceeds through the following steps, each handled by
helper functions:

- **Stopping rule check:** `looprun()` determines whether the iterative
  screen-and-select process should continue.

- **Leading variable identification:** `get_leadvars()` identifies
  leading variables; family-specific versions are `get_leadvars_LM`,
  `get_leadvars_GLM`, and `get_leadvars_SURV`.

- **Leading set identification:** `get_leadsets()` identifies the
  leading set for each leading variable.

- **Selection within leading sets:** `VS_method()` performs selection
  within leading sets; family-specific methods include `VS_method_LM()`,
  `VS_method_GLM()`, `VS_method_SURV()`, and `bridge_aft()` implements
  BRIDGE specifically for AFT models.

- **Aggregation of selected variables:** `select_vars()` retains
  promising variables as selected from an iteration.

- **Aggregation of non-selected variables:** (optional) `remove_vars()`
  removes variables deemed uninformative from future iterations (if no
  variable is selected in the current iteration by `select_vars()`.

- **Response update:** (optional) `update_y()` enables iterative
  response updates; family-specific variants include`update_y_LM()`
  and`update_y_GLM()`.

Together, these functions form a structured, iterative pipeline for
efficient variable screening and selection in high-dimensional
regression and survival analysis.

- **Prediction:** `pred_S3VS()` produces predictions using variables
  selected by S3VS, calling `pred_S3VS_LM()`, `pred_S3VS_GLM()`, or
  `pred_S3VS_SURV()` as appropriate.

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->
