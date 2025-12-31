## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  message  = FALSE,
  warning  = FALSE
)

## -----------------------------------------------------------------------------
library(S3VS)

## -----------------------------------------------------------------------------
set.seed(1)

n <- 200
p <- 300

rho <- 0.6
Sigma <- rho ^ abs(outer(1:p, 1:p, "-"))
X_lm <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
colnames(X_lm) <- paste0("X", seq_len(p))

beta <- rep(0, p)
active <- c(5, 17, 50, 120, 201)
beta[active] <- c(2.0, -1.5, 1.2, 1.8, -2.2)

y_lm <- as.numeric(X_lm %*% beta + rnorm(n, sd = 2))

## -----------------------------------------------------------------------------
fit_lm <- S3VS(
  y = y_lm,
  X = X_lm,
  family = "normal",
  method_xy = "percthresh",
  param_xy = list(thresh = 95),
  method_xx = "topk",
  param_xx = list(k = 10),
  vsel_method = "LASSO",
  method_sel = "conservative",
  method_rem = "conservative_begin",
  verbose = FALSE,
  seed = 1
)

## -----------------------------------------------------------------------------
fit_lm$selected

## -----------------------------------------------------------------------------
str(fit_lm$selected_iterwise)

## -----------------------------------------------------------------------------
Xsel <- X_lm[, fit_lm$selected, drop = FALSE]
pred_lm <- pred_S3VS(y = y_lm, X = Xsel, family = "normal", method = "LASSO")

head(pred_lm$y.pred)

## -----------------------------------------------------------------------------
set.seed(2)

n <- 250
p <- 400

rho <- 0.4
Sigma <- rho ^ abs(outer(1:p, 1:p, "-"))
X_glm <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
colnames(X_glm) <- paste0("X", seq_len(p))

beta <- rep(0, p)
active <- c(10, 25, 90, 200)
beta[active] <- c(1.2, -1.0, 0.9, -1.4)

eta <- as.numeric(X_glm %*% beta)
pr  <- 1 / (1 + exp(-eta))
y_glm   <- rbinom(n, 1, pr)

## -----------------------------------------------------------------------------
fit_glm <- S3VS(
  y = y_glm,
  X = X_glm,
  family = "binomial",
  method_xy = "_survtopk",
  param_xy = list(k = 2),
  method_xx = "percthresh",
  param_xx = list(thresh = 90),
  vsel_method = "SCAD",
  sel_regout = TRUE,
  rem_regout = TRUE,
  update_y_thresh = 0.5,
  verbose = FALSE,
  seed = 1
)

fit_glm$selected

## -----------------------------------------------------------------------------
set.seed(3)

n <- 300
p <- 200
X_surv <- matrix(rnorm(n*p), n, p)
colnames(X_surv) <- paste0("X", seq_len(p))

beta <- rep(0, p)
active <- c(3, 30, 77, 150)
beta[active] <- c(0.6, -0.5, 0.7, -0.8)

linpred <- as.numeric(X_surv %*% beta)
base_rate <- 0.05

time <- rexp(n, rate = base_rate * exp(linpred))
cens <- rexp(n, rate = 0.02)
status <- as.integer(time <= cens)
time <- pmin(time, cens)

y_surv <- list(time = time, status = status)

## -----------------------------------------------------------------------------
fit_cox <- S3VS(
  y = y_surv,
  X = X_surv,
  family = "survival",
  surv_model = "COX",
  method_xy = "topk",
  param_xy = list(k = 2),
  method_xx = "percthresh",
  param_xx = list(thresh = 90),
  vsel_method = "LASSO",
  verbose = FALSE,
  seed = 1
)

fit_cox$selected

## -----------------------------------------------------------------------------
fit_aft <- S3VS(
  y = y_surv,
  X = X_surv,
  family = "survival",
  surv_model = "AFT",
  method_xy = "topk",
  param_xy = list(k = 2),
  method_xx = "topk",
  param_xx = list(k = 2),
  vsel_method = "AFTGEE",
  verbose = FALSE,
  seed = 1,
  parallel = FALSE
)

fit_aft$selected

## -----------------------------------------------------------------------------
# For the earlier linear model example
leadvars <- get_leadvars(y = y_lm, X = X_lm, family = "normal", varsleft = colnames(X),
                        method = "topk", param = list(k = 20))
leadvars

# For the earlier generalized linear model example
leadvars <- get_leadvars_GLM(y = y_glm, X = X_glm, method = "percetasqthresh", 
                             param = list(thresh = 80))
leadvars

# For the earlier survival Cox model
leadvars <- get_leadvars_SURV(y = y_surv, X = X_surv, method = "topk", 
                              param = list(k = 5))
leadvars

## ----eval=FALSE---------------------------------------------------------------
# cor_xy <- abs(cor(y, X))
# fit <- S3VS(y = y, X = X, family = "normal", cor_xy = cor_xy, ...)

## ----eval=FALSE---------------------------------------------------------------
# leadsets <- get_leadsets(X = X, leadvars = leadvars,
#                          method_xx = "topk", param_xx = list(k = 40))

## ----eval=FALSE---------------------------------------------------------------
# sel <- VS_method_LM(y = y, X = X[, leadsets[[1]], drop = FALSE],
#                     vsel_method = "SCAD")
# sel$selected
# 
# sel <- VS_method_GLM(y = y, X = X[, leadsets[[1]], drop = FALSE],
#                     vsel_method = "MCP")
# sel$selected
# 
# sel <- VS_method_SURV(y = y, X = X[, leadsets[[1]], drop = FALSE], surv_model = "AFT",
#                     vsel_method = "AFTGEE")
# sel$selected

## ----eval=FALSE---------------------------------------------------------------
# agg <- select_vars(sel_list, method_sel = "conservative")
# agg

## ----eval=FALSE---------------------------------------------------------------
# rem <- remove_vars(varsleft, varsselected, method_rem = "conservative_begin")
# rem

## ----eval=FALSE---------------------------------------------------------------
# y_new <- update_y_LM(y = y, X = X[, varsselected, drop = FALSE])
# 
# y_new <- update_y_GLM(y = y, X = X[, varsselected, drop = FALSE], update_y_thresh = 0.4)

## ----eval=FALSE---------------------------------------------------------------
# looprun(varsselected, varsleft, max_nocollect = 2, m = 100, nskip = 3)

## ----eval=FALSE---------------------------------------------------------------
# library(future)
# plan(multisession)
# fit_aft <- S3VS(..., family="survival", surv_model="AFT", parallel=TRUE)

