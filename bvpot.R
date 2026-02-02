## ===============================================================
## Bivariate censored logistic POT in (y*, xi) parameterisation
## ---------------------------------------------------------------

library(evd)     # fpot, dgpd
library(copula)  # Gumbel / logistic EV copula for simulation

prep_censored_pot <- function(X, u) {
  X <- as.matrix(X)
  if (ncol(X) != 2L)
    stop("X must be an n x 2 matrix.")
  
  y1 <- X[, 1]
  y2 <- X[, 2]
  u1 <- u[1]
  u2 <- u[2]
  
  n <- length(y1)
  if (length(y2) != n) stop("X must have same length in both margins.")
  
  ## Complete cases only for the bivariate censored part
  cc   <- is.finite(y1) & is.finite(y2)
  y1cc <- y1[cc]
  y2cc <- y2[cc]
  n_cc <- sum(cc)
  
  ## Exceedance indicators vs thresholds
  exc1 <- y1cc > u1
  exc2 <- y2cc > u2
  
  ## Union of exceedances (at least one margin > threshold)
  in_union <- exc1 | exc2
  nn       <- sum(in_union)
  
  if (nn == 0L) {
    warning("No union exceedances in complete pairs; nothing for censored POT.")
    return(list(
      x1     = numeric(0L),
      x2     = numeric(0L),
      thdi   = numeric(0L),
      lambda = c(NA_real_, NA_real_),
      nn     = 0L,
      n      = n_cc,
      index  = which(cc)[0L]
    ))
  }
  
  ## Restrict to union exceedances
  idx_union <- which(cc)[in_union]   # indices in original X
  y1u       <- y1[idx_union]
  y2u       <- y2[idx_union]
  exc1u     <- y1u > u1
  exc2u     <- y2u > u2
  
  ## Threshold–distance (signed)
  x1 <- ifelse(exc1u, y1u - u1, u1 - y1u)
  x2 <- ifelse(exc2u, y2u - u2, u2 - y2u)
  
  ## thdi: quadrant indicator
  ##   < 1.5: margin 1 only
  ##   1.5–2.5: margin 2 only
  ##   ≥ 2.5: both
  thdi <- numeric(nn)
  thdi[ exc1u & !exc2u] <- 1
  thdi[!exc1u &  exc2u] <- 2
  thdi[ exc1u &  exc2u] <- 3
  
  ## lambda: empirical tail probs per margin over ALL complete pairs
  lambda1 <- sum(exc1) / n_cc
  lambda2 <- sum(exc2) / n_cc
  
  list(
    x1     = x1,
    x2     = x2,
    thdi   = thdi,
    lambda = c(lambda1, lambda2),
    nn     = nn,
    n      = n_cc,
    index  = idx_union
  )
}

## ---------------------------------------------------------------
## 1. Univariate helpers in (xi, y*) parametrisation
## ---------------------------------------------------------------

## Survival above u in (xi, y*) parameterisation:
##   S(y | Y>u) = ((y* - y)/(y* - u))^(-1/xi),  u < y < y*
gpd_survival_endpoint <- function(y, u, xi, y_star) {
  if (xi >= 0) stop("gpd_survival_endpoint: xi must be < 0 (Weibull domain).")
  if (y_star <= u) stop("gpd_survival_endpoint: y_star must be > u.")
  
  term <- (y_star - y) / (y_star - u)
  if (any(term <= 0, na.rm = TRUE)) {
    stop("gpd_survival_endpoint: y must satisfy u < y < y_star.")
  }
  term^(-1 / xi)
}

## Log-density of exceedances Y>u, still using implied σ = (y* - u)(-xi)
log_gpd_exceed_endpoint <- function(y, u, xi, y_star) {
  if (xi >= 0) stop("log_gpd_exceed_endpoint: xi must be < 0 (Weibull domain).")
  if (y_star <= u) stop("log_gpd_exceed_endpoint: y_star must be > u.")
  
  sigma <- (y_star - u) * (-xi)
  if (sigma <= 0) stop("log_gpd_exceed_endpoint: implied sigma <= 0.")
  
  excess <- y - u
  if (any(excess <= 0, na.rm = TRUE)) {
    stop("log_gpd_exceed_endpoint: all y must be > u for exceedance density.")
  }
  
  evd::dgpd(excess, loc = 0, scale = sigma, shape = xi, log = TRUE)
}

## GPD tail quantile in (xi,y*)
## For V ~ U(0,1), Y | Y>u has CDF:
##   F(y) = 1 - ((y* - y)/(y* - u))^(-1/xi)
## Inverse:
##   y(V) = y* - (y* - u) * (1 - V)^(-xi)
gpd_quantile_endpoint <- function(v, u, xi, y_star) {
  if (xi >= 0) stop("gpd_quantile_endpoint: xi must be < 0 (Weibull domain).")
  if (y_star <= u) stop("gpd_quantile_endpoint: y_star must be > u.")
  if (any(v <= 0 | v >= 1)) stop("gpd_quantile_endpoint: v must be in (0,1).")
  
  y_star - (y_star - u) * (1 - v)^(-xi)
}

## ---------------------------------------------------------------
## 2. Censored logistic neg log-likelihood in (y*, xi)
## ---------------------------------------------------------------
## Arguments:
##   cpot    : list from prep_censored_pot()
##   u       : thresholds c(u1,u2)
##   y1_star, xi1, y2_star, xi2 : endpoint + shape
##   dep     : logistic dependence in (0,1]
## ---------------------------------------------------------------
compute_negloglik_censored_logistic_yxi <- function(
    cpot,
    u,
    y1_star, xi1,
    y2_star, xi2,
    dep,
    penalty = 1e6
) {
  if (dep < 0.1 || dep > 1) return(penalty)
  
  u1 <- u[1]
  u2 <- u[2]
  
  ## Implied σ's (never treated as parameters)
  sigma1 <- (y1_star - u1) * (-xi1)
  sigma2 <- (y2_star - u2) * (-xi2)
  if (sigma1 <= 0 || sigma2 <= 0) return(penalty)
  
  z1_cpot   <- as.double(cpot$x1)
  z2_cpot   <- as.double(cpot$x2)
  censor_id <- as.double(cpot$thdi)
  lambda    <- as.double(cpot$lambda)
  
  n_union <- as.integer(cpot$nn)
  n_total <- as.integer(cpot$n)
  
  ## If no union exceedances: only "both below" term
  if (n_union == 0L) {
    lambda2  <- -1 / log(1 - lambda)
    lambda2  <- lambda2^(-1 / dep)
    sum_l2   <- lambda2[1] + lambda2[2]
    zdn      <- sum_l2^(dep - 1)
    zdn      <- -zdn * sum_l2
    return(-(n_total - n_union) * zdn)
  }
  
  ## 1. zdn term: "both below threshold" region
  lambda2  <- -1 / log(1 - lambda)
  lambda2  <- lambda2^(-1 / dep)
  sum_l2   <- lambda2[1] + lambda2[2]
  zdn      <- sum_l2^(dep - 1)
  zdn      <- -zdn * sum_l2
  
  nn <- n_union
  
  ## 2. GP transform with signed threshold distances / implied σ
  y1 <- z1_cpot / sigma1
  y2 <- z2_cpot / sigma2
  
  ## GP transform t_k (Coles "cpot" style)
  if (xi1 == 0) {
    t1_gp <- exp(-y1)
  } else {
    tmp1 <- 1 + xi1 * y1
    if (any(tmp1 <= 0)) return(penalty)
    t1_gp <- tmp1^(-1 / xi1)
  }
  
  if (xi2 == 0) {
    t2_gp <- exp(-y2)
  } else {
    tmp2 <- 1 + xi2 * y2
    if (any(tmp2 <= 0)) return(penalty)
    t2_gp <- tmp2^(-1 / xi2)
  }
  
  ## 3. Logistic transform of margins
  z1_log <- -1 / log(1 - lambda[1] * t1_gp)
  z2_log <- -1 / log(1 - lambda[2] * t2_gp)
  
  ## 4. Jacobians (y -> GP -> logistic) using implied σ's
  jac1 <- (z1_log^2) * (t1_gp^(1 + xi1)) / (1 - lambda[1] * t1_gp)
  jac1 <- lambda[1] * jac1 / sigma1
  
  jac2 <- (z2_log^2) * (t2_gp^(1 + xi2)) / (1 - lambda[2] * t2_gp)
  jac2 <- lambda[2] * jac2 / sigma2
  
  ## 5. Logistic exponent measure V and derivatives
  v1_raw  <- z1_log^(-1 / dep)
  v2_raw  <- z2_log^(-1 / dep)
  v12_pow <- (v1_raw + v2_raw)^(dep - 1)
  V_val   <- v12_pow * (v1_raw + v2_raw)   # (v1+v2)^dep
  
  dV_dz1  <- -(v1_raw / z1_log) * v12_pow
  dV_dz2  <- -(v2_raw / z2_log) * v12_pow
  
  d2V_dz1dz2 <- (1 - 1 / dep) * dV_dz1 * dV_dz2 / V_val
  
  ## 6. Observation-wise contributions
  dvec <- numeric(nn)
  
  idx_m1_only <- censor_id < 1.5
  idx_m2_only <- censor_id >= 1.5 & censor_id < 2.5
  idx_both    <- censor_id >= 2.5
  
  if (any(idx_m1_only)) {
    dvec[idx_m1_only] <-
      log(-dV_dz1[idx_m1_only]) +
      log(jac1[idx_m1_only]) -
      V_val[idx_m1_only]
  }
  
  if (any(idx_m2_only)) {
    dvec[idx_m2_only] <-
      log(-dV_dz2[idx_m2_only]) +
      log(jac2[idx_m2_only]) -
      V_val[idx_m2_only]
  }
  
  if (any(idx_both)) {
    joint_term <- dV_dz1[idx_both] * dV_dz2[idx_both] -
      d2V_dz1dz2[idx_both]
    dvec[idx_both] <-
      log(joint_term) +
      log(jac1[idx_both]) +
      log(jac2[idx_both]) -
      V_val[idx_both]
  }
  
  if (any(!is.finite(dvec))) return(penalty)
  
  neg_loglik <- -sum(dvec) - (n_total - n_union) * zdn
  if (!is.finite(neg_loglik)) neg_loglik <- penalty
  
  neg_loglik
}

## ---------------------------------------------------------------
## 3. Total nll() in (xi, y*) parametrisation
## ---------------------------------------------------------------
## p =
##   if !common_shape:
##     (xi1, y1_star, xi2, y2_star, dep)
##   if  common_shape:
##     (xi,  y1_star,      y2_star, dep)
## ---------------------------------------------------------------
nll <- function(p, X, u, common_shape = FALSE, penalty = 1e6) {
  X <- as.matrix(X)
  if (ncol(X) != 2)
    stop("X must have exactly 2 columns (two margins).")
  
  y1 <- X[, 1]
  y2 <- X[, 2]
  u1 <- u[1]
  u2 <- u[2]
  
  if (!common_shape) {
    xi1     <- p[1]
    y1_star <- p[2]
    xi2     <- p[3]
    y2_star <- p[4]
    dep     <- p[5]
  } else {
    xi      <- p[1]
    y1_star <- p[2]
    y2_star <- p[3]
    dep     <- p[4]
    xi1     <- xi
    xi2     <- xi
  }
  
  ## Basic parameter constraints (Weibull domain)
  if (xi1 >= 0 || xi2 >= 0) return(penalty)
  if (y1_star <= u1 || y2_star <= u2) return(penalty)
  
  ## Endpoints must lie above all observed values
  if (any(is.finite(y1) & y1 >= y1_star, na.rm = TRUE)) return(penalty)
  if (any(is.finite(y2) & y2 >= y2_star, na.rm = TRUE)) return(penalty)
  
  ## Observation patterns
  obs1 <- is.finite(y1)
  obs2 <- is.finite(y2)
  
  idx_complete       <- obs1 & obs2
  idx_y1_only_exceed <- obs1 & !obs2 & (y1 > u1)
  idx_y2_only_exceed <- !obs1 & obs2 & (y2 > u2)
  
  ## Bivariate censored component (complete pairs only)
  if (any(idx_complete)) {
    X_complete <- cbind(y1[idx_complete], y2[idx_complete])
    cpot <- prep_censored_pot(X_complete, u)
    
    neglog_biv <- compute_negloglik_censored_logistic_yxi(
      cpot    = cpot,
      u       = u,
      y1_star = y1_star,
      xi1     = xi1,
      y2_star = y2_star,
      xi2     = xi2,
      dep     = dep,
      penalty = penalty
    )
    
    loglik_biv <- -neglog_biv
  } else {
    loglik_biv <- 0
  }
  
  ## Univariate contributions (missing partners; Y_k > u_k)
  if (any(idx_y1_only_exceed)) {
    ll1 <- sum(
      log_gpd_exceed_endpoint(
        y      = y1[idx_y1_only_exceed],
        u      = u1,
        xi     = xi1,
        y_star = y1_star
      )
    )
  } else {
    ll1 <- 0
  }
  
  if (any(idx_y2_only_exceed)) {
    ll2 <- sum(
      log_gpd_exceed_endpoint(
        y      = y2[idx_y2_only_exceed],
        u      = u2,
        xi     = xi2,
        y_star = y2_star
      )
    )
  } else {
    ll2 <- 0
  }
  
  loglik_total <- loglik_biv + ll1 + ll2
  nll_value    <- -loglik_total
  
  if (!is.finite(nll_value)) nll_value <- penalty
  nll_value
}

## ---------------------------------------------------------------
## 4. fit_bvpot(): optimisation in (xi, y*) space
## ---------------------------------------------------------------
fit_bvpot <- function(
    X, u,
    common_shape = FALSE,
    start        = NULL,
    method       = "Nelder-Mead",
    control      = list(),
    penalty      = 1e6
) {
  X <- as.matrix(X)
  if (ncol(X) != 2L)
    stop("X must be an n x 2 matrix or data.frame with two margins.")
  
  y1 <- X[, 1]
  y2 <- X[, 2]
  u1 <- u[1]
  u2 <- u[2]
  
  ## Build (xi,y*) starts from univariate POT (σ,ξ) → y* = u - σ/ξ
  init_margin <- function(y_all, u_k, fit_k) {
    sigma_ml <- unname(fit_k$estimate["scale"])
    xi_ml    <- unname(fit_k$estimate["shape"])
    
    ## Force negative ξ if ML says otherwise
    xi_start <- if (is.finite(xi_ml) && xi_ml < -1e-3) xi_ml else -0.1
    y_star   <- u_k - sigma_ml / xi_start
    
    ymax <- max(y_all, na.rm = TRUE)
    if (!is.finite(y_star) || y_star <= ymax) {
      y_star <- ymax + 0.1 * max(1, ymax - u_k)
    }
    list(xi = xi_start, y_star = y_star)
  }
  
  ## 1. Starting values in (xi, y*) space
  if (is.null(start)) {
    y1_all <- y1[is.finite(y1)]
    y2_all <- y2[is.finite(y2)]
    
    fit1 <- evd::fpot(y1_all, threshold = u1, model = "gpd")
    fit2 <- evd::fpot(y2_all, threshold = u2, model = "gpd")
    
    init1 <- init_margin(y1_all, u1, fit1)
    init2 <- init_margin(y2_all, u2, fit2)
    
    dep_start <- 0.8
    
    if (!common_shape) {
      start <- c(
        xi1     = init1$xi,
        y1_star = init1$y_star,
        xi2     = init2$xi,
        y2_star = init2$y_star,
        dep     = dep_start
      )
    } else {
      xi_start <- mean(c(init1$xi, init2$xi), na.rm = TRUE)
      if (!is.finite(xi_start) || xi_start >= -1e-3) xi_start <- -0.1
      
      start <- c(
        xi      = xi_start,
        y1_star = init1$y_star,
        y2_star = init2$y_star,
        dep     = dep_start
      )
    }
  } else {
    start <- as.numeric(start)
  }
  
  start_named <- start
  if (is.null(names(start_named))) {
    if (!common_shape) {
      names(start_named) <- c("xi1", "y1_star", "xi2", "y2_star", "dep")
    } else {
      names(start_named) <- c("xi", "y1_star", "y2_star", "dep")
    }
  }
  
  ## 2. Optimisation in (xi, y*) space
  opt <- optim(
    par     = unname(start_named),
    fn      = nll,
    X       = X,
    u       = u,
    common_shape = common_shape,
    penalty = penalty,
    method  = method,
    control = control
  )
  
  ## 3. Decode (xi, y*) parameters
  if (!common_shape) {
    par_ystar <- c(
      xi1     = opt$par[1],
      y1_star = opt$par[2],
      xi2     = opt$par[3],
      y2_star = opt$par[4],
      dep     = opt$par[5]
    )
  } else {
    par_ystar <- c(
      xi      = opt$par[1],
      y1_star = opt$par[2],
      y2_star = opt$par[3],
      dep     = opt$par[4]
    )
  }
  
  endpoints_hat <- bvpot_endpoints(par_ystar, common_shape = common_shape)
  
  res <- list(
    par_ystar    = par_ystar,      # (xi,y*) parametrisation
    endpoints    = endpoints_hat,  # c(y1_star, y2_star)
    loglik       = -opt$value,
    nll          = opt$value,
    common_shape = common_shape,
    convergence  = opt$convergence,
    counts       = opt$counts,
    message      = opt$message,
    start_par    = start_named
  )
  
  ## Optional Hessian / vcov in (xi,y*) space
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    H <- try(
      numDeriv::hessian(
        func = nll,
        x    = opt$par,
        X    = X,
        u    = u,
        common_shape = common_shape,
        penalty      = penalty
      ),
      silent = TRUE
    )
    
    if (!inherits(H, "try-error") && all(is.finite(H))) {
      res$hessian_ystar <- H
      vc <- try(solve(H), silent = TRUE)
      if (!inherits(vc, "try-error") && all(is.finite(vc))) {
        res$vcov_ystar <- vc
      }
    }
  }
  
  class(res) <- "bvpot_fit"
  res
}

## ---------------------------------------------------------------
## 5. Endpoints helper from (xi,y*) vector
## ---------------------------------------------------------------
bvpot_endpoints <- function(par_ystar, common_shape = FALSE) {
  if (!common_shape) {
    c(y1_star = unname(par_ystar["y1_star"]),
      y2_star = unname(par_ystar["y2_star"]))
  } else {
    c(y1_star = unname(par_ystar["y1_star"]),
      y2_star = unname(par_ystar["y2_star"]))
  }
}

## ---------------------------------------------------------------
## 6. Parametric bootstrap in (xi,y*) parametrisation
##    - simulate under fitted (xi,y*)
##    - refit via fit_bvpot()
## ---------------------------------------------------------------
bvpot_bootstrap <- function(
    fit,
    X,
    u,
    B       = 1000L,
    seed    = NULL,
    method  = "Nelder-Mead",
    penalty = 1e6,
    verbose = TRUE
) {
  if (!inherits(fit, "bvpot_fit")) {
    stop("'fit' must be an object returned by fit_bvpot().")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.matrix(X)
  if (ncol(X) != 2L)
    stop("X must be an n x 2 matrix or data.frame with two margins.")
  
  y1 <- X[, 1]
  y2 <- X[, 2]
  u1 <- u[1]
  u2 <- u[2]
  
  common_shape <- isTRUE(fit$common_shape)
  par_yhat     <- fit$par_ystar
  
  if (!common_shape) {
    xi1_hat     <- par_yhat["xi1"]
    y1_star_hat <- par_yhat["y1_star"]
    xi2_hat     <- par_yhat["xi2"]
    y2_star_hat <- par_yhat["y2_star"]
    dep_hat     <- par_yhat["dep"]
  } else {
    xi_hat      <- par_yhat["xi"]
    xi1_hat     <- xi_hat
    xi2_hat     <- xi_hat
    y1_star_hat <- par_yhat["y1_star"]
    y2_star_hat <- par_yhat["y2_star"]
    dep_hat     <- par_yhat["dep"]
  }
  
  if (verbose) {
    cat("\n[bvpot_bootstrap] Using fitted (xi,y*) parameters:\n")
    print(par_yhat)
    cat("[bvpot_bootstrap] common_shape =", common_shape, "\n")
  }
  
  if (is.na(dep_hat) || dep_hat <= 0 || dep_hat > 1) {
    stop("Bootstrap requires a valid logistic dep parameter dep in (0,1].")
  }
  
  ## Missingness patterns
  obs1 <- is.finite(y1)
  obs2 <- is.finite(y2)
  
  idx_cc     <- obs1 & obs2
  idx_1only  <- obs1 & !obs2
  idx_2only  <- !obs1 & obs2
  
  n_cc     <- sum(idx_cc)
  n_1only  <- sum(idx_1only)
  n_2only  <- sum(idx_2only)
  
  if (n_cc == 0L) {
    stop("No complete pairs in X: bivariate logistic structure cannot be bootstrapped.")
  }
  
  ## Tail fractions per margin
  n1_total   <- sum(obs1)
  n2_total   <- sum(obs2)
  tail_frac1 <- sum(obs1 & (y1 > u1)) / n1_total
  tail_frac2 <- sum(obs2 & (y2 > u2)) / n2_total
  
  if (verbose) {
    cat("[bvpot_bootstrap] Tail fractions:\n")
    cat("  margin1: P(Y1 > u1) ≈", tail_frac1, "(u1 =", u1, ")\n")
    cat("  margin2: P(Y2 > u2) ≈", tail_frac2, "(u2 =", u2, ")\n")
  }
  
  ## Sub-threshold empirical samples for each pattern
  sub1_cc    <- y1[idx_cc    & y1 <= u1]
  sub2_cc    <- y2[idx_cc    & y2 <= u2]
  sub1_1only <- y1[idx_1only & y1 <= u1]
  sub2_2only <- y2[idx_2only & y2 <= u2]
  
  ## Logistic EV → Gumbel copula
  theta_gumbel <- 1 / dep_hat
  gumbel_cop   <- gumbelCopula(param = theta_gumbel, dim = 2)
  
  ## Output containers (in (xi,y*) space)
  if (!common_shape) {
    par_boot <- matrix(NA_real_, nrow = B, ncol = 5,
                       dimnames = list(NULL,
                                       c("xi1", "y1_star", "xi2", "y2_star", "dep")))
  } else {
    par_boot <- matrix(NA_real_, nrow = B, ncol = 4,
                       dimnames = list(NULL,
                                       c("xi", "y1_star", "y2_star", "dep")))
  }
  
  end_boot <- matrix(NA_real_, nrow = B, ncol = 2,
                     dimnames = list(NULL, c("y1_star", "y2_star")))
  conv     <- logical(B)
  
  ## Starting values for refits in (xi,y*)
  start_refit <- par_yhat
  
  if (verbose) {
    cat("[bvpot_bootstrap] Starting parametric bootstrap with B =", B, "...\n")
  }
  
  for (b in seq_len(B)) {
    ## 1) Simulate complete pairs via Gumbel copula + GPD tails in (xi,y*)
    U_cc  <- rCopula(n_cc, gumbel_cop)
    Y1_cc <- numeric(n_cc)
    Y2_cc <- numeric(n_cc)
    
    for (k in seq_len(n_cc)) {
      u1k <- U_cc[k, 1]
      u2k <- U_cc[k, 2]
      
      ## Margin 1
      if (u1k <= 1 - tail_frac1 || length(sub1_cc) == 0L) {
        Y1_cc[k] <- if (length(sub1_cc)) sample(sub1_cc, 1L, replace = TRUE) else u1
      } else {
        v1   <- (u1k - (1 - tail_frac1)) / tail_frac1  # ∈ (0,1)
        Y1_cc[k] <- gpd_quantile_endpoint(v1, u = u1,
                                          xi = xi1_hat,
                                          y_star = y1_star_hat)
      }
      
      ## Margin 2
      if (u2k <= 1 - tail_frac2 || length(sub2_cc) == 0L) {
        Y2_cc[k] <- if (length(sub2_cc)) sample(sub2_cc, 1L, replace = TRUE) else u2
      } else {
        v2   <- (u2k - (1 - tail_frac2)) / tail_frac2  # ∈ (0,1)
        Y2_cc[k] <- gpd_quantile_endpoint(v2, u = u2,
                                          xi = xi2_hat,
                                          y_star = y2_star_hat)
      }
    }
    
    ## 2) Margin 1 only
    if (n_1only > 0L) {
      U1_only <- runif(n_1only)
      Y1_only <- numeric(n_1only)
      for (k in seq_len(n_1only)) {
        u1k <- U1_only[k]
        if (u1k <= 1 - tail_frac1 || length(sub1_1only) == 0L) {
          Y1_only[k] <- if (length(sub1_1only)) sample(sub1_1only, 1L, replace = TRUE) else u1
        } else {
          v1   <- (u1k - (1 - tail_frac1)) / tail_frac1
          Y1_only[k] <- gpd_quantile_endpoint(v1, u = u1,
                                              xi = xi1_hat,
                                              y_star = y1_star_hat)
        }
      }
    } else {
      Y1_only <- numeric(0L)
    }
    
    ## 3) Margin 2 only
    if (n_2only > 0L) {
      U2_only <- runif(n_2only)
      Y2_only <- numeric(n_2only)
      for (k in seq_len(n_2only)) {
        u2k <- U2_only[k]
        if (u2k <= 1 - tail_frac2 || length(sub2_2only) == 0L) {
          Y2_only[k] <- if (length(sub2_2only)) sample(sub2_2only, 1L, replace = TRUE) else u2
        } else {
          v2   <- (u2k - (1 - tail_frac2)) / tail_frac2
          Y2_only[k] <- gpd_quantile_endpoint(v2, u = u2,
                                              xi = xi2_hat,
                                              y_star = y2_star_hat)
        }
      }
    } else {
      Y2_only <- numeric(0L)
    }
    
    ## 4) Assemble bootstrap X with same pattern
    n_boot <- n_cc + n_1only + n_2only
    X_boot <- matrix(NA_real_, nrow = n_boot, ncol = 2)
    colnames(X_boot) <- c("margin1", "margin2")
    
    if (n_cc > 0L) {
      X_boot[1:n_cc, ] <- cbind(Y1_cc, Y2_cc)
    }
    if (n_1only > 0L) {
      idx1 <- (n_cc + 1):(n_cc + n_1only)
      X_boot[idx1, 1] <- Y1_only
    }
    if (n_2only > 0L) {
      idx2 <- (n_cc + n_1only + 1):n_boot
      X_boot[idx2, 2] <- Y2_only
    }
    
    ## 5) Refit via fit_bvpot() (in (xi, y*))
    fit_b <- try(
      fit_bvpot(
        X            = X_boot,
        u            = u,
        common_shape = common_shape,
        start        = start_refit,
        method       = method,
        control      = list(),
        penalty      = penalty
      ),
      silent = TRUE
    )
    
    if (inherits(fit_b, "try-error") || fit_b$convergence != 0) {
      conv[b] <- FALSE
      next
    }
    
    conv[b]         <- TRUE
    par_b           <- fit_b$par_ystar
    par_boot[b, ]   <- par_b
    end_boot[b, ]   <- fit_b$endpoints
    
    if (verbose && (b %% max(1L, B %/% 10L) == 0L)) {
      cat("[bvpot_bootstrap] Done", b, "of", B, "\n")
    }
  }
  
  if (verbose) {
    cat("[bvpot_bootstrap] Finished. Successful refits:", sum(conv), "of", B, "\n")
  }
  
  res <- list(
    par_boot       = par_boot,      # (xi,y*) draws
    endpoints_boot = end_boot,      # y* draws
    convergence    = conv,
    common_shape   = common_shape,
    tail_frac      = c(margin1 = tail_frac1, margin2 = tail_frac2),
    B              = B,
    seed           = seed,
    fit            = fit
  )
  
  class(res) <- "bvpot_bootstrap"
  res
}

## ---------------------------------------------------------------
## 7. Usage sketch
## ---------------------------------------------------------------
## X_data <- cbind(log_SVL, log_TL)  # possibly with NAs
## u_vec  <- c(u1, u2)
##
## fit_eq <- fit_bvpot(X_data, u_vec, common_shape = TRUE)
## boot   <- bvpot_bootstrap(fit_eq, X_data, u_vec, B = 1000L, seed = 2026)
##
## par_draws   <- boot$par_boot        # (xi,y*)
## ystar_draws <- boot$endpoints_boot  # endpoints y*
## ok <- boot$convergence &
##       is.finite(par_draws[, if (fit_eq$common_shape) "xi" else "xi1"]) &
##       (par_draws[, if (fit_eq$common_shape) "xi" else "xi1"] < 0)
## ci_y1 <- quantile(ystar_draws[ok, "y1_star"], probs = 0.90, na.rm = TRUE)
## ci_y2 <- quantile(ystar_draws[ok, "y2_star"], probs = 0.90, na.rm = TRUE)

## ---------------------------------------------------------------
## 8. Joint profile likelihood in (y1*, y2*)
##    - nuisance = (xi1, xi2, dep) for free shapes
##    - nuisance = (xi, dep)       for common_shape = TRUE
## ---------------------------------------------------------------

# Profile at a single (y1_star, y2_star)
profile_ll_joint <- function(
    y1_star, y2_star,
    X, u,
    common_shape   = FALSE,
    start_nuisance = NULL,
    method         = "Nelder-Mead",
    control        = list(),
    penalty        = 1e6
) {
  X <- as.matrix(X)
  if (ncol(X) != 2L)
    stop("X must be an n x 2 matrix or data.frame with two margins.")
  
  u1 <- u[1]
  u2 <- u[2]
  
  # Endpoints must be above thresholds
  if (y1_star <= u1 || y2_star <= u2) {
    return(list(
      loglik      = NA_real_,
      nll         = NA_real_,
      par_nuis    = NULL,
      convergence = 1L,
      message     = "y*_k <= u_k"
    ))
  }
  
  # Reasonable default nuisance starts (if none supplied)
  if (is.null(start_nuisance)) {
    if (!common_shape) {
      start_nuisance <- c(xi1 = -0.1, xi2 = -0.1, dep = 0.8)
    } else {
      start_nuisance <- c(xi = -0.1, dep = 0.8)
    }
  }
  
  # Wrap nll() with fixed y1_star, y2_star
  if (!common_shape) {
    nll_wrap <- function(par_nuis) {
      xi1 <- par_nuis[1]
      xi2 <- par_nuis[2]
      dep <- par_nuis[3]
      p   <- c(xi1 = xi1,
               y1_star = y1_star,
               xi2 = xi2,
               y2_star = y2_star,
               dep = dep)
      nll(p, X = X, u = u, common_shape = FALSE, penalty = penalty)
    }
  } else {
    nll_wrap <- function(par_nuis) {
      xi  <- par_nuis[1]
      dep <- par_nuis[2]
      p   <- c(xi = xi,
               y1_star = y1_star,
               y2_star = y2_star,
               dep = dep)
      nll(p, X = X, u = u, common_shape = TRUE, penalty = penalty)
    }
  }
  
  opt <- try(
    optim(
      par     = unname(start_nuisance),
      fn      = nll_wrap,
      method  = method,
      control = control
    ),
    silent = TRUE
  )
  
  if (inherits(opt, "try-error") ||
      !is.finite(opt$value) ||
      opt$value >= penalty * 0.99) {
    return(list(
      loglik      = NA_real_,
      nll         = NA_real_,
      par_nuis    = NULL,
      convergence = 1L,
      message     = "optim failure or penalty region"
    ))
  }
  
  list(
    loglik      = -opt$value,
    nll         = opt$value,
    par_nuis    = opt$par,
    convergence = opt$convergence,
    message     = opt$message
  )
}

# Evaluate joint profile on a grid of (y1*, y2*)
profile_grid_joint <- function(
    fit,              # object from fit_bvpot()
    X, u,
    grid_y1, grid_y2, # vectors of endpoint candidates
    method       = "Nelder-Mead",
    control      = list(),
    penalty      = 1e6,
    verbose      = TRUE
) {
  if (!inherits(fit, "bvpot_fit"))
    stop("'fit' must be an object returned by fit_bvpot().")
  
  X <- as.matrix(X)
  if (ncol(X) != 2L)
    stop("X must be an n x 2 matrix or data.frame with two margins.")
  
  common_shape <- isTRUE(fit$common_shape)
  par_hat      <- fit$par_ystar
  
  # Initial nuisance start at MLE
  if (!common_shape) {
    start_nuis <- c(
      xi1 = par_hat["xi1"],
      xi2 = par_hat["xi2"],
      dep = par_hat["dep"]
    )
  } else {
    start_nuis <- c(
      xi  = par_hat["xi"],
      dep = par_hat["dep"]
    )
  }
  
  ny1 <- length(grid_y1)
  ny2 <- length(grid_y2)
  
  res <- data.frame(
    y1_star   = rep(grid_y1, each = ny2),
    y2_star   = rep(grid_y2, times = ny1),
    loglik    = NA_real_,
    nll       = NA_real_,
    conv      = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  idx <- 1L
  n_total <- ny1 * ny2
  
  for (i in seq_along(grid_y1)) {
    for (j in seq_along(grid_y2)) {
      y1s <- grid_y1[i]
      y2s <- grid_y2[j]
      
      prof <- profile_ll_joint(
        y1_star      = y1s,
        y2_star      = y2s,
        X            = X,
        u            = u,
        common_shape = common_shape,
        start_nuisance = start_nuis,
        method       = method,
        control      = control,
        penalty      = penalty
      )
      
      res$loglik[idx] <- prof$loglik
      res$nll[idx]    <- prof$nll
      res$conv[idx]   <- prof$convergence
      
      # Warm start: update nuisance start if optimisation succeeded
      if (!is.null(prof$par_nuis) && prof$convergence == 0L &&
          all(is.finite(prof$par_nuis))) {
        start_nuis <- prof$par_nuis
      }
      
      if (verbose && (idx %% max(1L, n_total %/% 10L) == 0L)) {
        cat("[profile_grid_joint] Done", idx, "of", n_total, "\n")
      }
      
      idx <- idx + 1L
    }
  }
  
  res
}

