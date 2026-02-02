# ================================================
# EVT on log-mass (allometry-first; subset-aware)
# ================================================
suppressPackageStartupMessages({
  library(MASSTIMATE)   # extants
  library(dplyr); library(tidyr); library(purrr)
  library(readr); library(ggplot2); library(scales)
})

# ---------- 1) Train allometries on extants ----------
data("extants")
extants_df <- extants %>%
  dplyr::transmute(
    lnM  = log(as.numeric(BM)),
    lnCF = log(as.numeric(FC)),
    lnCH = log(as.numeric(HC)),
    lnLF = log(as.numeric(Femur.Length)),
    lnLH = log(as.numeric(Humerus.Length))
  ) %>%
  dplyr::mutate(
    lnCFH = dplyr::if_else(is.finite(lnCF) & is.finite(lnCH),
                           log(exp(lnCF) + exp(lnCH)), NA_real_)
  )

fit_subset_ols <- function(train_df, predictors){
  df <- train_df %>% dplyr::select(lnM, dplyr::all_of(predictors)) %>% tidyr::drop_na()
  stopifnot(nrow(df) >= length(predictors) + 2)
  fm  <- stats::as.formula(paste("lnM ~", paste(predictors, collapse = " + ")))
  fit <- stats::lm(fm, data = df)
  X   <- stats::model.matrix(fit)
  list(
    fit   = fit,
    sigma = sqrt(sum(stats::residuals(fit)^2) / stats::df.residual(fit)),
    XTXi  = tryCatch(solve(t(X) %*% X), error = function(e) NULL),
    predictors = predictors
  )
}

choose_predictors <- function(avail){
  priority <- list(
    c("lnCFH","lnLF","lnLH"),
    c("lnCFH","lnLF"), c("lnCFH","lnLH"), c("lnCFH"),
    c("lnCF","lnLF","lnLH"), c("lnCH","lnLF","lnLH"),
    c("lnCF","lnLF"), c("lnCH","lnLH"),
    c("lnCF"), c("lnCH"),
    c("lnLF","lnLH"), c("lnLF"), c("lnLH")
  )
  for (p in priority) if (all(p %in% avail)) return(p)
  character(0)
}

ols_cache <- new.env(parent = emptyenv())
get_or_fit_model <- function(preds){
  key <- paste(sort(preds), collapse = "+")
  if (exists(key, envir = ols_cache)) return(get(key, envir = ols_cache))
  mod <- fit_subset_ols(extants_df, preds); assign(key, mod, envir = ols_cache); mod
}

predict_logmass_normal <- function(row_named){
  need <- c("lnCFH","lnCF","lnCH","lnLF","lnLH")
  x <- as.list(row_named); for (nm in need) if (is.null(x[[nm]])) x[[nm]] <- NA_real_
  if (!is.finite(x$lnCFH) && is.finite(x$lnCF) && is.finite(x$lnCH))
    x$lnCFH <- log(exp(x$lnCF) + exp(x$lnCH))
  avail <- need[is.finite(unlist(x[need]))]
  preds <- choose_predictors(avail)
  if (!length(preds)) return(list(mu = NA_real_, sd = NA_real_, used = character(0)))
  mod <- get_or_fit_model(preds)
  mm  <- stats::model.matrix(stats::as.formula(paste("~", paste(preds, collapse = "+"))),
                             as.data.frame(x))
  mu  <- as.numeric(cbind(1, mm) %*% stats::coef(mod$fit))
  if (is.null(mod$XTXi)) {
    Vb <- stats::vcov(mod$fit)
    v  <- as.numeric(cbind(1, mm) %*% Vb %*% t(cbind(1, mm)))
    sd <- sqrt(mod$sigma^2 + v)
  } else {
    v  <- as.numeric(cbind(1, mm) %*% mod$XTXi %*% t(cbind(1, mm)))
    sd <- sqrt(mod$sigma^2 * (1 + v))
  }
  list(mu = mu, sd = sd, used = preds)
}

# ---------- 2) Load fossils (explicit mapping) ----------
parse_mm <- function(x) readr::parse_number(as.character(x))

fossils_raw <- readr::read_csv(
  "Data/DEmic23_updated_Supplemental_Data_withPubYear_withAvailability.csv",
  show_col_types = FALSE
)

fossils <- fossils_raw %>%
  dplyr::transmute(
    Species = `genus and species`,
    year    = `publication year`,
    CF      = parse_mm(`femur circ (mm)`),
    CH      = parse_mm(`humerus circ (mm)`),
    LF      = parse_mm(`femur L`),
    LH      = parse_mm(`humerus L`),
    CFH_sum = parse_mm(`hum+fem circ (mm)`)
  ) %>%
  dplyr::mutate(
    lnCF  = log(CF),
    lnCH  = log(CH),
    lnLF  = log(LF),
    lnLH  = log(LH),
    lnCFH = dplyr::case_when(
      is.finite(CFH_sum) ~ log(CFH_sum),
      is.finite(lnCF) & is.finite(lnCH) ~ log(exp(lnCF) + exp(lnCH)),
      TRUE ~ NA_real_
    )
  )

# ---------- 3) Per-fossil log-mass posteriors ----------
cols_pred <- c("lnCFH","lnCF","lnCH","lnLF","lnLH")
logmass_post <- fossils %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    tmp  = list(predict_logmass_normal(as.list(dplyr::pick(dplyr::all_of(cols_pred))))),
    mu   = tmp$mu,
    sd   = tmp$sd,
    used_predictors = paste(tmp$used, collapse = "+")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-tmp) %>%
  dplyr::filter(is.finite(mu), is.finite(sd))

# ---------- 4) Multiple imputation (soft → crisp) ----------
set.seed(42)
S_impute <- 2000L
draws_tbl <- logmass_post %>%
  dplyr::mutate(draws = purrr::map2(mu, sd, ~rnorm(S_impute, .x, .y)))
mi_sample <- draws_tbl %>%
  dplyr::transmute(Species, year, logM = draws) %>%
  tidyr::unnest(logM)

# ---------- 5) Threshold selection (log-mass) ----------
cand_u <- quantile(mi_sample$logM, probs = seq(0.70, 0.95, by = 0.01), na.rm = TRUE) %>% as.numeric()

negloglik_gpd <- function(par, y, u){
  xi <- par[1]; sig <- par[2]
  if (sig <= 0 || xi >= 0) return(Inf)
  y_exc <- y[y > u]; if (length(y_exc) < 10) return(Inf)
  if (u - sig/xi <= max(y_exc)) return(Inf)
  z <- y_exc - u; t <- 1 + xi * z / sig; if (any(t <= 0)) return(Inf)
  n <- length(z)
  -( -n * log(sig) + (-1/xi - 1) * sum(log(t)) )
}
fit_gpd_weibull <- function(y, u){
  y_exc <- y[y > u]; if (length(y_exc) < 10) return(NULL)
  init <- c(xi = -0.2, sigma = stats::sd(y_exc - u))
  fit  <- tryCatch(stats::optim(init, negloglik_gpd, y = y, u = u, method = "L-BFGS-B",
                                lower = c(-5, 1e-8), upper = c(-1e-3, Inf)),
                   error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$value)) return(NULL)
  xi <- fit$par[1]; sig <- fit$par[2]
  list(u = u, xi = xi, sigma_u = sig, adj_scale = sig - xi * u,
       ystar = u - sig/xi, n = sum(y > u))
}
stab <- purrr::map(cand_u, ~fit_gpd_weibull(mi_sample$logM, .x)) %>%
  purrr::discard(is.null) %>% dplyr::bind_rows()
stopifnot(nrow(stab) >= 2)
stab <- stab %>%
  dplyr::mutate(dx = c(NA, diff(xi)), ds = c(NA, diff(adj_scale)),
                score = abs(dx) + abs(ds)) %>%
  dplyr::filter(n >= 25)
u_star  <- stab$u[which.min(stab$score)]
fit_star <- fit_gpd_weibull(mi_sample$logM, u_star)

# ---------- 6) Exceedance probabilities per specimen ----------
exceed_tbl <- logmass_post %>%
  dplyr::mutate(p_exc = 1 - pnorm((u_star - mu) / sd)) %>%
  dplyr::arrange(dplyr::desc(p_exc))

# ---------- 7) Figures ----------
theme_science <- ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = element_blank(),
                 panel.grid.major = element_line(color = "gray85"),
                 panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
                 axis.ticks = element_line(color = "black", linewidth = 0.3))

# (a) Top 10 densities + threshold
top10 <- exceed_tbl %>% dplyr::slice_max(p_exc, n = 10) %>% dplyr::pull(Species)
dens_df <- draws_tbl %>% dplyr::filter(Species %in% top10) %>%
  tidyr::unnest(draws) %>% dplyr::rename(logM = draws)

p_top10 <- ggplot2::ggplot(dens_df, aes(x = logM, y = after_stat(scaled))) +
  ggplot2::geom_density(fill = "steelblue", alpha = 0.35, linewidth = 0.6, color = "steelblue") +
  ggplot2::geom_vline(xintercept = u_star, linetype = "dashed", color = "firebrick", linewidth = 1) +
  ggplot2::facet_wrap(~ Species, scales = "free_y") +
  ggplot2::labs(x = "log(Mass)", y = "Scaled density",
                title = "Top 10 sauropods: log-mass distributions (dashed = threshold)") +
  theme_science

# (b) All specimens: exceedance bars
p_excbars <- exceed_tbl %>%
  dplyr::mutate(Species = reorder(Species, p_exc)) %>%
  ggplot2::ggplot(aes(Species, p_exc)) +
  ggplot2::geom_col(fill = "gray40") + ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = NULL, y = "Pr[log(M) > threshold]",
                title = "Exceedance probabilities for all sauropod specimens") +
  theme_science

dir.create("Figures_v4", showWarnings = FALSE)
ggplot2::ggsave("Figures_v4/top10_logmass_densities.png", p_top10, width = 10, height = 8, dpi = 300)
ggplot2::ggsave("Figures_v4/exceedance_bars.png", p_excbars, width = 8, height = 12, dpi = 300)

# ---------- 8) Endpoint report ----------
endpoint_tons <- exp(fit_star$ystar) / 1e6
cat("\n--- EVT Endpoint (MLE, Weibull GPD) ---\n",
    "Threshold u*:        ", round(u_star, 4), " (log-mass)\n",
    "xi (shape):          ", round(fit_star$xi, 4), "\n",
    "sigma_u (scale):     ", round(fit_star$sigma_u, 4), "\n",
    "Endpoint y* (log M): ", round(fit_star$ystar, 4), "\n",
    "Endpoint M* (tons):  ", round(endpoint_tons, 1), " t\n", sep = "")
