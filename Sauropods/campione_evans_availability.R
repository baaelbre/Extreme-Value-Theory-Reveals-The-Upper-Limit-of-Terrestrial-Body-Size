#################################################
# Title: Hierarchical mass allometry 
# Description: Mass-from-bones with latent lnM,
# using the STAN programming language.
# Author: Bastiaan Aelbrecht
# Date: 15/08/2025
#################################################

library(cmdstanr)
library(dplyr)
library(tidyr)
library(posterior)
library(MASSTIMATE)

trait_levels <- c("lnCF","lnCH","lnLF","lnLH")  # J=4

# Priors (Alexander 1979)
ah_CF <- log(0.43); ah_CH <- log(0.35); ah_LF <- log(5.24); ah_LH <- log(4.24)
bh_CF <- log(0.36); bh_CH <- log(0.38); bh_LF <- log(0.36); bh_LH <- log(0.36)

alpha_hat <- c(ah_CF, ah_CH, ah_LF, ah_LH)
beta_hat  <- c(bh_CF, bh_CH, bh_LF, bh_LH)
sd_alpha <- 0.05
sd_beta  <- 0.02

# Inv-Gamma hyperparams
a0 <- 2.1
b0 <- 0.1

# Prior for lnM when latent (fossils without trustworthy volumetric estimates)
mu_M  <- 10
tau_M <- 5

###############################################
# Data preparation:
# Why the long format: 
# 1) missing bones are absent rows
# 2) we can always add more traits by 
# adding more rows (and changing J accordingly)
###############################################
data("extants")

extants_long <- extants %>%
  mutate(
    lnM  = log(as.numeric(BM)),
    lnCF = log(as.numeric(FC)), 
    lnCH = log(as.numeric(HC)),
    lnLF = log(as.numeric(Femur.Length)), 
    lnLH = log(as.numeric(Humerus.Length))
  ) %>%
  select(Species, lnM, lnCF, lnCH, lnLF, lnLH) %>%
  pivot_longer(
    cols      = starts_with("ln"),
    names_to  = "trait",
    values_to = "val"
  ) %>%
  mutate(
    is_mass  = (trait == "lnM"),
    trait_id = match(trait, trait_levels)  # NA for lnM
  )


# Keep only bone rows (likelihood part)
bone_rows <- extants_long %>%
  filter(!is_mass, is.finite(val), !is.na(Species), trait %in% trait_levels)

# Set Species IDs for those Speciess that have at least one bone
spec_levels <- unique(bone_rows$Species)
bone_rows   <- mutate(bone_rows, spec_id = match(Species, spec_levels))

# Observed lnM per included Species
lnM_tab <- extants_long %>%
  filter(is_mass, is.finite(val), Species %in% spec_levels) %>%
  distinct(Species, .keep_all = TRUE) %>%
  transmute(Species, lnM = val, spec_id = match(Species, spec_levels)) %>%
  filter(!is.na(spec_id))

N_spec   <- length(spec_levels)
lnM_data <- rep(0, N_spec)
has_mass <- rep(0L, N_spec)

if (nrow(lnM_tab) > 0) {
  lnM_data[lnM_tab$spec_id] <- lnM_tab$lnM
  has_mass[lnM_tab$spec_id] <- 1L
}

# Demote bad lnM values to latent
#bad <- which(has_mass == 1L & !is.finite(lnM_data))
#if (length(bad)) {
#  has_mass[bad] <- 0L
#  lnM_data[bad] <- 0
#}

# Assemble the Stan LIST
train_data <- list(
  J = length(trait_levels),
  N_spec = N_spec,
  N_obs  = nrow(bone_rows),
  obs_Species = bone_rows$spec_id,
  obs_trait    = match(bone_rows$trait, trait_levels),
  y            = bone_rows$val,
  
  has_mass = has_mass,
  lnM_data = lnM_data,
  
  has_vol_prior = rep(0L, N_spec),
  vol_mu = rep(0.0, N_spec),
  vol_sd = rep(1.0, N_spec),
  
  alpha_hat = alpha_hat,
  beta_hat  = beta_hat,
  sd_alpha  = sd_alpha,
  sd_beta   = sd_beta,
  
  a0 = a0, b0 = b0,
  mu_M = mu_M, tau_M = tau_M,
  
  N_pred = 0L,
  pred_Species = integer(0),
  pred_trait    = integer(0),
  
  .spec_levels = spec_levels,
  .bone_rows   = bone_rows
)

# Compile STAN model
mod <- cmdstan_model("mass_allometry.stan")

clean_stan_data <- function(lst) {
  # drop helpers (names starting with ".") and NULLs
  lst <- lst[!startsWith(names(lst), ".")]
  lst <- Filter(Negate(is.null), lst)
  
  # coerce integer-ish things to integer
  int_names <- c("J","N_spec","N_obs","obs_Species","obs_trait",
                 "has_mass","has_vol_prior","N_pred","pred_Species","pred_trait")
  for (nm in intersect(int_names, names(lst))) {
    lst[[nm]] <- as.integer(lst[[nm]])
  }
  lst
}

fit_train <- mod$sample(
  data = clean_stan_data(train_data),
  seed = 2025, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000
)


# Summary diagnostics: Rhat, ESS
posterior::summarise_draws(
  fit_train$draws(), rhat, ess_bulk, ess_tail
)

color_scheme_set("blue")

# nice labels
trait_labels <- c("lnCF","lnCH","lnLF","lnLH")
pars_alpha <- paste0("alpha[", 1:4, "]")
pars_beta  <- paste0("beta[",  1:4, "]")
pars_sig2  <- paste0("sigma2[",1:4, "]")

dr_train <- fit_train$draws()         

# α_j posterior
mcmc_areas(dr_train, pars = pars_alpha, prob = 0.8) +
  ggtitle("Posterior of intercepts α_j") +
  scale_y_discrete(labels = trait_labels)
# β_j posterior
mcmc_areas(dr_train, pars = pars_beta, prob = 0.8) +
  ggtitle("Posterior of slopes β_j") +
  scale_y_discrete(labels = trait_labels)
# σ_j posterior
mcmc_areas(dr_train, pars = paste0("sigma[", 1:4, "]"), prob = 0.8) +
  ggtitle("Posterior of residual SD σ_j") +
  scale_y_discrete(labels = trait_labels)

# -------- 3) Intervals (caterpillars) --------
mcmc_intervals(dr_train, pars = c(pars_alpha, pars_beta)) +
  ggtitle("α and β intervals (80%/95%)")

# -------- 4) Trace & pairs for one trait (e.g., CF = [1]) --------
mcmc_trace(dr_train, pars = c("alpha[1]","beta[1]","sigma[1]")) +
  ggtitle("Trace plots (trait 1 = lnCF)")

mcmc_pairs(
  dr_train,
  pars = c("alpha[1]","beta[1]","sigma[1]"),
  off_diag_args = list(size = 1, alpha = 0.5)
) + ggtitle("Pairs (α, β, σ) for lnCF")

# Optional: rank plots to assess mixing
mcmc_rank_overlay(dr_train, pars = c("alpha[1]","beta[1]","sigma[1]")) +
  ggtitle("Rank histograms (lnCF)")

# Grouped by trait
group_lab <- factor(train_data$obs_trait, levels = 1:4, labels = trait_labels)
ppc_dens_overlay_grouped(
  y    = train_data$y,
  yrep = yrep_mat[1:min(200, nrow(yrep_mat)), ],
  group = group_lab
) + ggtitle("PPC: density by trait")

# Intervals by trait (central tendency across observations)
ppc_intervals_grouped(
  y    = train_data$y,
  yrep = yrep_mat[1:min(200, nrow(yrep_mat)), ],
  x    = as.numeric(group_lab),
  group = group_lab
) + ggtitle("PPC: intervals by trait")



########################
# Prediction
########################
# fossils_df: columns can be CF, CH, LF, LH (mm), Species name; NAs allowed.
build_prediction_data <- function(fossils_df,
                                  volumetric_tbl = NULL,  # data.frame: Species, meanlog, sdlog (on ln scale)
                                  predict_all_pairs = TRUE) {
  
  # Normalize fossil columns (flexible names)
  fossils_df <- fossils_df %>%
    rename_with(~"Species", .cols = matches("^spec(imen)?$", ignore.case = TRUE)) %>%
    rename_with(~"CF", .cols = matches("^cf$|^femur.?circ", ignore.case = TRUE)) %>%
    rename_with(~"CH", .cols = matches("^ch$|^humerus.?circ", ignore.case = TRUE)) %>%
    rename_with(~"LF", .cols = matches("^lf$|^femur.?length", ignore.case = TRUE)) %>%
    rename_with(~"LH", .cols = matches("^lh$|^humerus.?length", ignore.case = TRUE))
  
  fos_long <- fossils_df %>%
    mutate(
      lnCF = log(CF), lnCH = log(CH),
      lnLF = log(LF), lnLH = log(LH)
    ) %>%
    select(Species, lnCF, lnCH, lnLF, lnLH) %>%
    pivot_longer(starts_with("ln"), names_to = "trait", values_to = "val") %>%
    filter(trait %in% trait_levels) %>%
    mutate(trait_id = match(trait, trait_levels)) %>%
    filter(is.finite(val))
  
  spec_levels <- unique(fossils_df$Species)
  fos_long <- mutate(fos_long, spec_id = match(Species, spec_levels))
  
  N_spec <- length(spec_levels)
  
  # Volumetric priors mapping (optional)
  has_vol_prior <- rep(0L, N_spec)
  vol_mu <- rep(0.0, N_spec)
  vol_sd <- rep(1.0, N_spec)
  if (!is.null(volumetric_tbl) && nrow(volumetric_tbl) > 0) {
    vm <- match(spec_levels, volumetric_tbl$Species)
    idx <- which(!is.na(vm))
    if (length(idx) > 0) {
      has_vol_prior[idx] <- 1L
      vol_mu[idx] <- volumetric_tbl$meanlog[vm[idx]]
      vol_sd[idx] <- volumetric_tbl$sdlog[vm[idx]]
    }
  }
  
  # Prediction grid
  if (predict_all_pairs) {
    pred_spec <- rep(seq_len(N_spec), each = length(trait_levels))
    pred_trait <- rep(seq_len(length(trait_levels)), times = N_spec)
  } else {
    # only predict for missing pairs: those not in fos_long
    all_grid <- expand.grid(
      spec_id = seq_len(N_spec),
      trait_id = seq_len(length(trait_levels))
    )
    seen <- fos_long %>% distinct(spec_id, trait_id) %>% mutate(seen = TRUE)
    miss <- all_grid %>%
      left_join(seen, by = c("spec_id","trait_id")) %>%
      filter(is.na(seen))
    pred_spec <- miss$spec_id
    pred_trait <- miss$trait_id
  }
  
  list(
    J = length(trait_levels),
    N_spec = N_spec,
    N_obs  = nrow(fos_long),
    obs_Species = fos_long$spec_id,
    obs_trait    = fos_long$trait_id,
    y            = fos_long$val,
    
    has_mass = rep(0L, N_spec),
    lnM_data = rep(0.0, N_spec),
    
    has_vol_prior = has_vol_prior,
    vol_mu = vol_mu,
    vol_sd = vol_sd,
    
    alpha_hat = alpha_hat,
    beta_hat  = beta_hat,
    sd_alpha  = sd_alpha,
    sd_beta   = sd_beta,
    
    a0 = a0, b0 = b0,
    mu_M = mu_M, tau_M = tau_M,
    
    N_pred = length(pred_spec),
    pred_Species = if (length(pred_spec) == 0) integer(0) else pred_spec,
    pred_trait    = if (length(pred_trait) == 0) integer(0) else pred_trait,
    
    .spec_levels = spec_levels,
    .fos_long    = fos_long
  )
}

summarise_draws_df <- function(x, probs = c(0.025, 0.5, 0.975)) {
  tibble::as_tibble(posterior::summarise_draws(x, mean, sd, ~quantile2(.x, probs))) |>
    rename_with(~sub("q\\.", "q", .x))
}

# =========================
# 6) Example prediction run (replace with your fossils)
# =========================
fossils_df <- tibble::tibble(
  Species = c("Fossil_A","Fossil_B","Fossil_C"),
  CF = c(320, NA, 420),
  CH = c(290, 310, NA),
  LF = c(1350, 1200, NA),
  LH = c(1280, NA, 1100)
)

# Optional volumetric priors (ln-mass mean & sd)
vol_tbl <- tibble::tibble(
  Species = c("Fossil_A"),
  meanlog  = c(12.5),
  sdlog    = c(0.30)
)

pred_data <- build_prediction_data(
  fossils_df,
  volumetric_tbl = vol_tbl,
  predict_all_pairs = TRUE  # also predict bones not preserved
)

fit_pred <- mod$sample(
  data = clean_stan_data(pred_data),
  seed = 2026, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000
)

draws_pred <- fit_pred$draws()

#### TODO! ####
# =========================
#Posterior masses for fossils (g and tonnes)
# =========================
lnM_idx <- grep("^lnM\\[", colnames(draws_pred))
lnM_draws <- as_draws_matrix(draws_pred[, lnM_idx])
colnames(lnM_draws) <- pred_data$.spec_levels

lnM_summ <- apply(lnM_draws, 2, function(v) {
  c(mean = mean(v),
    q025 = quantile(v, 0.025),
    q50  = quantile(v, 0.5),
    q975 = quantile(v, 0.975))
}) |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Species")

mass_summ <- lnM_summ |>
  mutate(
    mean_g = exp(mean),
    lo_g   = exp(q025),
    med_g  = exp(q50),
    hi_g   = exp(q975),
    mean_t = mean_g / 1e6,
    lo_t   = lo_g / 1e6,
    med_t  = med_g / 1e6,
    hi_t   = hi_g / 1e6
  ) |>
  select(Species, mean_g, lo_g, med_g, hi_g, mean_t, lo_t, med_t, hi_t)

print(mass_summ)

# =========================
# 8) Posterior predictive bones (mm) for requested pairs
# =========================
N_pred <- fit_pred$metadata()$data()$N_pred
if (N_pred > 0) {
  ypred_idx <- grep("^y_pred\\[", colnames(draws_pred))
  ypred_draws <- as_draws_matrix(draws_pred[, ypred_idx])
  
  lab <- tibble::tibble(
    k = seq_len(N_pred),
    spec_id = fit_pred$metadata()$data()$pred_Species,
    trait_id = fit_pred$metadata()$data()$pred_trait
  ) |>
    mutate(
      Species = pred_data$.spec_levels[spec_id],
      trait = c("lnCF","lnCH","lnLF","lnLH")[trait_id]
    )
  
  bone_pred <- lapply(seq_len(N_pred), function(k) {
    x <- ypred_draws[, k]
    tibble::tibble(
      k = k,
      mean_mm = exp(mean(x)),
      lo_mm   = exp(quantile(x, 0.025)),
      med_mm  = exp(quantile(x, 0.5)),
      hi_mm   = exp(quantile(x, 0.975))
    )
  }) |>
    bind_rows() |>
    left_join(lab, by = "k") |>
    select(Species, trait, mean_mm, lo_mm, med_mm, hi_mm) |>
    arrange(Species, trait)
  
  print(bone_pred)
}

# =========================
# Notes
# - Trained on MASSTIMATE::extants; logs are NATURAL logs.
# - Priors match manuscript, including Inv-Gamma(a0,b0) on sigma^2.
# - Training: lnM observed where available (tight normal).
# - Prediction: lnM latent; volumetric priors optional per Species.
# =========================
