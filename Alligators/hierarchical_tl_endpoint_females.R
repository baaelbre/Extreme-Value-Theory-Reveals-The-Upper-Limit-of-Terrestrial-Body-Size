############################################################
# Hierarchical EVT for log TL exceedances (two datasets) — FEMALES
# - Sergio (CSV, females) + Allan/Woodward (Excel, females)
# - Hierarchical Weibull-domain POT in Stan
# - TL* -> Mass via female-only allometry (from CSV)
# Author: Bastiaan Aelbrecht
# Date: 16/10/2025
############################################################

## ---------------------------
## 0) Libraries & global setup
## ---------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readxl)
  library(pracma);  library(grid)
  library(rstan);   library(MASS)
  library(bayesplot); library(posterior)
})
options(stringsAsFactors = FALSE)
set.seed(42)
rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

## Figure dir and plot theme
fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

theme_science_polished <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.length= unit(0.20, "cm"),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    plot.margin      = margin(5, 5, 5, 5),
    legend.position  = "right"
  )

## ---------------------------
## 1) CONFIG: thresholds, priors
## ---------------------------
# Threshold quantiles per dataset (tune if female sample sizes differ)
q_s <- 0.985  # Sergio (CSV, females)
q_a <- 0.90  # Allan/Woodward (Excel, females)

# Optional reference lines (set NA to hide)
ref_tl_cm <- 322      # e.g. 350
ref_w_kg  <- 170      # e.g. a known record female mass

# CI level for “upper bound” line
ci_level <- 0.95

# Hyperprior centers (on population μ_y = log TL*, and μ_ξ)
mu_y0   <- log(320);  sd_y0   <- 0.20   # center female TL* around ~3.2 m (tune as needed)
mu_xi0  <- -0.25;     sd_xi0  <- 0.10
hc_scale_y  <- 0.10                     # half-Cauchy scale for τ_y
hc_scale_xi <- 0.10                     # half-Cauchy scale for τ_ξ

## ---------------------------
## 2) Data ingest (FEMALES ONLY)
## ---------------------------
## 2a) SERGIO BR dataset (CSV) — used for EVT & allometry (female-only)
sergio_path <- "Data/CaptureData_Gator_allometry_paper_2024.csv"
dat_s_raw   <- read.csv(sergio_path, stringsAsFactors = FALSE)

dat_s <- dat_s_raw %>%
  mutate(
    SVL   = suppressWarnings(as.numeric(SV.length)),
    TL    = suppressWarnings(as.numeric(Total.Length)),
    WTkg  = suppressWarnings(as.numeric(Weight)),
    Sex   = toupper(trimws(Sex)),
    Deform = case_when(
      is.na(Deformities_Notes) ~ 0L,
      nchar(trimws(Deformities_Notes)) == 0 ~ 0L,
      grepl("^none$", trimws(tolower(Deformities_Notes))) ~ 0L,
      TRUE ~ 1L
    )
  ) %>%
  filter(Sex == "F") %>%   # FEMALES ONLY
  mutate(
    SVL  = ifelse(is.finite(SVL)  & SVL  > 0, SVL,  NA_real_),
    TL   = ifelse(is.finite(TL)   & TL   > 0, TL,   NA_real_),
    WTkg = ifelse(is.finite(WTkg) & WTkg > 0, WTkg, NA_real_)
  )

## 2b) ALLAN WOODWARD dataset (Excel) — used for EVT (female-only)
allan_path <- "Data/experimental_alligator_harvest_woodward.xlsx"
dat_a_raw  <- readxl::read_excel(allan_path)

dat_a <- dat_a_raw %>%
  mutate(
    Sex = tolower(trimws(as.character(Sex))),
    Sex = case_when(
      Sex %in% c("f","female","fem","vrouw","f?","female?") ~ "f",
      TRUE ~ Sex
    ),
    SVL   = suppressWarnings(as.numeric(SVL)),
    TL    = suppressWarnings(as.numeric(TL)),
    WTkg  = suppressWarnings(as.numeric(if ("WTkg" %in% names(.)) WTkg else NA_real_)),
    Deform = case_when(is.na(Deform) ~ 0L, TRUE ~ as.integer(Deform))
  ) %>%
  filter(Sex == "f") %>%   # FEMALES ONLY
  mutate(
    TL   = ifelse(is.finite(TL)   & TL   > 0, TL,   NA_real_),
    WTkg = ifelse(is.finite(WTkg) & WTkg > 0, WTkg, NA_real_),
    SVL  = ifelse(is.finite(SVL)  & SVL  > 0, SVL,  NA_real_)
  )

stopifnot(sum(is.finite(dat_s$TL)) >= 50, sum(is.finite(dat_a$TL)) >= 50)

## ---------------------------
## 3) Build log TL series, impute via logSVL if needed (per set)
## ---------------------------
# Regress logTL ~ logSVL on *non-deformed* females within each dataset
reg_s <- dat_s %>% filter(Deform==0L, is.finite(TL), is.finite(SVL)) %>%
  transmute(logTL = log(TL), logSVL = log(SVL))
stopifnot(nrow(reg_s) >= 10)
lm_s <- lm(logTL ~ logSVL, data = reg_s)

dat_s <- dat_s %>%
  mutate(
    logSVL    = ifelse(is.finite(SVL), log(SVL), NA_real_),
    logTL_obs = ifelse(is.finite(TL),  log(TL),  NA_real_),
    logTL_pred= ifelse(is.finite(logSVL),
                       as.numeric(predict(lm_s, newdata = data.frame(logSVL=logSVL))),
                       NA_real_),
    logTL     = ifelse(is.finite(logTL_obs), logTL_obs, logTL_pred)
  )
y_s <- dat_s$logTL[is.finite(dat_s$logTL)]

reg_a <- dat_a %>% filter(Deform==0L, is.finite(TL), is.finite(SVL)) %>%
  transmute(logTL = log(TL), logSVL = log(SVL))
stopifnot(nrow(reg_a) >= 10)
lm_a <- lm(logTL ~ logSVL, data = reg_a)

dat_a <- dat_a %>%
  mutate(
    logSVL    = ifelse(is.finite(SVL), log(SVL), NA_real_),
    logTL_obs = ifelse(is.finite(TL),  log(TL),  NA_real_),
    logTL_pred= ifelse(is.finite(logSVL),
                       as.numeric(predict(lm_a, newdata = data.frame(logSVL=logSVL))),
                       NA_real_),
    logTL     = ifelse(is.finite(logTL_obs), logTL_obs, logTL_pred)
  )
y_a <- dat_a$logTL[is.finite(dat_a$logTL)]

## ---------------------------
## 4) Thresholds & exceedances (Weibull domain)
## ---------------------------
u_s <- as.numeric(quantile(y_s, q_s))
u_a <- as.numeric(quantile(y_a, q_a))
ex_s <- y_s[y_s > u_s]
ex_a <- y_a[y_a > u_a]

stopifnot(length(ex_s) >= 20, length(ex_a) >= 20)

## ---------------------------
## 5) Stan setup (two-group hierarchical Weibull POT)
## ---------------------------
stan_file <- "stan/hier_gp_weibull.stan"  # same as your original, but now fed with female data

N1 <- length(ex_s); y1 <- as.vector(ex_s); u1 <- u_s; ymax1 <- max(ex_s)
N2 <- length(ex_a); y2 <- as.vector(ex_a); u2 <- u_a; ymax2 <- max(ex_a)

stan_data <- list(
  N1 = N1, y1 = y1, u1 = u1, ymax1 = ymax1,
  N2 = N2, y2 = y2, u2 = u2, ymax2 = ymax2,
  mu_y0 = mu_y0,    sd_y0 = sd_y0,
  mu_xi0 = mu_xi0,  sd_xi0 = sd_xi0,
  hc_scale_y = hc_scale_y,
  hc_scale_xi = hc_scale_xi
)

fit <- rstan::stan(
  file   = stan_file,
  data   = stan_data,
  seed   = 42,
  chains = 4,
  iter   = 4000,
  warmup = 1500,
  control = list(adapt_delta = 0.98, max_treedepth = 12)
)

print(fit, pars = c("ystar1","ystar2","xi1","xi2","mu_y","tau_y","mu_xi","tau_xi",
                    "TLstar1_cm","TLstar2_cm","TLstar_pop_cm"),
      probs = c(0.05,0.5,0.95))

## ---------------------------
## 6) Diagnostics: R-hat, ESS, traces, NUTS
## ---------------------------
color_scheme_set("brightblue")
params_main    <- c("ystar1","ystar2","xi1","xi2","mu_y","tau_y","mu_xi","tau_xi")
params_derived <- c("TLstar1_cm","TLstar2_cm","TLstar_pop_cm")

drws <- posterior::as_draws_array(fit)
summ_tbl <- posterior::summarise_draws(
  drws,
  mean, sd, ~quantile(.x, probs = c(0.05,0.5,0.95)),
  rhat, ess_bulk, ess_tail
)
diag_tbl <- subset(summ_tbl, variable %in% c(params_main, params_derived))
write.csv(as.data.frame(diag_tbl), file.path(fig_dir, "stan_diagnostics_female.csv"), row.names = FALSE)

posterior_array <- rstan::extract(fit, permuted = FALSE)
p_trace_core <- bayesplot::mcmc_trace(posterior_array, pars = params_main, facet_args = list(ncol = 2))
ggsave(file.path(fig_dir, "trace_core_params_female.png"), p_trace_core, dpi = 300, w = 12, h = 8, units = "in")

np <- bayesplot::nuts_params(fit)
p_energy <- bayesplot::mcmc_nuts_energy(np)
ggsave(file.path(fig_dir, "nuts_energy_female.png"), p_energy, dpi = 300, w = 8, h = 5, units = "in")

## ---------------------------
## 7) Extract draws
## ---------------------------
post <- as.data.frame(fit)
ystar1_draw <- post$ystar1; ystar2_draw <- post$ystar2
xi1_draw    <- post$xi1;    xi2_draw    <- post$xi2
mu_y_draw   <- post$mu_y;   tau_y_draw  <- post$tau_y
mu_xi_draw  <- post$mu_xi;  tau_xi_draw <- post$tau_xi

TLstar1_cm_draw     <- exp(ystar1_draw)      # Sergio (females)
TLstar2_cm_draw     <- exp(ystar2_draw)      # Allan  (females)
TLstar_pop_cm_draw  <- exp(mu_y_draw)        # population μ_y (females)

## ---------------------------
## 8) Helper density & MAP
## ---------------------------
kde_df <- function(x, from=200, to=600, n=5012) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
mode_from_kde <- function(x, from=200, to=600) {
  d <- kde_df(x, from, to)
  d$t[ which.max(d$d) ]
}

## ---------------------------
## 9) Posterior densities (Sergio/Allan/Population TL*) — FEMALES
## ---------------------------
dens_s <- kde_df(TLstar1_cm_draw, from = 200, to = 600)
dens_a <- kde_df(TLstar2_cm_draw, from = 200, to = 600)
dens_p <- kde_df(TLstar_pop_cm_draw, from = 200, to = 600)

sergio_map <- mode_from_kde(TLstar1_cm_draw, 200, 600)
allan_map  <- mode_from_kde(TLstar2_cm_draw, 200, 600)
pop_map    <- mode_from_kde(TLstar_pop_cm_draw, 200, 600)

sergio_mean   <- mean(TLstar1_cm_draw); sergio_median <- median(TLstar1_cm_draw)
allan_mean    <- mean(TLstar2_cm_draw); allan_median  <- median(TLstar2_cm_draw)
pop_mean      <- mean(TLstar_pop_cm_draw); pop_median <- median(TLstar_pop_cm_draw)

sergio_up  <- as.numeric(quantile(TLstar1_cm_draw, probs = ci_level))
allan_up   <- as.numeric(quantile(TLstar2_cm_draw, probs = ci_level))
pop_up     <- as.numeric(quantile(TLstar_pop_cm_draw, probs = ci_level))

make_single_plot_prior_df <- function(kde_df_obj, prior_df, title_lab, map_val, up_val, map_label=NULL) {
  if (is.null(map_label)) map_label <- sprintf("%.0f cm", map_val)
  ggplot() +
    geom_area(data = prior_df, aes(t, d), fill = "grey70", alpha = 0.18) +
    geom_line(data = prior_df, aes(t, d), color = "grey50", linewidth = 0.8, alpha = 0.9) +
    geom_line(data = kde_df_obj, aes(t, d), color = "darkblue", linewidth = 1.2) +
    geom_vline(xintercept = map_val, color = "purple",  linetype = "dashed",  linewidth = 1.05) +
    geom_vline(xintercept = up_val,  color = "orange",  linetype = "dotdash", linewidth = 1.0) +
    { if (is.finite(ref_tl_cm)) geom_vline(xintercept = ref_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) } +
    annotate("text", x = map_val, y = max(kde_df_obj$d, na.rm = TRUE) * 0.08,
             angle = 90, vjust = -0.4, label = map_label, size = 4.2,
             fontface = "bold", color = "purple") +
    labs(x = "Female TL* (cm)", y = "Density", title = title_lab,
         caption = "Grey: prior on TL* (log-normal). Blue: posterior (Stan).") +
    theme_science_polished + xlim(200,600)
}

t_prior_min_cm <- 200; t_prior_max_cm <- 600; n_prior_points <- 4000
t_prior <- seq(t_prior_min_cm, t_prior_max_cm, length.out = n_prior_points)
prior_df <- data.frame(t = t_prior, d = dlnorm(t_prior, meanlog = mu_y0, sdlog = sd_y0))

p_sergio <- make_single_plot_prior_df(dens_s, prior_df,
                                      "Female TL* — Sergio dataset (Stan)", sergio_map, sergio_up,
                                      sprintf("%.0f cm", sergio_map))
p_sergio
ggsave(file.path(fig_dir, "female_Sergio_TLstar.png"), p_sergio, dpi = 600, w = 7, h = 5)

p_allan <- make_single_plot_prior_df(dens_a, prior_df,
                                     "Female TL* — Allan dataset (Stan)", allan_map, allan_up,
                                     sprintf("%.0f cm", allan_map))
p_allan
ggsave(file.path(fig_dir, "female_Allan_TLstar.png"), p_allan, dpi = 600, w = 7, h = 5)

p_pop <- ggplot() +
  geom_line(data = dens_p, aes(t, d), color = "darkblue", linewidth = 1.2) +
  geom_vline(xintercept = pop_map, color = "purple", linetype = "dashed", linewidth = 1.05) +
  geom_vline(xintercept = pop_up,  color = "orange", linetype = "dotdash", linewidth = 1.0) +
  { if (is.finite(ref_tl_cm)) geom_vline(xintercept = ref_tl_cm, color = "firebrick", linetype = "dashed", linewidth = 1.0) } +
  annotate("text", x = pop_map, y = max(dens_p$d, na.rm = TRUE) * 0.08,
           angle = 90, vjust = -0.4, label = sprintf("%.0f cm", pop_map),
           size = 4.2, fontface = "bold", color = "purple") +
  labs(x = "Female TL* (cm)", y = "Density") +
  theme_science_polished + xlim(200,600)
p_pop
ggsave(file.path(fig_dir, "female_pop_tl_endpoint.png"), p_pop, dpi = 600, w = 7, h = 5)

post_both_df <- bind_rows(mutate(dens_s, dataset = "Sergio"),
                          mutate(dens_a, dataset = "Allan"))
p_both <- ggplot() +
  geom_area(data = prior_df, aes(t, d), fill = "grey75", alpha = 0.22, show.legend = FALSE) +
  geom_line(data = prior_df, aes(t, d), color = "grey50", linewidth = 0.9, show.legend = FALSE) +
  geom_line(data = post_both_df, aes(t, d, color = dataset), linewidth = 1.25, show.legend = FALSE) +
  scale_color_manual(values = c("Sergio"="#1f77b4","Allan"="#2ca02c"), guide = "none") +
  coord_cartesian(xlim = c(200, 600)) +
  labs(x="Female TL* (cm)", y="Density") +
  theme_science_polished + theme(legend.position = "none")
p_both
ggsave(file.path(fig_dir, "female_sergio_allan_tl_endpoint.png"),
       p_both, dpi = 600, w = 7.0, h = 5.0)

cat("\n--- Hierarchical Female TL* (Stan) summaries ---\n")
cat(sprintf("Sergio (F):   MAP = %.1f cm | median = %.1f cm | mean = %.1f cm | %d%% upper = %.1f cm\n",
            sergio_map, sergio_median, sergio_mean, round(100*ci_level), sergio_up))
cat(sprintf("Allan  (F):   MAP = %.1f cm | median = %.1f cm | mean = %.1f cm | %d%% upper = %.1f cm\n",
            allan_map,  allan_median,  allan_mean,  round(100*ci_level), allan_up))
cat(sprintf("Population(F):MAP = %.1f cm | median = %.1f cm | mean = %.1f cm | %d%% upper = %.1f cm\n\n",
            pop_map,    pop_median,    pop_mean,    round(100*ci_level), pop_up))

## =========================================================
## 10) Mass propagation from TL* (Sergio, Allan, Population) — FEMALES
## =========================================================
# Fit female-only allometry on CSV females: log(W) = a + b * log(TL)
w_fit_df_s <- dat_s %>%
  filter(is.finite(WTkg), is.finite(TL)) %>%
  transmute(logW = log(WTkg), logTL = log(TL))
stopifnot(nrow(w_fit_df_s) >= 10)

lm_w_tl_s <- lm(logW ~ logTL, data = w_fit_df_s)
beta_s    <- coef(lm_w_tl_s)
Vbeta_s   <- vcov(lm_w_tl_s)
sigma_e_s <- summary(lm_w_tl_s)$sigma

mass_from_tlstar <- function(tlstar_cm_draw, beta_mu, beta_V, sigma_e,
                             ndraw = 20000, seed = 42) {
  set.seed(seed)
  idx <- sample.int(length(tlstar_cm_draw), size = min(ndraw, length(tlstar_cm_draw)), replace = TRUE)
  tl_log <- log(tlstar_cm_draw[idx])
  betas  <- MASS::mvrnorm(n = length(idx), mu = beta_mu, Sigma = beta_V)
  eps    <- rnorm(length(idx), mean = 0, sd = sigma_e)
  logM   <- betas[,1] + betas[,2] * tl_log + eps
  exp(logM)  # kg
}

mass_sergio <- mass_from_tlstar(TLstar1_cm_draw, beta_s, Vbeta_s, sigma_e_s)
mass_allan  <- mass_from_tlstar(TLstar2_cm_draw, beta_s, Vbeta_s, sigma_e_s)

set.seed(42)
nd   <- min(20000, length(mu_y_draw))
idp  <- sample.int(length(mu_y_draw), size = nd, replace = TRUE)
betp <- MASS::mvrnorm(n = nd, mu = beta_s, Sigma = Vbeta_s)
epsp <- rnorm(nd, 0, sigma_e_s)
logM_pop <- betp[,1] + betp[,2] * mu_y_draw[idp] + epsp   # μ_y already log(TL*)
mass_pop <- exp(logM_pop)

kde_df_mass <- function(x, from = 20, to = 300, n = 4096) {
  d <- density(x, from = from, to = to, n = n, kernel = "gaussian", adjust = 1.0)
  data.frame(t = d$x, d = d$y)
}
mass_mode_kg <- function(x, from = 20, to = 300) { d <- kde_df_mass(x, from, to); d$t[which.max(d$d)] }
summ_mass <- function(x, name, ci = 0.95) {
  data.frame(
    which   = name,
    MAP     = mass_mode_kg(x),
    mean    = mean(x),
    median  = median(x),
    up95    = as.numeric(quantile(x, probs = ci))
  )
}

sum_sergio <- summ_mass(mass_sergio, "Sergio (F)")
sum_allan  <- summ_mass(mass_allan,  "Allan (F)")
sum_pop    <- summ_mass(mass_pop,    "Population (F)")
print(rbind(sum_sergio, sum_allan, sum_pop), row.names = FALSE)

cat("\n--- Female mass endpoints (predictive; TL* -> mass via female allometry) ---\n")
cat(sprintf("Sergio (F):   MAP = %.1f kg | median = %.1f kg | mean = %.1f kg | %d%% upper = %.1f kg\n",
            sum_sergio$MAP, sum_sergio$median, sum_sergio$mean, round(100*ci_level), sum_sergio$up95))
cat(sprintf("Allan  (F):   MAP = %.1f kg | median = %.1f kg | mean = %.1f kg | %d%% upper = %.1f kg\n",
            sum_allan$MAP,  sum_allan$median,  sum_allan$mean,  round(100*ci_level), sum_allan$up95))
cat(sprintf("Population(F):MAP = %.1f kg | median = %.1f kg | mean = %.1f kg | %d%% upper = %.1f kg\n\n",
            sum_pop$MAP,    sum_pop$median,    sum_pop$mean,    round(100*ci_level), sum_pop$up95))

## 10b) Mass plots (95% upper; no legends; optional reference mass line)
dens_m_sergio <- kde_df_mass(mass_sergio, from = 20, to = 300)
dens_m_allan  <- kde_df_mass(mass_allan,  from = 20, to = 300)
dens_m_pop    <- kde_df_mass(mass_pop,    from = 20, to = 300)

mass_map_sergio <- mass_mode_kg(mass_sergio, 20, 300)
mass_map_allan  <- mass_mode_kg(mass_allan,  20, 300)
mass_map_pop    <- mass_mode_kg(mass_pop,    20, 300)

mass_up95_sergio <- as.numeric(quantile(mass_sergio, probs = ci_level))
mass_up95_allan  <- as.numeric(quantile(mass_allan,  probs = ci_level))
mass_up95_pop    <- as.numeric(quantile(mass_pop,    probs = ci_level))

col_sergio <- "#1f77b4"; col_allan <- "#2ca02c"; col_pop <- "darkblue"

p_mass_sergio <- ggplot(dens_m_sergio, aes(t, d)) +
  geom_line(linewidth = 1.2, color = col_sergio) +
  geom_vline(xintercept = mass_map_sergio,  color = "purple",   linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = mass_up95_sergio, color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  { if (is.finite(ref_w_kg)) geom_vline(xintercept = ref_w_kg, color = "firebrick", linetype = "dashed", linewidth = 1.0) } +
  annotate("text", x = mass_map_sergio, y = max(dens_m_sergio$d)*0.10,
           label = sprintf("%.0f kg", mass_map_sergio),
           angle = 90, vjust = -0.5, size = 4.0, fontface = "bold", color = "purple") +
  labs(x = "Female Weight* (kg)", y = "Density", title = "Mass endpoint — Sergio (F, predictive at TL*)") +
  coord_cartesian(xlim = c(20, 300)) +
  theme_science_polished + theme(legend.position = "none")
ggsave(file.path(fig_dir, "female_mass_endpoint_sergio.png"), p_mass_sergio, dpi = 600, w = 7, h = 5)

p_mass_allan <- ggplot(dens_m_allan, aes(t, d)) +
  geom_line(linewidth = 1.2, color = col_allan) +
  geom_vline(xintercept = mass_map_allan,   color = "purple",   linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = mass_up95_allan,  color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  { if (is.finite(ref_w_kg)) geom_vline(xintercept = ref_w_kg, color = "firebrick", linetype = "dashed", linewidth = 1.0) } +
  annotate("text", x = mass_map_allan, y = max(dens_m_allan$d)*0.10,
           label = sprintf("%.0f kg", mass_map_allan),
           angle = 90, vjust = -0.5, size = 4.0, fontface = "bold", color = "purple") +
  labs(x = "Female Weight* (kg)", y = "Density", title = "Mass endpoint — Allan (F, predictive at TL*)") +
  coord_cartesian(xlim = c(20, 300)) +
  theme_science_polished + theme(legend.position = "none")
ggsave(file.path(fig_dir, "female_mass_endpoint_allan.png"), p_mass_allan, dpi = 600, w = 7, h = 5)

p_mass_pop <- ggplot(dens_m_pop, aes(t, d)) +
  geom_line(linewidth = 1.2, color = col_pop) +
  geom_vline(xintercept = mass_map_pop,     color = "purple",   linetype = "dashed",  linewidth = 1.0) +
  geom_vline(xintercept = mass_up95_pop,    color = "orange",   linetype = "dotdash", linewidth = 1.0) +
  { if (is.finite(ref_w_kg)) geom_vline(xintercept = ref_w_kg, color = "firebrick", linetype = "dashed", linewidth = 1.0) } +
  annotate("text", x = mass_map_pop, y = max(dens_m_pop$d)*0.10,
           label = sprintf("%.0f kg", mass_map_pop),
           angle = 90, vjust = -0.5, size = 4.0, fontface = "bold", color = "purple") +
  labs(x = "Female Weight* (kg)", y = "Density") +
  coord_cartesian(xlim = c(20, 300)) +
  theme_science_polished + theme(legend.position = "none")
ggsave(file.path(fig_dir, "female_mass_endpoint_population.png"), p_mass_pop, dpi = 600, w = 7, h = 5)

both_mass_df <- bind_rows(mutate(dens_m_sergio, set = "Sergio (F)"),
                          mutate(dens_m_allan,  set = "Allan (F)"))
p_mass_both <- ggplot(both_mass_df, aes(t, d, color = set)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Sergio (F)" = col_sergio, "Allan (F)" = col_allan)) +
  labs(x = "Female Weight* (kg)", y = "Density", color = NULL) +
  coord_cartesian(xlim = c(20, 300)) +
  theme_science_polished + theme(legend.position = "none")
ggsave(file.path(fig_dir, "female_mass_endpoint_sergio_allan.png"),
       p_mass_both, dpi = 600, w = 8, h = 5.5)

## =========================================================
## 11) Tail probabilities (Stan draws) — dataset & population (FEMALES)
## =========================================================
S_weibull_draws <- function(target_cm, u_log, ystar_draw, xi_draw) {
  y_log <- log(target_cm)
  ratio <- (ystar_draw - y_log) / (ystar_draw - u_log)
  out   <- ifelse(ystar_draw > y_log & is.finite(ratio) & ratio > 0,
                  ratio ^ (-1/xi_draw), 0.0)
  pmin(pmax(out, 0), 1)
}

tail_eval_dataset_stan <- function(target_cm, u_log, ystar_draw, xi_draw, Nu, N,
                                   qlo = 0.05, qhi = 0.95) {
  SW <- S_weibull_draws(target_cm, u_log, ystar_draw, xi_draw)
  data.frame(
    target_cm = target_cm,
    p_cond    = mean(SW),
    p_cond_lo = quantile(SW, probs = qlo),
    p_cond_hi = quantile(SW, probs = qhi),
    p_full    = (Nu/N) * mean(SW),
    p_full_lo = (Nu/N) * quantile(SW, probs = qlo),
    p_full_hi = (Nu/N) * quantile(SW, probs = qhi)
  )
}

rtnorm1 <- function(mu, sd, a, b) {
  u <- runif(1, pnorm(a, mu, sd), pnorm(b, mu, sd))
  qnorm(u, mu, sd)
}

tail_eval_population_stan <- function(target_cm, q_pop = max(q_s, q_a),
                                      qlo = 0.05, qhi = 0.95) {
  y_pool <- c(y_s, y_a); N <- length(y_pool)
  u_pop  <- as.numeric(quantile(y_pool, q_pop))
  Nu     <- sum(y_pool > u_pop)
  
  B <- length(mu_y_draw)
  idx <- sample.int(B, size = min(20000, B), replace = TRUE)
  ystar_pop <- rnorm(length(idx), mean = mu_y_draw[idx], sd = pmax(tau_y_draw[idx], 1e-8))
  xi_pop    <- vapply(seq_along(idx),
                      function(i) rtnorm1(mu_xi_draw[idx[i]], pmax(tau_xi_draw[idx[i]], 1e-8), -1, 0),
                      numeric(1))
  SW <- S_weibull_draws(target_cm, u_log = u_pop, ystar_draw = ystar_pop, xi_draw = xi_pop)
  
  data.frame(
    target_cm = target_cm,
    p_cond    = mean(SW),
    p_cond_lo = quantile(SW, probs = qlo),
    p_cond_hi = quantile(SW, probs = qhi),
    p_full    = (Nu/N) * mean(SW),
    p_full_lo = (Nu/N) * quantile(SW, probs = qlo),
    p_full_hi = (Nu/N) * quantile(SW, probs = qhi)
  )
}

targets <- c(300, 325, 350, 375)  # female-relevant TLs; change as desired
res_sergio <- do.call(rbind, lapply(targets, function(t)
  tail_eval_dataset_stan(t, u_log = u_s, ystar_draw = ystar1_draw, xi_draw = xi1_draw,
                         Nu = length(ex_s), N = length(y_s))))
res_allan  <- do.call(rbind, lapply(targets, function(t)
  tail_eval_dataset_stan(t, u_log = u_a, ystar_draw = ystar2_draw, xi_draw = xi2_draw,
                         Nu = length(ex_a), N = length(y_a))))
res_pop    <- do.call(rbind, lapply(targets, function(t)
  tail_eval_population_stan(t)))

print(res_sergio); print(res_allan); print(res_pop)

if (is.finite(ref_tl_cm)) {
  cat(sprintf("\nSergio (F): at %.0f cm  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n",
              ref_tl_cm,
              subset(res_sergio, target_cm==ref_tl_cm)$p_cond,
              subset(res_sergio, target_cm==ref_tl_cm)$p_cond_lo,
              subset(res_sergio, target_cm==ref_tl_cm)$p_cond_hi,
              subset(res_sergio, target_cm==ref_tl_cm)$p_full))
  cat(sprintf("Allan  (F): at %.0f cm  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n",
              ref_tl_cm,
              subset(res_allan, target_cm==ref_tl_cm)$p_cond,
              subset(res_allan, target_cm==ref_tl_cm)$p_cond_lo,
              subset(res_allan, target_cm==ref_tl_cm)$p_cond_hi,
              subset(res_allan, target_cm==ref_tl_cm)$p_full))
  cat(sprintf("Popul. (F): at %.0f cm  p_cond=%.3e [%.3e, %.3e] ; p_full=%.3e\n\n",
              ref_tl_cm,
              subset(res_pop, target_cm==ref_tl_cm)$p_cond,
              subset(res_pop, target_cm==ref_tl_cm)$p_cond_lo,
              subset(res_pop, target_cm==ref_tl_cm)$p_cond_hi,
              subset(res_pop, target_cm==ref_tl_cm)$p_full))
}
