#############################################
# Alligator EVT analysis
# Bastiaan A. Van Velthoven
#############################################

# Run this script from the main analysis directory.
rm(list = ls())

library(dplyr)
library(ggplot2)
library(ismev)
library(patchwork)
library(readxl)
library(scales)

options(scipen = 999, dplyr.summarise.inform = FALSE)
set.seed(42)

source("scripts/functions.R")

data_dir <- "data"
main_figure_dir <- file.path("figures", "manuscript")
supplement_figure_dir <- file.path("figures", "supplement")
results_dir <- "results"

# Bootstrap sizes 
B_main <- 2000
B_stability <- 500
B_holdout <- 1000
B_gof <- 2000

threshold_cm <- 325
threshold_y <- log(threshold_cm)
target_cm <- c(350, 375, 400)
min_exceedances <- 20
threshold_probabilities <- seq(0.85, 0.97, by = 0.005)
return_m_max <- 1e6
max_holdout <- 10

dir.create(main_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_figure_dir, recursive = TRUE, showWarnings = FALSE)
result_dir <- file.path(results_dir, "Alligators")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

woodward_file <- file.path(data_dir, "alligators_woodward.xlsx")
harvest_file <- file.path(
  data_dir,
  "alligator_harvest_3_lakes_combined.csv"
)

###################
# data preprocessing
#######################
raw <- readxl::read_excel(woodward_file)
check_columns(raw, c("Sex", "TL"), "Woodward workbook")

alligators <- raw |>
  mutate(
    TL = suppressWarnings(as.numeric(TL)),
    sex = tolower(trimws(as.character(Sex))),
    sex = case_when(
      sex %in% c("m", "male") | grepl("^m", sex) ~ "male",
      sex %in% c("f", "female") | grepl("^f", sex) ~ "female",
      TRUE ~ "unknown"
    )
  )
if ("Deform" %in% names(alligators)) {
  alligators <- alligators |> filter(is.na(Deform) | Deform == 0)
}
alligators <- alligators |> filter(is.finite(TL), TL > 0)

males <- alligators |> filter(sex == "male")
y <- log(males$TL)

harvest <- read.csv(harvest_file, stringsAsFactors = FALSE)
check_columns(harvest, c("amu", "year", "length_cm"), "FWC harvest CSV")
harvest <- harvest |>
  mutate(
    amu = suppressWarnings(as.numeric(amu)),
    year = suppressWarnings(as.numeric(year)),
    length_cm = suppressWarnings(as.numeric(length_cm))
  )
harvest <- harvest |>
  filter(
    amu %in% c(722, 723, 724),
    year >= 2000, year <= 2025,
    is.finite(length_cm), length_cm > 0
  )

annual_harvest <- harvest |>
  group_by(year) |>
  summarise(
    n_records = n(),
    # The fitted GPD uses strict exceedances above 325 cm.
    n_gt_325 = sum(length_cm > 325),
    n_ge_350 = sum(length_cm >= 350),
    n_ge_375 = sum(length_cm >= 375),
    n_ge_400 = sum(length_cm >= 400),
    max_cm = max(length_cm),
    .groups = "drop"
  )

data_audit <- data.frame(
  item = c(
    "usable alligator records",
    "male records used in the EVT fit",
    "male exceedances above 325 cm",
    "usable FWC harvest records",
    "harvest years"
  ),
  value = c(
    nrow(alligators),
    nrow(males),
    sum(males$TL > threshold_cm),
    nrow(harvest),
    n_distinct(annual_harvest$year)
  )
)

######
# fit
#####

fit <- fit_gpd_y(y, threshold_y, min_exceedances)
if (is.null(fit)) stop("The main alligator GPD fit failed.")

# fixed seed of 42
stability <- bind_rows(lapply(threshold_probabilities, function(probability) {
  u <- as.numeric(quantile(y, probability, type = 8, names = FALSE))
  fit_u <- fit_gpd_y(y, u, min_exceedances)
  if (is.null(fit_u)) return(NULL)
  boot_u <- parametric_bootstrap(fit_u, B_stability)
  if (!nrow(boot_u)) return(NULL)
  endpoints <- exp(boot_u$endpoint_y)
  finite_endpoints <- endpoints[is.finite(endpoints)]
  me <- mean_excess_row(y, u)
  data.frame(
    threshold_probability = probability,
    threshold_cm = exp(u),
    n_exc = fit_u$n_exc,
    mean_excess = me$mean_excess,
    mean_excess_se = me$se,
    xi = fit_u$xi,
    xi_lower95 = q_finite(boot_u$xi, 0.025),
    xi_upper95 = q_finite(boot_u$xi, 0.975),
    endpoint_cm = exp(fit_u$endpoint_y),
    endpoint_finite_lower95_cm = q_finite(finite_endpoints, 0.025),
    endpoint_finite_upper95_cm = q_finite(finite_endpoints, 0.975),
    proportion_unbounded = mean(boot_u$xi >= 0)
  )
}))

# The Anderson--Darling null distribution is obtained parametrically. 
gof <- goodness_of_fit(y, fit, B_gof)
gof_bootstrap <- gof$draws
observed_ad <- gof$summary$statistic
ad_p <- gof$summary$p_value

# bootstrapping
bootstrap <- nonparametric_bootstrap(
  y, threshold_y, B_main, min_exceedances
)

endpoint_cm <- exp(fit$endpoint_y)
bootstrap$endpoint_cm <- exp(bootstrap$endpoint_y)
cov_fit <- fit$raw$cov
se_xi <- if (
  is.matrix(cov_fit) && all(dim(cov_fit) == c(2, 2)) &&
    is.finite(cov_fit[2, 2]) && cov_fit[2, 2] > 0
) sqrt(cov_fit[2, 2]) else NA_real_

fit_summary <- data.frame(
  population = "Male Woodward alligators",
  n = fit$n,
  threshold_cm = threshold_cm,
  n_exc = fit$n_exc,
  exceedance_fraction = fit$p_u,
  sigma = fit$sigma,
  xi = fit$xi,
  se_xi = se_xi,
  xi_lower95 = q_finite(bootstrap$xi, 0.025),
  xi_upper95 = q_finite(bootstrap$xi, 0.975),
  proportion_unbounded = mean(bootstrap$xi >= 0),
  endpoint_cm = endpoint_cm,
  endpoint_upper_one_sided95_cm = q_endpoint(bootstrap$endpoint_cm, 0.95),
  observed_max_cm = max(males$TL),
  anderson_darling = observed_ad,
  anderson_darling_p = ad_p,
  bootstrap_requested = B_main,
  bootstrap_successful = nrow(bootstrap),
  bootstrap_failed = B_main - nrow(bootstrap)
)

# return levels
m_grid <- unique(round(10^seq(
  log10(ceiling(1 / fit$p_u)), log10(return_m_max), length.out = 240
)))
return_matrix <- return_level_matrix(bootstrap, m_grid, threshold_y, exp)
return_interval <- summarise_curve(return_matrix)
return_curve <- data.frame(
  m = m_grid,
  estimate_cm = exp(return_level_y(m_grid, fit)),
  lower95_cm = return_interval$lower95,
  median_cm = return_interval$median,
  upper95_cm = return_interval$upper95
)

empirical_return <- males |>
  arrange(desc(TL)) |>
  mutate(rank = row_number(), m = (nrow(males) + 1) / rank) |>
  filter(TL > threshold_cm) |>
  transmute(m, TL_cm = TL, rank)

qq <- qq_data(y, fit)

survival_upper_cm <- if (is.finite(endpoint_cm)) {
  endpoint_cm * (1 - 1e-8)
} else {
  1.25 * max(males$TL)
}
survival_cm <- seq(threshold_cm, survival_upper_cm, length.out = 220)
survival_y <- log(survival_cm)
survival_draws <- survival_matrix(bootstrap, survival_y, threshold_y)
survival_interval <- summarise_curve(survival_draws)
survival_curve <- data.frame(
  TL_cm = survival_cm,
  estimate = tail_survival_y(survival_y, fit),
  lower95 = survival_interval$lower95,
  upper95 = survival_interval$upper95
)
empirical_survival <- sort(males$TL[males$TL > threshold_cm]) |>
  (function(z) data.frame(
    TL_cm = z,
    survival = vapply(z, function(value) mean(males$TL > value), numeric(1))
  ))() |>
  filter(survival > 0)

# comparison with FWC records
target_table <- bind_rows(lapply(target_cm, function(target) {
  point_probability <- gpd_conditional_survival(
    log(target) - threshold_y, fit$sigma, fit$xi
  )
  bootstrap_probability <- vapply(seq_len(nrow(bootstrap)), function(i) {
    gpd_conditional_survival(
      log(target) - threshold_y,
      bootstrap$sigma[i],
      bootstrap$xi[i]
    )
  }, numeric(1))
  data.frame(
    target_cm = target,
    probability = point_probability,
    lower95 = q_finite(bootstrap_probability, 0.025),
    upper95 = q_finite(bootstrap_probability, 0.975),
    empirical_probability =
      mean(males$TL >= target) / mean(males$TL > threshold_cm)
  )
}))

annual_validation <- bind_rows(lapply(target_cm, function(target) {
  q_point <- target_table$probability[target_table$target_cm == target]
  q_boot <- vapply(seq_len(nrow(bootstrap)), function(i) {
    gpd_conditional_survival(
      log(target) - threshold_y,
      bootstrap$sigma[i],
      bootstrap$xi[i]
    )
  }, numeric(1))

  annual_base <- annual_harvest |>
    transmute(
      year,
      target_cm = target,
      n_gt_reference = n_gt_325,
      observed = .data[[paste0("n_ge_", target)]],
      expected = n_gt_325 * q_point
    )

  bind_rows(lapply(seq_len(nrow(annual_base)), function(i) {
    predictive <- rbinom(
      length(q_boot), annual_base$n_gt_reference[i], q_boot
    )
    annual_base[i, ] |>
      mutate(
        lower95 = q_endpoint(predictive, 0.025),
        upper95 = q_endpoint(predictive, 0.975)
      )
  }))
}))

validation_summary <- annual_validation |>
  mutate(
    inside_interval = observed >= lower95 & observed <= upper95,
    below_interval = observed < lower95,
    above_interval = observed > upper95
  ) |>
  group_by(target_cm) |>
  summarise(
    n_years = n(),
    inside_interval = sum(inside_interval),
    below_interval = sum(below_interval),
    above_interval = sum(above_interval),
    .groups = "drop"
  ) |>
  bind_rows(
    annual_validation |>
      summarise(
        target_cm = NA_real_,
        n_years = n(),
        inside_interval = sum(observed >= lower95 & observed <= upper95),
        below_interval = sum(observed < lower95),
        above_interval = sum(observed > upper95)
      )
  ) |>
  mutate(summary = if_else(is.na(target_cm), "All panels", "By target"))

############
# leave the top out analysis
##################

ordered_lengths <- sort(males$TL, decreasing = TRUE)
holdout <- bind_rows(lapply(0:max_holdout, function(k) {
  training_lengths <- if (k == 0) {
    ordered_lengths
  } else {
    ordered_lengths[-seq_len(k)]
  }
  heldout_lengths <- if (k == 0) numeric() else ordered_lengths[seq_len(k)]

  fit_k <- fit_gpd_y(log(training_lengths), threshold_y, min_exceedances)
  if (is.null(fit_k)) return(NULL)

  bootstrap_k <- nonparametric_bootstrap(
    log(training_lengths), threshold_y,
    B_holdout, min_exceedances
  )
  if (!nrow(bootstrap_k)) return(NULL)

  endpoints_k <- exp(bootstrap_k$endpoint_y)
  upper_k <- q_endpoint(endpoints_k, 0.95)

  data.frame(
    k_removed = k,
    n_train = length(training_lengths),
    n_exc = fit_k$n_exc,
    xi = fit_k$xi,
    xi_lower95 = q_finite(bootstrap_k$xi, 0.025),
    xi_upper95 = q_finite(bootstrap_k$xi, 0.975),
    proportion_unbounded = mean(bootstrap_k$xi >= 0),
    endpoint_cm = exp(fit_k$endpoint_y),
    endpoint_upper_one_sided95_cm = upper_k,
    largest_training_cm = max(training_lengths),
    largest_heldout_cm = if (k == 0) NA_real_ else max(heldout_lengths),
    heldout_below_upper_bound = if (k == 0) {
      NA
    } else {
      max(heldout_lengths) <= upper_k
    },
    bootstrap_requested = B_holdout,
    bootstrap_successful = nrow(bootstrap_k),
    bootstrap_failed = B_holdout - nrow(bootstrap_k)
  )
}))

##############
# figures
############
figure2 <- ggplot(alligators, aes(TL)) +
  geom_histogram(
    bins = 55, fill = COL$blue, colour = "white", linewidth = 0.2
  ) +
  geom_vline(xintercept = threshold_cm, colour = COL$reference, linetype = 2) +
  scale_x_log10(
    breaks = c(100, 150, 200, 300, 400, 500, 600),
    limits = c(NA, 600),
    labels = label_number(big.mark = ",")
  ) +
  labs(x = "Total length (cm, log scale)", y = "Number of alligators") +
  theme(legend.position = "none")
save_figure(
  figure2, "Figure2_alligator_length_distribution",
  main_figure_dir, 7.6, 5.0
)

y_max <- max(
  500,
  if (is.finite(fit_summary$endpoint_upper_one_sided95_cm)) {
    25 * ceiling(fit_summary$endpoint_upper_one_sided95_cm / 25)
  } else 500
)
p_return <- ggplot(return_curve, aes(m, estimate_cm)) +
  geom_ribbon(
    aes(ymin = pmin(lower95_cm, y_max), ymax = pmin(upper95_cm, y_max)),
    fill = COL$blue, alpha = 0.14
  ) +
  geom_line(colour = COL$blue, linewidth = 1) +
  geom_point(
    data = empirical_return, aes(m, TL_cm),
    shape = 21, fill = "white", colour = COL$blue, size = 1.6
  ) +
  scale_x_log10(
    limits = c(10, return_m_max),
    labels = label_number(big.mark = ",")
  ) +
  coord_cartesian(ylim = c(threshold_cm, y_max)) +
  labs(x = "Number of individuals, m", y = "Return level (cm)")

p_qq <- ggplot(qq, aes(theoretical, empirical)) +
  geom_abline(slope = 1, intercept = 0, colour = COL$reference, linetype = 2) +
  geom_point(colour = COL$blue, alpha = 0.7, size = 1.5) +
  coord_equal() +
  labs(x = "Theoretical GPD excess quantile", y = "Empirical excess quantile")

figure3 <- (p_return | p_qq) + plot_annotation(tag_levels = "a")
save_figure(
  figure3, "Figure3_alligator_return_levels_and_qq",
  main_figure_dir, 12, 5.2
)

figure4_data <- annual_validation |>
  mutate(target = factor(
    paste0("TL ≥ ", target_cm, " cm"),
    levels = paste0("TL ≥ ", target_cm |> unique() |> sort(), " cm")
  ))
figure4 <- ggplot(figure4_data, aes(year)) +
  geom_ribbon(
    aes(ymin = lower95, ymax = upper95, fill = "95% predictive interval"),
    alpha = 0.65
  ) +
  geom_col(
    aes(y = observed, fill = "Observed harvest count"),
    width = 0.75, alpha = 0.72
  ) +
  geom_line(
    aes(y = expected, colour = "Expected from fitted tail"), linewidth = 0.9
  ) +
  geom_point(
    aes(y = expected, colour = "Expected from fitted tail"), size = 1.5
  ) +
  facet_wrap(~target, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Observed harvest count" = COL$grey,
    "95% predictive interval" = COL$blue_band
  )) +
  scale_colour_manual(values = c("Expected from fitted tail" = COL$blue)) +
  scale_x_continuous(breaks = breaks_pretty(6)) +
  scale_y_continuous(breaks = breaks_pretty(5), expand = expansion(c(0, 0.06))) +
  labs(x = "Harvest year", y = "Annual number of harvested alligators") +
  theme(panel.spacing.x = grid::unit(12, "pt"))
save_figure(
  figure4, "Figure4_alligator_harvest_validation",
  main_figure_dir, 11.5, 4.8
)

p_mean_excess <- ggplot(stability, aes(threshold_cm, mean_excess)) +
  geom_ribbon(
    aes(ymin = mean_excess - 1.96 * mean_excess_se,
        ymax = mean_excess + 1.96 * mean_excess_se),
    fill = COL$blue_band, alpha = 0.65
  ) +
  geom_line(colour = COL$blue) +
  geom_point(colour = COL$blue, size = 1.3) +
  geom_vline(xintercept = threshold_cm, colour = COL$orange, linetype = 2) +
  labs(x = "Threshold (cm)", y = "Mean log-excess")

p_shape <- ggplot(stability, aes(threshold_cm, xi)) +
  geom_hline(yintercept = 0, colour = COL$reference, linetype = 2) +
  geom_ribbon(
    aes(ymin = xi_lower95, ymax = xi_upper95),
    fill = COL$blue_band, alpha = 0.65
  ) +
  geom_line(colour = COL$blue) +
  geom_point(colour = COL$blue, size = 1.3) +
  geom_vline(xintercept = threshold_cm, colour = COL$orange, linetype = 2) +
  labs(x = "Threshold (cm)", y = expression("GPD shape " * hat(xi)))

p_endpoint <- stability |>
  filter(is.finite(endpoint_cm)) |>
  ggplot(aes(threshold_cm, endpoint_cm)) +
  geom_ribbon(
    aes(ymin = endpoint_finite_lower95_cm, ymax = endpoint_finite_upper95_cm),
    fill = COL$blue_band, alpha = 0.65
  ) +
  geom_line(colour = COL$blue) +
  geom_point(colour = COL$blue, size = 1.3) +
  geom_vline(xintercept = threshold_cm, colour = COL$orange, linetype = 2) +
  labs(x = "Threshold (cm)", y = "Implied endpoint (cm)")

figure_s1 <- (p_mean_excess | p_shape | p_endpoint) +
  plot_annotation(tag_levels = "A")
save_figure(
  figure_s1, "Figure_S1_alligator_threshold_diagnostics",
  supplement_figure_dir, 13, 4.4
)

n_unbounded <- sum(bootstrap$xi >= 0)
unbounded_label <- sprintf(
  "%.2f%% (%s/%s) unbounded",
  100 * n_unbounded / nrow(bootstrap),
  format(n_unbounded, big.mark = ",", scientific = FALSE),
  format(nrow(bootstrap), big.mark = ",", scientific = FALSE)
)

p_xi <- ggplot(bootstrap, aes(xi)) +
  geom_histogram(bins = 40, fill = COL$grey, colour = "white") +
  geom_vline(xintercept = 0, colour = COL$reference) +
  geom_vline(xintercept = fit$xi, colour = COL$blue, linetype = 2) +
  labs(
    x = expression("Sample-bootstrap shape " * xi), y = "Fits",
    subtitle = unbounded_label
  )

finite_endpoints <- bootstrap |> filter(is.finite(endpoint_cm))
p_endpoint_boot <- ggplot(finite_endpoints, aes(endpoint_cm)) +
  geom_histogram(bins = 40, fill = COL$grey, colour = "white") +
  geom_vline(xintercept = endpoint_cm, colour = COL$blue, linetype = 2) +
  labs(x = "Finite endpoint (cm)", y = "Fits")

p_ad <- ggplot(gof_bootstrap |> filter(is.finite(anderson_darling)),
               aes(anderson_darling)) +
  geom_histogram(bins = 40, fill = COL$grey, colour = "white") +
  geom_vline(xintercept = observed_ad, colour = COL$orange, linewidth = 0.9) +
  labs(
    x = "Anderson-Darling statistic", y = "Simulations",
    subtitle = paste0(
      "Parametric-bootstrap p = ", format(round(ad_p, 3), nsmall = 3)
    )
  )

# Panels A--B use the nonparametric full-sample bootstrap. Panel C uses the
# parametric bootstrap distribution under the fitted GPD null model.
figure_s2 <- (p_xi | p_endpoint_boot | p_ad) +
  plot_annotation(tag_levels = "A")
save_figure(
  figure_s2, "Figure_S2_alligator_fit_and_bootstrap",
  supplement_figure_dir, 13, 4.6
)

######## save all output ##########
write_result(data_audit, result_dir, "data_audit.csv")
write_result(fit_summary, result_dir, "fit_summary.csv")
write_result(
  gof$summary,
  result_dir,
  "goodness_of_fit.csv"
)
write_result(stability, result_dir, "threshold_stability.csv")
write_result(target_table, result_dir, "target_probabilities.csv")
write_result(annual_validation, result_dir, "annual_harvest_validation.csv")
write_result(validation_summary, result_dir, "harvest_validation_summary.csv")
write_result(holdout, result_dir, "top_k_holdout.csv")
write_result(return_curve, result_dir, "return_curve.csv")
write_result(survival_curve, result_dir, "fitted_survival_curve.csv")

results <- list(
  settings = list(
    woodward_file = woodward_file,
    harvest_file = harvest_file,
    threshold_cm = threshold_cm,
    target_cm = target_cm,
    bootstrap_sizes = c(
      main = B_main,
      threshold_stability = B_stability,
      top_k = B_holdout,
      goodness_of_fit = B_gof
    )
  ),
  data = list(alligators = alligators, males = males, harvest = harvest),
  data_audit = data_audit,
  fit = fit,
  fit_summary = fit_summary,
  bootstrap = bootstrap,
  goodness_of_fit = gof,
  threshold_stability = stability,
  qq = qq,
  return_curve = return_curve,
  empirical_return = empirical_return,
  survival_curve = survival_curve,
  target_probabilities = target_table,
  annual_validation = annual_validation,
  validation_summary = validation_summary,
  top_k_holdout = holdout
)
saveRDS(results, file.path(result_dir, "alligator_evt_results.rds"))
writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "sessionInfo.txt")
)

cat("\nALLIGATORS\n")
print(fit_summary)
print(target_table)
print(validation_summary)
invisible(results)
