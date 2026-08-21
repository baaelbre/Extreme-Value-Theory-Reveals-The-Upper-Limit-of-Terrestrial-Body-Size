#############################################
# Mammal EVT analysis
# Bastiaan A. Van Velthoven
#############################################

# Run this script from the main analysis directory.
rm(list = ls())

library(dplyr)
library(ggplot2)
library(ismev)
library(patchwork)
library(scales)

options(scipen = 999, dplyr.summarise.inform = FALSE)
set.seed(42)

source("scripts/functions.R")
data_dir <- "data"
main_figure_dir <- file.path("figures", "manuscript")
supplement_figure_dir <- file.path("figures", "supplement")
results_dir <- "results"
B_main <- 2000
B_stability <- 500
B_rarefaction <- 10000
B_gof <- 2000

threshold_kg <- 50
threshold_y <- log(threshold_kg * 1000)  # natural-log grams
min_exceedances <- 8
threshold_probabilities <- seq(0.70, 0.99, by = 0.01)
diagnostic_max_kg <- 200
richness_targets <- c(sauropod_scale = 200, full_fauna = 4530)
return_m_max <- 50000

dir.create(main_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_figure_dir, recursive = TRUE, showWarnings = FALSE)
result_dir <- file.path(results_dir, "Mammals")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

data_file <- file.path(data_dir, "PHYLACINE_1.2.1_Trait_data.csv")

raw <- read.csv(
  data_file, stringsAsFactors = FALSE, check.names = FALSE,
  fileEncoding = "UTF-8"
)
required <- c(
  "Binomial.1.2", "Order.1.2", "Family.1.2", "Terrestrial",
  "Marine", "Aerial", "Mass.g", "Mass.Method", "IUCN.Status.1.2"
)
check_columns(raw, required, "PHYLACINE trait table")

full <- raw |>
  filter(
    Terrestrial == 1, Marine == 0, Aerial == 0,
    is.finite(Mass.g), Mass.g > 0
  ) |>
  transmute(
    species = gsub("_", " ", Binomial.1.2),
    order = as.character(Order.1.2),
    family = as.character(Family.1.2),
    mass_method = trimws(as.character(Mass.Method)),
    # PHYLACINE labels phylogenetically filled body masses exactly as
    # "Imputed" in Mass.Method. Other allometric estimates are retained.
    mass_is_imputed = !is.na(mass_method) &
      tolower(mass_method) == "imputed",
    status = as.character(IUCN.Status.1.2),
    extinct = status %in% c("EP", "EX", "EW"),
    mass_g = as.numeric(Mass.g),
    mass_kg = mass_g / 1000,
    mass_t = mass_g / 1e6,
    y = log(mass_g)
  ) |>
  arrange(mass_g)
extant <- full |> filter(!extinct)

datasets <- list(
  "Late-Quaternary fauna" = full,
  "Extant fauna" = extant
)

fits <- lapply(datasets, function(data) {
  fit_gpd_y(data$y, threshold_y, min_exceedances)
})
if (any(vapply(fits, is.null, logical(1)))) stop("A main mammal fit failed.")

full_fit <- fits[["Late-Quaternary fauna"]]
gof <- goodness_of_fit(full$y, full_fit, B_gof)

bootstraps <- lapply(datasets, function(data) {
  nonparametric_bootstrap(
    data$y, threshold_y, B_main, min_exceedances
  )
})
if (any(!vapply(bootstraps, nrow, integer(1)))) {
  stop("A mammal taxon bootstrap produced no successful fits.")
}

full_boot <- bootstraps[["Late-Quaternary fauna"]]

fit_summary <- bind_rows(lapply(names(datasets), function(label) {
  fit <- fits[[label]]
  boot <- bootstraps[[label]]
  data.frame(
    dataset = label,
    n = fit$n,
    threshold_kg = threshold_kg,
    n_exc = fit$n_exc,
    exceedance_fraction = fit$p_u,
    sigma = fit$sigma,
    xi = fit$xi,
    xi_lower95 = q_finite(boot$xi, 0.025),
    xi_upper95 = q_finite(boot$xi, 0.975),
    proportion_unbounded = mean(boot$xi >= 0),
    endpoint_t = exp(fit$endpoint_y) / 1e6,
    endpoint_upper_one_sided95_t = q_endpoint(exp(boot$endpoint_y) / 1e6, 0.95),
    observed_max_t = max(datasets[[label]]$mass_t),
    anderson_darling = if (label == "Late-Quaternary fauna") {
      gof$summary$statistic
    } else NA_real_,
    anderson_darling_p = if (label == "Late-Quaternary fauna") {
      gof$summary$p_value
    } else NA_real_,
    bootstrap_requested = B_main,
    bootstrap_successful = nrow(boot),
    bootstrap_failed = B_main - nrow(boot)
  )
}))

data_audit <- data.frame(
  item = c(
    "terrestrial late-Quaternary species",
    "late-Quaternary exceedances above 50 kg",
    "extant terrestrial species",
    "extant exceedances above 50 kg",
    "phylogenetically imputed masses in the full dataset"
  ),
  value = c(
    nrow(full),
    full_fit$n_exc,
    nrow(extant),
    fits[["Extant fauna"]]$n_exc,
    sum(full$mass_is_imputed)
  )
)
####### return levels (richness )
return_curves <- bind_rows(lapply(names(datasets), function(label) {
  fit <- fits[[label]]
  boot <- bootstraps[[label]]
  m <- unique(round(10^seq(
    log10(ceiling(1 / fit$p_u)), log10(return_m_max), length.out = 190
  )))
  matrix <- return_level_matrix(
    boot, m, threshold_y, function(value) exp(value) / 1e6
  )
  interval <- summarise_curve(matrix)
  data.frame(
    dataset = label,
    m = m,
    estimate_t = exp(return_level_y(m, fit)) / 1e6,
    lower95_t = interval$lower95,
    median_t = interval$median,
    upper95_t = interval$upper95
  )
}))

empirical_return <- bind_rows(lapply(names(datasets), function(label) {
  data <- datasets[[label]]
  data |>
    arrange(desc(mass_t)) |>
    mutate(rank = row_number(), m = (nrow(data) + 1) / rank) |>
    filter(mass_kg > threshold_kg) |>
    transmute(dataset = label, m, mass_t, rank)
}))

richness_return_levels <- bind_rows(lapply(names(datasets), function(label) {
  fit <- fits[[label]]
  boot <- bootstraps[[label]]
  m <- unname(richness_targets)
  matrix <- return_level_matrix(
    boot, m, threshold_y, function(value) exp(value) / 1e6
  )
  interval <- summarise_curve(matrix)
  data.frame(
    dataset = label,
    richness = names(richness_targets),
    m = m,
    return_level_t = exp(return_level_y(m, fit)) / 1e6,
    lower95_t = interval$lower95,
    median_t = interval$median,
    upper95_t = interval$upper95
  )
}))

return_markers <- richness_return_levels |>
  filter(m == richness_targets[["full_fauna"]])

############ threshold diagnostics ############

# As in the alligator analysis, uncertainty is recomputed at every threshold.
threshold_grid_kg <- sort(unique(c(
  as.numeric(quantile(
    full$mass_kg, threshold_probabilities, type = 1, names = FALSE
  )),
  threshold_kg,
  diagnostic_max_kg
)))
threshold_grid_kg <- threshold_grid_kg[
  threshold_grid_kg <= diagnostic_max_kg
]

stability <- bind_rows(lapply(threshold_grid_kg, function(u_kg) {
  u <- log(u_kg * 1000)
  fit_u <- fit_gpd_y(full$y, u, min_exceedances)
  if (is.null(fit_u)) return(NULL)
  boot_u <- parametric_bootstrap(fit_u, B_stability)
  if (!nrow(boot_u)) return(NULL)
  endpoints_t <- exp(boot_u$endpoint_y) / 1e6
  finite_endpoints_t <- endpoints_t[is.finite(endpoints_t)]
  me <- mean_excess_row(full$y, u)
  data.frame(
    threshold_probability = mean(full$mass_kg <= u_kg),
    threshold_kg = u_kg,
    n_exc = fit_u$n_exc,
    mean_excess = me$mean_excess,
    mean_excess_se = me$se,
    xi = fit_u$xi,
    xi_lower95 = q_finite(boot_u$xi, 0.025),
    xi_upper95 = q_finite(boot_u$xi, 0.975),
    endpoint_t = exp(fit_u$endpoint_y) / 1e6,
    endpoint_finite_lower95_t = q_finite(finite_endpoints_t, 0.025),
    endpoint_finite_upper95_t = q_finite(finite_endpoints_t, 0.975),
    proportion_unbounded = mean(boot_u$xi >= 0)
  )
}))

qq <- qq_data(full$y, full_fit)

survival_max_t <- if (is.finite(full_fit$endpoint_y)) {
  exp(full_fit$endpoint_y) / 1e6 * (1 - 1e-8)
} else {
  1.25 * max(full$mass_t)
}
survival_t <- exp(seq(log(threshold_kg / 1000), log(survival_max_t), length.out = 220))
survival_y <- log(survival_t * 1e6)
survival_draws <- survival_matrix(full_boot, survival_y, threshold_y)
survival_interval <- summarise_curve(survival_draws)
survival_curve <- data.frame(
  mass_t = survival_t,
  estimate = tail_survival_y(survival_y, full_fit),
  lower95 = survival_interval$lower95,
  upper95 = survival_interval$upper95
)
empirical_survival <- sort(full$mass_t[full$mass_kg > threshold_kg]) |>
  (function(z) data.frame(
    mass_t = z,
    survival = vapply(z, function(value) mean(full$mass_t > value), numeric(1))
  ))() |>
  filter(survival > 0)

#### some sensitivity checks ; what if we only have 200 mammals, randomly sampled ?

rarefaction_m <- 200
rarefied_maxima <- replicate(
  B_rarefaction,
  max(sample(full$mass_t, rarefaction_m, replace = FALSE))
)
rarefaction_summary <- data.frame(
  m = rarefaction_m,
  mean_t = mean(rarefied_maxima),
  median_t = median(rarefied_maxima),
  lower95_t = q_finite(rarefied_maxima, 0.025),
  upper95_t = q_finite(rarefied_maxima, 0.975)
)

no_imputation <- full |> filter(!mass_is_imputed)
no_imputation_fit <- fit_gpd_y(no_imputation$y, threshold_y, min_exceedances)
no_imputation_boot <- if (!is.null(no_imputation_fit)) {
  nonparametric_bootstrap(
    no_imputation$y, threshold_y, B_main, min_exceedances
  )
} else data.frame()

# what if we leave out the proboscidea?
extant_no_proboscidea <- extant |> filter(order != "Proboscidea")
no_proboscidea_fit <- fit_gpd_y(
  extant_no_proboscidea$y, threshold_y, min_exceedances
)
no_proboscidea_boot <- if (!is.null(no_proboscidea_fit)) {
  nonparametric_bootstrap(
    extant_no_proboscidea$y,
    threshold_y,
    B_main,
    min_exceedances
  )
} else data.frame()

imputation_sensitivity <- bind_rows(
  data.frame(
    analysis = "All PHYLACINE masses",
    n = full_fit$n,
    n_exc = full_fit$n_exc,
    xi = full_fit$xi,
    xi_lower95 = q_finite(full_boot$xi, 0.025),
    xi_upper95 = q_finite(full_boot$xi, 0.975),
    proportion_unbounded = mean(full_boot$xi >= 0),
    endpoint_t = exp(full_fit$endpoint_y) / 1e6,
    endpoint_upper_one_sided95_t = q_endpoint(
      exp(full_boot$endpoint_y) / 1e6, 0.95
    ),
    return_level_m200_t = exp(return_level_y(200, full_fit)) / 1e6,
    lower95_m200_t = q_finite(
      return_level_matrix(
        full_boot, 200, threshold_y, function(value) exp(value) / 1e6
      ), 0.025
    ),
    upper95_m200_t = q_finite(
      return_level_matrix(
        full_boot, 200, threshold_y, function(value) exp(value) / 1e6
      ), 0.975
    ),
    bootstrap_requested = B_main,
    bootstrap_successful = nrow(full_boot),
    bootstrap_failed = B_main - nrow(full_boot)
  ),
  if (!is.null(no_imputation_fit) && nrow(no_imputation_boot)) {
    no_imp_return <- return_level_matrix(
      no_imputation_boot, 200, threshold_y,
      function(value) exp(value) / 1e6
    )
    data.frame(
      analysis = "Exclude Mass.Method == 'Imputed'",
      n = no_imputation_fit$n,
      n_exc = no_imputation_fit$n_exc,
      xi = no_imputation_fit$xi,
      xi_lower95 = q_finite(no_imputation_boot$xi, 0.025),
      xi_upper95 = q_finite(no_imputation_boot$xi, 0.975),
      proportion_unbounded = mean(no_imputation_boot$xi >= 0),
      endpoint_t = exp(no_imputation_fit$endpoint_y) / 1e6,
      endpoint_upper_one_sided95_t = q_endpoint(
        exp(no_imputation_boot$endpoint_y) / 1e6, 0.95
      ),
      return_level_m200_t = exp(return_level_y(200, no_imputation_fit)) / 1e6,
      lower95_m200_t = q_finite(no_imp_return, 0.025),
      upper95_m200_t = q_finite(no_imp_return, 0.975),
      bootstrap_requested = B_main,
      bootstrap_successful = nrow(no_imputation_boot),
      bootstrap_failed = B_main - nrow(no_imputation_boot)
    )
  } else NULL
)

proboscidea_curves <- bind_rows(
  {
    m <- unique(round(10^seq(
      log10(ceiling(1 / fits[["Extant fauna"]]$p_u)),
      log10(return_m_max), length.out = 190
    )))
    data.frame(
      analysis = "All extant species", m = m,
      return_level_t = exp(return_level_y(m, fits[["Extant fauna"]])) / 1e6
    )
  },
  if (!is.null(no_proboscidea_fit)) {
    m <- unique(round(10^seq(
      log10(ceiling(1 / no_proboscidea_fit$p_u)),
      log10(return_m_max), length.out = 190
    )))
    data.frame(
      analysis = "Extant, excluding Proboscidea", m = m,
      return_level_t = exp(return_level_y(m, no_proboscidea_fit)) / 1e6
    )
  } else NULL
)

sensitivity_summary <- bind_rows(
  imputation_sensitivity,
  if (!is.null(no_proboscidea_fit) && nrow(no_proboscidea_boot)) {
    no_proboscidea_return <- return_level_matrix(
      no_proboscidea_boot,
      200,
      threshold_y,
      function(value) exp(value) / 1e6
    )

    data.frame(
      analysis = "Extant, excluding Proboscidea",
      n = no_proboscidea_fit$n,
      n_exc = no_proboscidea_fit$n_exc,
      xi = no_proboscidea_fit$xi,
      xi_lower95 = q_finite(no_proboscidea_boot$xi, 0.025),
      xi_upper95 = q_finite(no_proboscidea_boot$xi, 0.975),
      proportion_unbounded = mean(no_proboscidea_boot$xi >= 0),
      endpoint_t = exp(no_proboscidea_fit$endpoint_y) / 1e6,
      endpoint_upper_one_sided95_t = q_endpoint(
        exp(no_proboscidea_boot$endpoint_y) / 1e6, 0.95
      ),
      return_level_m200_t = exp(return_level_y(200, no_proboscidea_fit)) / 1e6,
      lower95_m200_t = q_finite(no_proboscidea_return, 0.025),
      upper95_m200_t = q_finite(no_proboscidea_return, 0.975),
      bootstrap_requested = B_main,
      bootstrap_successful = nrow(no_proboscidea_boot),
      bootstrap_failed = B_main - nrow(no_proboscidea_boot)
    )
  } else NULL
)


####### figures
p_distribution <- ggplot(full, aes(mass_kg)) +
  geom_histogram(bins = 60, fill = COL$blue, colour = "white") +
  geom_vline(xintercept = threshold_kg, colour = COL$reference, linetype = 2) +
  scale_x_log10(labels = label_number(big.mark = ",")) +
  labs(x = "Characteristic species mass (kg, log scale)", y = "Species")

p_qq <- ggplot(qq, aes(theoretical, empirical)) +
  geom_abline(slope = 1, intercept = 0, colour = COL$reference, linetype = 2) +
  geom_point(colour = COL$blue, alpha = 0.7, size = 1.5) +
  coord_equal() +
  labs(x = "Theoretical GPD excess quantile", y = "Empirical excess quantile")

figure5 <- (p_distribution | p_qq) + plot_annotation(tag_levels = "a")
save_figure(
  figure5, "Figure5_mammal_distribution_and_qq",
  main_figure_dir, 12, 5.2
)

mammal_colours <- c(
  "Late-Quaternary fauna" = COL$blue,
  "Extant fauna" = COL$orange
)
figure6 <- ggplot(
  return_curves,
  aes(m, estimate_t, colour = dataset, fill = dataset)
) +
  geom_ribbon(aes(ymin = lower95_t, ymax = upper95_t), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(
    data = empirical_return,
    aes(m, mass_t, colour = dataset),
    inherit.aes = FALSE, shape = 21, fill = "white", size = 1.25, alpha = 0.55
  ) +
  geom_point(
    data = return_markers,
    aes(m, return_level_t, colour = dataset),
    inherit.aes = FALSE, shape = 23, fill = "white", size = 3.2, stroke = 0.9
  ) +
  scale_colour_manual(values = mammal_colours) +
  scale_fill_manual(values = mammal_colours) +
  scale_x_log10(labels = label_number(big.mark = ",")) +
  scale_y_log10(labels = label_number()) +
  labs(x = "Number of represented species, m", y = "Return level (tonnes)")
save_figure(
  figure6, "Figure6_mammal_return_levels",
  main_figure_dir, 7.8, 5.2
)

stability_plot <- stability |>
  filter(threshold_kg <= diagnostic_max_kg)
stability_x_limits <- c(min(stability_plot$threshold_kg), diagnostic_max_kg)

p_mean_excess <- ggplot(stability_plot, aes(threshold_kg, mean_excess)) +
  geom_ribbon(
    aes(ymin = mean_excess - 1.96 * mean_excess_se,
        ymax = mean_excess + 1.96 * mean_excess_se),
    fill = COL$blue_band, alpha = 0.65
  ) +
  geom_line(colour = COL$blue) +
  geom_point(colour = COL$blue, size = 1.3) +
  geom_vline(xintercept = threshold_kg, colour = COL$orange, linetype = 2) +
  scale_x_log10(limits = stability_x_limits, labels = label_number()) +
  labs(x = "Threshold (kg, log scale)", y = "Mean log-excess")

p_shape <- ggplot(stability_plot, aes(threshold_kg, xi)) +
  geom_hline(yintercept = 0, colour = COL$reference, linetype = 2) +
  geom_ribbon(
    aes(ymin = xi_lower95, ymax = xi_upper95),
    fill = COL$blue_band, alpha = 0.65
  ) +
  geom_line(colour = COL$blue) +
  geom_point(colour = COL$blue, size = 1.3) +
  geom_vline(xintercept = threshold_kg, colour = COL$orange, linetype = 2) +
  scale_x_log10(limits = stability_x_limits, labels = label_number()) +
  labs(x = "Threshold (kg, log scale)", y = expression("GPD shape " * hat(xi)))

p_endpoint <- stability_plot |>
  filter(is.finite(endpoint_t)) |>
  ggplot(aes(threshold_kg, endpoint_t)) +
  geom_ribbon(
    aes(
      ymin = endpoint_finite_lower95_t,
      ymax = endpoint_finite_upper95_t
    ),
    fill = COL$blue_band, alpha = 0.65
  ) +
  geom_line(colour = COL$blue) +
  geom_point(colour = COL$blue, size = 1.3) +
  geom_vline(xintercept = threshold_kg, colour = COL$orange, linetype = 2) +
  scale_x_log10(limits = stability_x_limits, labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  labs(x = "Threshold (kg, log scale)", y = "Implied endpoint (tonnes)")

figure_s3 <- (p_mean_excess | p_shape | p_endpoint) +
  plot_annotation(tag_levels = "A")
save_figure(
  figure_s3, "Figure_S3_mammal_threshold_diagnostics",
  supplement_figure_dir, 13, 4.5
)

full_boot_plot <- full_boot |>
  mutate(endpoint_t = exp(endpoint_y) / 1e6)
p_xi <- ggplot(full_boot_plot, aes(xi)) +
  geom_histogram(bins = 40, fill = COL$grey, colour = "white") +
  geom_vline(xintercept = 0, colour = COL$reference) +
  geom_vline(xintercept = full_fit$xi, colour = COL$blue, linetype = 2) +
  labs(
    x = expression("Taxon-bootstrap shape " * xi), y = "Fits",
    subtitle = paste0(round(100 * mean(full_boot$xi >= 0), 1), "% unbounded")
  )

p_endpoint_boot <- full_boot_plot |>
  filter(is.finite(endpoint_t)) |>
  ggplot(aes(endpoint_t)) +
  geom_histogram(bins = 40, fill = COL$grey, colour = "white") +
  geom_vline(
    xintercept = exp(full_fit$endpoint_y) / 1e6,
    colour = COL$blue, linetype = 2
  ) +
  scale_x_log10(labels = label_number()) +
  labs(x = "Finite endpoint (tonnes, log scale)", y = "Fits")

p_ad <- gof$draws |>
  filter(is.finite(anderson_darling)) |>
  ggplot(aes(anderson_darling)) +
  geom_histogram(bins = 40, fill = COL$grey, colour = "white") +
  geom_vline(
    xintercept = gof$summary$statistic,
    colour = COL$orange, linewidth = 0.9
  ) +
  labs(
    x = "Anderson-Darling statistic", y = "Simulations",
    subtitle = paste0(
      "Parametric-bootstrap p = ",
      format(round(gof$summary$p_value, 3), nsmall = 3)
    )
  )


figure_s4 <- (p_xi | p_endpoint_boot | p_ad) +
  plot_annotation(tag_levels = "A")
save_figure(
  figure_s4, "Figure_S4_mammal_fit_and_bootstrap",
  supplement_figure_dir, 13, 4.6
)

#######
# save output
######
mass_method_counts <- full |>
  count(mass_method, mass_is_imputed, name = "n_species") |>
  arrange(desc(n_species))

write_result(data_audit, result_dir, "data_audit.csv")
write_result(fit_summary, result_dir, "fit_summary.csv")
write_result(gof$summary, result_dir, "goodness_of_fit.csv")
write_result(richness_return_levels, result_dir, "richness_return_levels.csv")
write_result(return_curves, result_dir, "return_curves.csv")
write_result(stability, result_dir, "threshold_stability.csv")
write_result(sensitivity_summary, result_dir, "sensitivity_summary.csv")
write_result(imputation_sensitivity, result_dir, "imputation_sensitivity.csv")
write_result(mass_method_counts, result_dir, "mass_method_counts.csv")
write_result(proboscidea_curves, result_dir, "proboscidea_return_curves.csv")
write_result(survival_curve, result_dir, "fitted_survival_curve.csv")
write_result(rarefaction_summary, result_dir, "rarefaction_m200.csv")
write_result(
  full |> arrange(desc(mass_t)) |>
    dplyr::select(species, order, family, mass_t, mass_method, status) |>
    slice_head(n = 30),
  result_dir, "top_30_species.csv"
)

results <- list(
  settings = list(
    data_file = data_file,
    threshold_kg = threshold_kg,
    diagnostic_max_kg = diagnostic_max_kg,
    richness_targets = richness_targets,
    bootstrap_sizes = c(
      main = B_main,
      threshold_stability = B_stability,
      rarefaction = B_rarefaction,
      goodness_of_fit = B_gof
    )
  ),
  data = list(full = full, extant = extant),
  data_audit = data_audit,
  mass_method_counts = mass_method_counts,
  fits = fits,
  fit_summary = fit_summary,
  bootstraps = bootstraps,
  goodness_of_fit = gof,
  threshold_stability = stability,
  qq = qq,
  return_curves = return_curves,
  empirical_return = empirical_return,
  richness_return_levels = richness_return_levels,
  rarefied_maxima = rarefied_maxima,
  rarefaction_summary = rarefaction_summary,
  no_imputation_data = no_imputation,
  no_imputation_bootstrap = no_imputation_boot,
  imputation_sensitivity = imputation_sensitivity,
  no_proboscidea_bootstrap = no_proboscidea_boot,
  proboscidea_return_curves = proboscidea_curves,
  sensitivity_summary = sensitivity_summary,
  survival_curve = survival_curve
)
saveRDS(results, file.path(result_dir, "mammal_evt_results.rds"))
writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "sessionInfo.txt")
)

cat("\nMAMMALS\n")
print(fit_summary)
print(richness_return_levels)
invisible(results)
