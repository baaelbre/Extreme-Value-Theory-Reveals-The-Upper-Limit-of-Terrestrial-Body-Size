#############################################
# Sauropod EVT analysis
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
B_holdout <- 1000
B_combined <- 5000
B_gof <- 2000

threshold_t <- 22
threshold_y <- log(threshold_t)  # natural-log tonnes
min_exceedances <- 8
campione_see_log10 <- 0.136
threshold_grid_t <- seq(1, 30, by = 1)
richness_targets <- c(current_record = 200, mammal_scale = 4530)
return_m_min <- 10
return_m_max <- 50000
max_holdout <- 10

dir.create(main_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_figure_dir, recursive = TRUE, showWarnings = FALSE)
result_dir <- file.path(results_dir, "Sauropods")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

data_file <- file.path(
  data_dir,
  "DEmic23_updated_Supplemental_Data_withPubYear_withAvailability.csv"
)
if (!file.exists(data_file)) stop("Missing input file: ", data_file)

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a")] <- NA_character_
  x
}
safe_numeric <- function(x) suppressWarnings(as.numeric(x))
draw_mass <- function(mass_t) {
  10^(log10(mass_t) + rnorm(length(mass_t), 0, campione_see_log10))
}

bootstrap_masses <- function(mass_t, B, mode, threshold = threshold_t) {
  n <- length(mass_t)
  bind_rows(lapply(seq_len(B), function(b) {
    draw <- switch(
      mode,
      sampling = sample(mass_t, n, replace = TRUE),
      reconstruction = draw_mass(mass_t),
      combined = sample(draw_mass(mass_t), n, replace = TRUE),
      stop("Unknown sauropod bootstrap mode: ", mode)
    )
    row <- fit_draw_row(
      fit_gpd_y(log(draw), log(threshold), min_exceedances), b
    )
    if (!is.null(row)) row$uncertainty <- mode
    row
  }))
}
# read data
raw <- read.csv(
  data_file, check.names = FALSE, stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)
names(raw) <- trimws(names(raw))
required <- c(
  "genus and species", "name for phylogenetic tree", "specimen number",
  "taxon", "period", "stage", "Ma", "landmass",
  "Campione est mass quadratic (kg)",
  "quadratic mass est low (kg)", "quadratic mass est high (kg)",
  "source", "notes", "publication year"
)
check_columns(raw, required, "D'Emic sauropod table")

all_records <- raw |>
  transmute(
    species = clean_character(`genus and species`),
    tree_name = clean_character(`name for phylogenetic tree`),
    specimen_number = clean_character(`specimen number`),
    clade = clean_character(taxon),
    period = clean_character(period),
    stage = clean_character(stage),
    age_ma = safe_numeric(Ma),
    landmass = clean_character(landmass),
    mass_kg = safe_numeric(`Campione est mass quadratic (kg)`),
    mass_lower_kg = safe_numeric(`quadratic mass est low (kg)`),
    mass_upper_kg = safe_numeric(`quadratic mass est high (kg)`),
    source = clean_character(source),
    notes = clean_character(notes),
    publication_year = safe_numeric(`publication year`),
    taxon_id = coalesce(tree_name, species),
    mass_t = mass_kg / 1000,
    mass_lower_t = mass_lower_kg / 1000,
    mass_upper_t = mass_upper_kg / 1000
  )

specimens <- all_records |>
  filter(
    !is.na(species), !is.na(taxon_id),
    is.finite(mass_t), mass_t > 0
  )
taxa <- specimens |>
  group_by(taxon_id) |>
  arrange(desc(mass_t), .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  arrange(desc(mass_t))


figure1_data <- taxa |>
  filter(is.finite(publication_year))

argentinosaurus <- figure1_data |>
  filter(grepl("^Argentinosaurus\\b", species, ignore.case = TRUE)) |>
  slice_max(mass_t, n = 1, with_ties = FALSE)
tail_band_upper_t <- argentinosaurus$mass_upper_t[[1]]


top10 <- figure1_data |>
  slice_max(mass_t, n = 10, with_ties = FALSE) |>
  arrange(desc(mass_t)) |>
  mutate(highlight = factor(species, levels = species))

highlight_colours <- c(
  "#E76F51", "#C99A00", "#94A900", "#39AF35", "#2DBE6C",
  "#3CB6C4", "#40A4EC", "#8668E8", "#C346DD", "#E83AAE"
)
names(highlight_colours) <- levels(top10$highlight)
italic_labels <- as.expression(lapply(
  levels(top10$highlight), function(label) bquote(italic(.(label)))
))

figure1_upper_t <- 5 * tail_band_upper_t # to draw the phyliopic silhouettes.

figure1 <- ggplot() +
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = threshold_t, ymax = tail_band_upper_t,
    fill = "grey85", alpha = 0.68
  ) +
  geom_errorbar(
    data = figure1_data,
    aes(
      x = publication_year,
      ymin = mass_lower_t, ymax = mass_upper_t
    ),
    width = 0, colour = "grey78", linewidth = 0.55, alpha = 0.78
  ) +
  geom_point(
    data = figure1_data,
    aes(publication_year, mass_t),
    colour = "grey78", size = 2.25, alpha = 0.78
  ) +
  geom_errorbar(
    data = top10,
    aes(
      x = publication_year,
      ymin = mass_lower_t, ymax = mass_upper_t,
      colour = highlight
    ),
    width = 0, linewidth = 0.75
  ) +
  geom_point(
    data = top10,
    aes(publication_year, mass_t, colour = highlight),
    size = 2.8
  ) +
  scale_colour_manual(
    values = highlight_colours,
    breaks = levels(top10$highlight),
    labels = italic_labels,
    guide = guide_legend(
      ncol = 2, byrow = FALSE,
      override.aes = list(linewidth = 0, size = 3)
    )
  ) +
  scale_x_continuous(
    breaks = seq(1880, 2020, by = 20),
    expand = expansion(mult = c(0.045, 0.045))
  ) +
  scale_y_log10(
    limits = c(0.1, figure1_upper_t),
    breaks = c(0.1, 1, 10, 100),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.04, 0.03))
  ) +
  labs(x = "Publication year", y = "Body mass (tonnes)") +
  theme(
    legend.position = c(0.04, 0.04),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key.height = grid::unit(0.40, "cm"),
    legend.key.width = grid::unit(0.45, "cm"),
    legend.text = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.5)
  )

save_figure(
  figure1, "Figure1_sauropod_mass_by_publication_year",
  main_figure_dir, 8.5, 6.3
)

audit <- data.frame(
  item = c(
    "CSV rows", "valid specimen records", "represented taxa",
    paste0("taxa above ", threshold_t, " tonnes"),
    "largest published mass (tonnes)"
  ),
  value = c(
    nrow(all_records), nrow(specimens), nrow(taxa),
    sum(taxa$mass_t > threshold_t), max(taxa$mass_t)
  )
)

# fit
fit <- fit_gpd_y(log(taxa$mass_t), threshold_y, min_exceedances)
gof <- goodness_of_fit(log(taxa$mass_t), fit, B_gof)

sampling_boot <- bootstrap_masses(
  taxa$mass_t, B_main, "sampling"
)
reconstruction_boot <- bootstrap_masses(
  taxa$mass_t, B_main, "reconstruction"
)
combined_boot <- bootstrap_masses(
  taxa$mass_t, B_combined, "combined"
)


summarise_bootstrap <- function(draws, uncertainty, B_requested) {
  endpoint_t <- exp(draws$endpoint_y)
  data.frame(
    uncertainty = uncertainty,
    n = fit$n,
    threshold_t = threshold_t,
    n_exc = fit$n_exc,
    exceedance_fraction = fit$p_u,
    sigma = fit$sigma,
    xi = fit$xi,
    xi_lower95 = q_finite(draws$xi, 0.025),
    xi_upper95 = q_finite(draws$xi, 0.975),
    proportion_unbounded = mean(draws$xi >= 0),
    endpoint_t = exp(fit$endpoint_y),
    endpoint_upper_one_sided95_t = q_endpoint(endpoint_t, 0.95),
    observed_max_t = max(taxa$mass_t),
    anderson_darling = gof$summary$statistic,
    anderson_darling_p = gof$summary$p_value,
    bootstrap_requested = B_requested,
    bootstrap_successful = nrow(draws),
    bootstrap_failed = B_requested - nrow(draws)
  )
}
fit_summary <- bind_rows(
  summarise_bootstrap(
    sampling_boot,
    "Taxon sampling; published masses fixed",
    B_main
  ),
  summarise_bootstrap(
    reconstruction_boot,
    "Allometric residual uncertainty only",
    B_main
  ),
  summarise_bootstrap(
    combined_boot,
    "Taxon sampling + allometric residual uncertainty",
    B_combined
  )
)
# return levels
make_return_curve <- function(draws) {
  m <- unique(round(10^seq(
    log10(max(return_m_min, ceiling(1 / fit$p_u))),
    log10(return_m_max), length.out = 190
  )))
  matrix <- return_level_matrix(draws, m, threshold_y, exp)
  interval <- summarise_curve(matrix)
  data.frame(
    m = m,
    estimate_t = exp(return_level_y(m, fit)),
    lower95_t = interval$lower95,
    median_t = interval$median,
    upper95_t = interval$upper95
  )
}
return_fixed <- make_return_curve(sampling_boot)
return_combined <- make_return_curve(combined_boot)

empirical_return <- taxa |>
  arrange(desc(mass_t)) |>
  mutate(rank = row_number(), m = (nrow(taxa) + 1) / rank) |>
  filter(mass_t > threshold_t) |>
  transmute(m, mass_t, rank)

richness_return_levels <- bind_rows(lapply(
  list(
    "Taxon sampling; published masses fixed" = sampling_boot,
    "Allometric residual uncertainty only" = reconstruction_boot,
    "Taxon sampling + allometric residual uncertainty" = combined_boot
  ),
  function(draws) {
    uncertainty <- unique(draws$uncertainty)
    m <- unname(richness_targets)
    matrix <- return_level_matrix(draws, m, threshold_y, exp)
    interval <- summarise_curve(matrix)
    data.frame(
      uncertainty = uncertainty,
      richness = names(richness_targets),
      m = m,
      return_level_t = exp(return_level_y(m, fit)),
      lower95_t = interval$lower95,
      median_t = interval$median,
      upper95_t = interval$upper95
    )
  }
))

markers <- richness_return_levels |>
  filter(uncertainty == "sampling") |>
  dplyr::select(richness, m, return_level_t)

survival_max_t <- if (is.finite(fit$endpoint_y)) {
  exp(fit$endpoint_y) * (1 - 1e-8)
} else {
  1.25 * max(taxa$mass_t)
}
survival_t <- exp(seq(log(threshold_t), log(survival_max_t), length.out = 220))
survival_y <- log(survival_t)
survival_draws <- survival_matrix(sampling_boot, survival_y, threshold_y)
survival_interval <- summarise_curve(survival_draws)
survival_curve <- data.frame(
  mass_t = survival_t,
  estimate = tail_survival_y(survival_y, fit),
  lower95 = survival_interval$lower95,
  upper95 = survival_interval$upper95
)
empirical_survival <- sort(taxa$mass_t[taxa$mass_t >= threshold_t]) |>
  (function(z) data.frame(
    mass_t = z,
    survival = vapply(z, function(value) mean(taxa$mass_t > value), numeric(1))
  ))() |>
  filter(survival > 0)

# threshold sensitivity
stability <- bind_rows(lapply(threshold_grid_t, function(u_t) {
  fit_u <- fit_gpd_y(log(taxa$mass_t), log(u_t), min_exceedances)
  if (is.null(fit_u)) return(NULL)
  fixed_u <- bootstrap_masses(
    taxa$mass_t, B_stability, "sampling", u_t
  )
  combined_u <- bootstrap_masses(
    taxa$mass_t, B_stability, "combined", u_t
  )
  if (!nrow(fixed_u) || !nrow(combined_u)) return(NULL)

  fixed_endpoints <- exp(fixed_u$endpoint_y)
  combined_endpoints <- exp(combined_u$endpoint_y)
  me <- mean_excess_row(log(taxa$mass_t), log(u_t))
  data.frame(
    threshold_t = u_t,
    n_exc = fit_u$n_exc,
    mean_excess = me$mean_excess,
    mean_excess_se = me$se,
    xi = fit_u$xi,
    xi_lower95_fixed = q_finite(fixed_u$xi, 0.025),
    xi_upper95_fixed = q_finite(fixed_u$xi, 0.975),
    xi_lower95_combined = q_finite(combined_u$xi, 0.025),
    xi_upper95_combined = q_finite(combined_u$xi, 0.975),
    endpoint_t = exp(fit_u$endpoint_y),
    endpoint_finite_lower95_fixed_t = q_finite(fixed_endpoints, 0.025),
    endpoint_finite_upper95_fixed_t = q_finite(fixed_endpoints, 0.975),
    endpoint_finite_lower95_combined_t = q_finite(combined_endpoints, 0.025),
    endpoint_finite_upper95_combined_t = q_finite(combined_endpoints, 0.975),
    proportion_unbounded_fixed = mean(fixed_u$xi >= 0),
    proportion_unbounded_combined = mean(combined_u$xi >= 0)
  )
}))

qq <- qq_data(log(taxa$mass_t), fit)


# Remove the largest taxon-level estimates in order; ie truncate the tail.

holdout <- bind_rows(lapply(0:max_holdout, function(k) {
  train_mass <- if (k == 0) taxa$mass_t else taxa$mass_t[-seq_len(k)]
  heldout_mass <- if (k == 0) numeric() else taxa$mass_t[seq_len(k)]
  fit_k <- fit_gpd_y(log(train_mass), threshold_y, min_exceedances)
  if (is.null(fit_k)) return(NULL)
  draws_k <- bootstrap_masses(
    train_mass, B_holdout, "combined", threshold_t
  )
  if (!nrow(draws_k)) return(NULL)
  endpoints <- exp(draws_k$endpoint_y)
  finite_endpoints <- endpoints[is.finite(endpoints)]
  data.frame(
    k_removed = k,
    n_train = length(train_mass),
    n_exc = fit_k$n_exc,
    xi = fit_k$xi,
    xi_lower95 = q_finite(draws_k$xi, 0.025),
    xi_upper95 = q_finite(draws_k$xi, 0.975),
    proportion_unbounded = mean(draws_k$xi >= 0),
    endpoint_t = exp(fit_k$endpoint_y),
    endpoint_finite_lower95_t = q_finite(finite_endpoints, 0.025),
    endpoint_finite_upper95_t = q_finite(finite_endpoints, 0.975),
    endpoint_upper_one_sided95_t = q_endpoint(endpoints, 0.95),
    largest_train_t = max(train_mass),
    heldout_max_t = if (k) max(heldout_mass) else NA_real_,
    heldout_below_upper_bound = if (k) {
      max(heldout_mass) <= q_endpoint(endpoints, 0.95)
    } else {
      NA
    },
    bootstrap_requested = B_holdout,
    bootstrap_successful = nrow(draws_k),
    bootstrap_failed = B_holdout - nrow(draws_k)
  )
}))

figure7 <- ggplot() +
  geom_ribbon(
    data = return_combined,
    aes(m, ymin = lower95_t, ymax = upper95_t,
        fill = "Taxon resampling + allometric residual uncertainty"),
    alpha = 0.22, colour = NA
  ) +
  geom_ribbon(
    data = return_fixed,
    aes(m, ymin = lower95_t, ymax = upper95_t,
        fill = "Taxon-resampling uncertainty"),
    alpha = 0.30, colour = NA
  ) +
  geom_line(
    data = return_combined, aes(m, estimate_t),
    colour = COL$blue, linewidth = 1.05
  ) +
  geom_point(
    data = empirical_return, aes(m, mass_t),
    shape = 21, fill = "white", colour = COL$blue,
    size = 1.65, alpha = 0.75
  ) +
  geom_point(
    data = markers |> filter(richness == "current_record"),
    aes(m, return_level_t),
    shape = 23, fill = COL$blue, colour = COL$blue, size = 3.3
  ) +
  geom_point(
    data = markers |> filter(richness == "mammal_scale"),
    aes(m, return_level_t),
    shape = 22, fill = COL$yellow, colour = COL$yellow, size = 3.3
  ) +
  scale_fill_manual(
    values = c(
      "Taxon-resampling uncertainty" = COL$blue,
      "Taxon resampling + allometric residual uncertainty" = COL$yellow
    ),
    breaks = c(
      "Taxon-resampling uncertainty",
      "Taxon resampling + allometric residual uncertainty"
    )
  ) +
  scale_x_log10(
    limits = c(return_m_min, return_m_max),
    labels = label_number(big.mark = ",")
  ) +
  scale_y_continuous(labels = label_number()) +
  labs(x = "Number of represented taxa, m", y = "Return level (tonnes)")
save_figure(
  figure7, "Figure7_sauropod_return_levels",
  main_figure_dir, 7.5, 5.2
)

p_distribution <- ggplot(taxa, aes(mass_t)) +
  geom_histogram(bins = 30, fill = COL$grey, colour = "white") +
  geom_vline(xintercept = threshold_t, colour = COL$orange, linetype = 2) +
  scale_x_log10(labels = label_number()) +
  labs(
    x = "Body mass (tonnes, log scale)",
    y = "Number of represented taxa"
  )

p_mean_excess <- ggplot(stability, aes(threshold_t, mean_excess)) +
  geom_ribbon(
    aes(
      ymin = mean_excess - 1.96 * mean_excess_se,
      ymax = mean_excess + 1.96 * mean_excess_se
    ),
    fill = COL$blue_band, alpha = 0.65
  ) +
  geom_line(colour = COL$blue) +
  geom_point(colour = COL$blue, size = 1.2) +
  geom_vline(
    xintercept = threshold_t,
    colour = COL$orange,
    linetype = 2
  ) +
  scale_x_log10(
    limits = c(1, 30),
    breaks = c(1, 2, 5, 10, 20, 30),
    labels = label_number()
  ) +
  labs(
    x = "Threshold (tonnes, log scale)",
    y = "Mean log-excess"
  )

p_shape <- ggplot(stability, aes(threshold_t, xi)) +
  geom_hline(yintercept = 0, colour = COL$reference, linetype = 2) +
  geom_ribbon(
    aes(ymin = xi_lower95_combined, ymax = xi_upper95_combined,
        fill = "Taxon + allometric"),
    alpha = 0.20
  ) +
  geom_ribbon(
    aes(ymin = xi_lower95_fixed, ymax = xi_upper95_fixed,
        fill = "Taxon sampling"),
    alpha = 0.28
  ) +
  geom_line(colour = COL$blue) + geom_point(colour = COL$blue, size = 1.2) +
  geom_vline(xintercept = threshold_t, colour = COL$orange, linetype = 2) +
  scale_x_continuous(limits = c(1, 30), breaks = c(1, 5, 10, 15, 20, 25, 30)) +
  scale_fill_manual(values = c(
    "Taxon sampling" = COL$blue,
    "Taxon + allometric" = COL$yellow
  )) +
  labs(x = "Threshold (tonnes)", y = expression("GPD shape " * hat(xi)))

p_endpoint <- stability |>
  filter(is.finite(endpoint_t)) |>
  ggplot(aes(threshold_t, endpoint_t)) +
  geom_line(colour = COL$blue) + geom_point(colour = COL$blue, size = 1.2) +
  geom_vline(xintercept = threshold_t, colour = COL$orange, linetype = 2) +
  scale_x_continuous(limits = c(1, 30), breaks = c(1, 5, 10, 15, 20, 25, 30)) +
  scale_y_log10(labels = label_number()) +
  labs(x = "Threshold (tonnes)", y = "Implied endpoint (tonnes, log scale)")

figure_s5 <- ((p_distribution | p_mean_excess) / (p_shape | p_endpoint)) +
  plot_annotation(tag_levels = "A")
save_figure(
  figure_s5, "Figure_S5_sauropod_threshold_diagnostics",
  supplement_figure_dir, 11, 8.3
)

p_qq <- ggplot(qq, aes(theoretical, empirical)) +
  geom_abline(slope = 1, intercept = 0, colour = COL$reference, linetype = 2) +
  geom_point(
    shape = 21, fill = "white", colour = COL$blue,
    stroke = 0.65, size = 1.8
  ) +
  coord_equal() +
  labs(x = "Theoretical GPD excess quantile", y = "Empirical excess quantile")

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

figure_s6 <- (p_qq | p_ad) + plot_annotation(tag_levels = "A")
save_figure(
  figure_s6, "Figure_S6_sauropod_fit_diagnostics",
  supplement_figure_dir, 9.5, 4.6
)

# save output
write_result(audit, result_dir, "data_audit.csv")
write_result(
  taxa |> dplyr::select(
    taxon_id, species, specimen_number, mass_t, publication_year, source
  ),
  result_dir, "analysis_data.csv"
)
write_result(fit_summary, result_dir, "fit_summary.csv")
write_result(gof$summary, result_dir, "goodness_of_fit.csv")
write_result(richness_return_levels, result_dir, "richness_return_levels.csv")
write_result(stability, result_dir, "threshold_stability.csv")
write_result(holdout, result_dir, "top_k_holdout.csv")
write_result(return_fixed, result_dir, "return_curve_taxon_resampling.csv")
write_result(return_combined, result_dir, "return_curve_combined_uncertainty.csv")
write_result(survival_curve, result_dir, "fitted_survival_curve.csv")

results <- list(
  settings = list(
    data_file = data_file,
    threshold_t = threshold_t,
    campione_see_log10 = campione_see_log10,
    richness_targets = richness_targets,
    bootstrap_sizes = c(
      main = B_main,
      threshold_stability = B_stability,
      top_k = B_holdout,
      combined = B_combined,
      goodness_of_fit = B_gof
    ),
    interpretation = paste(
      "One maximum-known mass estimate per represented taxon;",
      "allometric draws are a separate residual-error sensitivity analysis."
    )
  ),
  audit = audit,
  specimen_data = specimens,
  taxon_data = taxa,
  figure1 = list(
    data = figure1_data,
    highlighted_taxa = top10,
    threshold_t = threshold_t,
    argentinosaurus_upper_mass_bound_t = tail_band_upper_t
  ),
  fit = fit,
  fit_summary = fit_summary,
  goodness_of_fit = gof,
  sampling_bootstrap = sampling_boot,
  reconstruction_bootstrap = reconstruction_boot,
  combined_bootstrap = combined_boot,
  threshold_stability = stability,
  qq = qq,
  return_fixed = return_fixed,
  return_combined = return_combined,
  empirical_return = empirical_return,
  richness_return_levels = richness_return_levels,
  survival_curve = survival_curve,
  holdout = holdout
)
saveRDS(results, file.path(result_dir, "sauropod_evt_results.rds"))
writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "sessionInfo.txt")
)

cat("\nSAUROPODS\n")
print(audit)
print(fit_summary)
print(richness_return_levels)
invisible(results)
