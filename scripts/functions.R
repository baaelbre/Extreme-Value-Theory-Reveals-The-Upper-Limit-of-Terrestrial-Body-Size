# some shared functions for the three scripts, to be sourced
# Bastiaan A. Van Velthoven

# Figure style and output

COL <- list(
  blue = "#0072B2",
  orange = "#D55E00",
  yellow = "#E69F00",
  grey = "#3B4045",
  reference = "#6B6B6B",
  blue_band = "#B9DDF0"
)

theme_set(
  theme_classic(base_size = 12, base_family = "sans") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = grid::unit(2.5, "pt"),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      plot.tag = element_text(face = "bold", size = 13),
      plot.title.position = "plot",
      plot.margin = margin(7, 9, 7, 7)
    )
)

check_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

save_figure <- function(plot, filename, directory, width, height, dpi = 400) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  print(plot)

  ggsave(
    file.path(directory, paste0(filename, ".png")),
    plot = plot, width = width, height = height,
    dpi = dpi, bg = "white"
  )
  ggsave(
    file.path(directory, paste0(filename, ".pdf")),
    plot = plot, width = width, height = height,
    bg = "white"
  )

  invisible(plot)
}

write_result <- function(x, directory, filename) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file.path(directory, filename), row.names = FALSE)
}

# Generalized Pareto calculations

q_finite <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probability, type = 8, names = FALSE))
}

 
q_endpoint <- function(x, probability = 0.95) { # infinite when xi>=0
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  sort(x)[max(1, ceiling(probability * length(x)))]
}

# Fit a GPD to exceedances of y above a fixed threshold. In all three analyses,
# y is log(size), so finite endpoints are transformed back afterwards.
fit_gpd_y <- function(y, threshold, min_exceedances = 8) {
  y <- y[is.finite(y)]
  n <- length(y)
  n_exc <- sum(y > threshold)

  if (n_exc < min_exceedances) return(NULL)
  # use coles' package for gpd fits
  raw_fit <-ismev::gpd.fit(y, threshold = threshold, show = FALSE)

  sigma <- unname(raw_fit$mle[1])
  xi <- unname(raw_fit$mle[2])
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(xi)) return(NULL)

  list(
    raw = raw_fit,
    n = n,
    n_exc = n_exc,
    p_u = n_exc / n,
    threshold = threshold,
    sigma = sigma,
    xi = xi,
    endpoint_y = if (xi < 0) threshold - sigma / xi else Inf
  )
}

fit_draw_row <- function(fit, draw) {
  if (is.null(fit)) return(NULL)

  data.frame(
    draw = draw,
    n = fit$n,
    n_exc = fit$n_exc,
    p_u = fit$p_u,
    sigma = fit$sigma,
    xi = fit$xi,
    endpoint_y = fit$endpoint_y
  )
}

fit_from_draw <- function(draw, threshold) {
  list(
    n = draw$n[[1]],
    n_exc = draw$n_exc[[1]],
    p_u = draw$p_u[[1]],
    threshold = threshold,
    sigma = draw$sigma[[1]],
    xi = draw$xi[[1]],
    endpoint_y = draw$endpoint_y[[1]]
  )
}

# Random generation, quantiles, and survival for a GPD excess.
rgpd_excess <- function(n, sigma, xi) {
  p <- runif(n, .Machine$double.eps, 1 - .Machine$double.eps)

  if (abs(xi) < 1e-8) {
    -sigma * log1p(-p)
  } else {
    sigma / xi * ((1 - p)^(-xi) - 1)
  }
}

qgpd_excess <- function(p, sigma, xi) {
  p <- pmin(pmax(as.numeric(p), 0), 1 - .Machine$double.eps)

  if (abs(xi) < 1e-8) {
    -sigma * log1p(-p)
  } else {
    sigma / xi * ((1 - p)^(-xi) - 1)
  }
}

gpd_conditional_survival <- function(excess, sigma, xi) {
  excess <- as.numeric(excess)
  answer <- rep(1, length(excess))
  positive <- excess > 0

  if (!any(positive)) return(answer)

  if (abs(xi) < 1e-8) {
    answer[positive] <- exp(-excess[positive] / sigma)
  } else {
    support <- 1 + xi * excess[positive] / sigma
    answer[positive] <- 0
    valid <- support > 0
    answer[which(positive)[valid]] <- support[valid]^(-1 / xi)
  }

  pmin(pmax(answer, 0), 1)
}

tail_survival_y <- function(y_value, fit) {
  fit$p_u * gpd_conditional_survival(
    y_value - fit$threshold, fit$sigma, fit$xi
  )
}

# The return level is the size with marginal exceedance probability 1 / m.
return_level_y <- function(m, fit) {
  m <- as.numeric(m)
  valid <- m * fit$p_u >= 1
  answer <- rep(NA_real_, length(m))

  if (abs(fit$xi) < 1e-8) {
    answer[valid] <- fit$threshold +
      fit$sigma * log(m[valid] * fit$p_u)
  } else {
    answer[valid] <- fit$threshold +
      fit$sigma / fit$xi * ((m[valid] * fit$p_u)^fit$xi - 1)
  }

  answer
}

gpd_cdf_excess <- function(excess, sigma, xi) {
  pmin(
    pmax(1 - gpd_conditional_survival(excess, sigma, xi), 0),
    1
  )
}

anderson_darling_statistic <- function(excess, sigma, xi) {
  excess <- sort(excess[is.finite(excess) & excess >= 0])
  n <- length(excess)
  if (!n) return(NA_real_)

  eps <- sqrt(.Machine$double.eps)
  u <- sort(pmin(
    pmax(gpd_cdf_excess(excess, sigma, xi), eps),
    1 - eps
  ))
  i <- seq_len(n)

  -n - sum((2 * i - 1) * (log(u) + log1p(-rev(u)))) / n
}


# Simulate excesses from the fitted GPD and refit the model.
parametric_bootstrap <- function(fit, B, include_ad = FALSE) {
  rows <- lapply(seq_len(B), function(b) {
    excess <- rgpd_excess(fit$n_exc, fit$sigma, fit$xi)

    fit_b <- fit_gpd_y(
      fit$threshold + excess,
      fit$threshold,
      min_exceedances = fit$n_exc
    )
    if (is.null(fit_b)) return(NULL)

    fit_b$n <- fit$n
    fit_b$n_exc <- fit$n_exc
    fit_b$p_u <- fit$p_u

    row <- fit_draw_row(fit_b, b)
    if (include_ad) {
      row$anderson_darling <- anderson_darling_statistic(
        excess, fit_b$sigma, fit_b$xi
      )
    }
    row
  })

  dplyr::bind_rows(rows)
}

# Resample allobservations, keep the threshold fixed,
# recalculate threshold membership, and refit 
nonparametric_bootstrap <- function(y, threshold, B, min_exceedances = 8) {
  n <- length(y)

  dplyr::bind_rows(lapply(seq_len(B), function(b) {
    fit_draw_row(
      fit_gpd_y(
        sample(y, n, replace = TRUE),
        threshold,
        min_exceedances
      ),
      b
    )
  }))
}

goodness_of_fit <- function(y, fit, B) {
  excess <- y[y > fit$threshold] - fit$threshold
  observed <- anderson_darling_statistic(excess, fit$sigma, fit$xi)
  draws <- parametric_bootstrap(fit, B, include_ad = TRUE)


  valid <- is.finite(draws$anderson_darling)
  p_value <- (1 + sum(draws$anderson_darling[valid] >= observed)) /
    (1 + sum(valid))

  list(
    summary = data.frame(
      threshold_y = fit$threshold,
      n_exc = fit$n_exc,
      statistic = observed,
      B_requested = B,
      B_successful = sum(valid),
      p_value = p_value,
      reject_at_0_05 = p_value < 0.05
    ),
    draws = draws
  )
}

qq_data <- function(y, fit) {
  empirical <- sort(y[y > fit$threshold] - fit$threshold)
  probability <- (seq_along(empirical) - 0.5) / length(empirical)

  data.frame(
    theoretical = qgpd_excess(probability, fit$sigma, fit$xi),
    empirical = empirical
  )
}

mean_excess_row <- function(y, threshold) {
  excess <- y[y > threshold] - threshold

  data.frame(
    n_exc = length(excess),
    mean_excess = if (length(excess)) mean(excess) else NA_real_,
    se = if (length(excess) > 1) {
      sd(excess) / sqrt(length(excess))
    } else {
      NA_real_
    }
  )
}

return_level_matrix <- function(draws, m, threshold, inverse) {
  t(vapply(seq_len(nrow(draws)), function(i) {
    inverse(return_level_y(m, fit_from_draw(draws[i, ], threshold)))
  }, numeric(length(m))))
}

survival_matrix <- function(draws, y_grid, threshold) {
  t(vapply(seq_len(nrow(draws)), function(i) {
    tail_survival_y(y_grid, fit_from_draw(draws[i, ], threshold))
  }, numeric(length(y_grid))))
}

summarise_curve <- function(matrix) {
  data.frame(
    lower95 = apply(matrix, 2, q_finite, probability = 0.025),
    median = apply(matrix, 2, q_finite, probability = 0.50),
    upper95 = apply(matrix, 2, q_finite, probability = 0.975)
  )
}
