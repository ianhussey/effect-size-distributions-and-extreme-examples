# Beta-bound SMD analysis
#
# Given a mean M and SD on a bounded scale [scale_min, scale_max], fit a Beta
# distribution by method of moments and use the constant-precision assumption
# to derive bounds on plausible SMDs.
#
# Core functions:
#   bhatia_davis_max_sd() — parameter-free SD ceiling
#   beta_fit_mom()        — method-of-moments Beta fit
#   predict_sd_at_mean()  — SD predicted at a new mean (constant phi)
#   smd_bounds()          — naive vs Beta-implied SMD ceilings
#   smd_curves()          — full SMD-vs-target-mean curves under three conventions
#
# Usage examples are at the bottom.

# ---------------------------------------------------------------- core math --

#' Bhatia-Davis upper bound on SD given a mean and bounds.
#' For ANY distribution on [scale_min, scale_max] with the given mean,
#' SD <= sqrt((M - scale_min) * (scale_max - M)).
bhatia_davis_max_sd <- function(M, scale_min, scale_max) {
  stopifnot(scale_max > scale_min, M >= scale_min, M <= scale_max)
  sqrt((M - scale_min) * (scale_max - M))
}

#' Method-of-moments Beta fit on a bounded scale.
#' Returns rescaled mean (mu), rescaled SD (sigma_tilde), precision (phi),
#' and shape parameters (alpha, beta). NA values + warning if SD exceeds
#' the Bhatia-Davis ceiling (i.e. summary stats are impossible).
beta_fit_mom <- function(M, SD, scale_min = 0, scale_max = 1) {
  stopifnot(scale_max > scale_min, SD > 0)
  if (M <= scale_min || M >= scale_max) {
    stop("Mean must lie strictly inside (scale_min, scale_max).")
  }
  range <- scale_max - scale_min
  mu <- (M - scale_min) / range
  sigma_tilde <- SD / range
  bd_ceiling <- bhatia_davis_max_sd(M, scale_min, scale_max)
  if (SD > bd_ceiling) {
    warning(sprintf(
      "SD = %.4f exceeds the Bhatia-Davis ceiling %.4f for mean %.4f on [%.4f, %.4f]. Summary statistics are impossible.",
      SD, bd_ceiling, M, scale_min, scale_max
    ))
    return(list(mu = mu, sigma_tilde = sigma_tilde, phi = NA_real_,
                alpha = NA_real_, beta = NA_real_,
                bd_ceiling = bd_ceiling, impossible = TRUE))
  }
  phi <- mu * (1 - mu) / sigma_tilde^2 - 1
  list(mu = mu, sigma_tilde = sigma_tilde, phi = phi,
       alpha = mu * phi, beta = (1 - mu) * phi,
       bd_ceiling = bd_ceiling, impossible = FALSE)
}

#' Predicted SD at a candidate mean, holding precision phi constant.
#' Returns SD on the raw scale.
predict_sd_at_mean <- function(M_target, phi, scale_min = 0, scale_max = 1) {
  range <- scale_max - scale_min
  mu_t <- (M_target - scale_min) / range
  sqrt(mu_t * (1 - mu_t) / (1 + phi)) * range
}

# ---------------------------------------------------------------- SMD bounds --

#' Limiting SMDs at the scale endpoints under the three denominator conventions.
#' "Naive" uses the observed SD. "Pooled" and "post" use the Beta-implied
#' post-intervention SD via constant phi.
smd_bounds <- function(M, SD, scale_min = 0, scale_max = 1, eps = 5e-4) {
  fit <- beta_fit_mom(M, SD, scale_min, scale_max)
  range <- scale_max - scale_min

  naive_up   <- (scale_max - M) / SD
  naive_down <- (scale_min - M) / SD

  if (fit$impossible) {
    return(list(fit = fit,
                naive   = c(up = naive_up, down = naive_down),
                pooled  = c(up = NA_real_, down = NA_real_),
                post    = c(up = NA_real_, down = NA_real_)))
  }

  # Approach the bounds with a tiny epsilon
  M_up   <- scale_min + (1 - eps) * range
  M_down <- scale_min + eps * range
  sd_up   <- predict_sd_at_mean(M_up,   fit$phi, scale_min, scale_max)
  sd_down <- predict_sd_at_mean(M_down, fit$phi, scale_min, scale_max)

  pooled_up   <- (M_up   - M) / sqrt((SD^2 + sd_up^2)   / 2)
  pooled_down <- (M_down - M) / sqrt((SD^2 + sd_down^2) / 2)
  post_up     <- (M_up   - M) / sd_up
  post_down   <- (M_down - M) / sd_down

  list(fit = fit,
       naive  = c(up = naive_up,  down = naive_down),
       pooled = c(up = pooled_up, down = pooled_down),
       post   = c(up = post_up,   down = post_down))
}

#' Full SMD-vs-target-mean curves under three denominator conventions.
#' Returns a data.frame for plotting.
smd_curves <- function(M, SD, scale_min = 0, scale_max = 1, n = 200) {
  fit <- beta_fit_mom(M, SD, scale_min, scale_max)
  range <- scale_max - scale_min
  mu_grid <- seq(1 / n, 1 - 1 / n, length.out = n)
  M_grid  <- scale_min + mu_grid * range
  sd_post <- if (fit$impossible) rep(NA_real_, n) else
    sqrt(mu_grid * (1 - mu_grid) / (1 + fit$phi)) * range
  diff    <- M_grid - M
  data.frame(
    target_mean = M_grid,
    sd_post     = sd_post,
    smd_naive   = diff / SD,
    smd_post    = ifelse(is.na(sd_post) | sd_post < 1e-9, NA_real_, diff / sd_post),
    smd_pooled  = ifelse(is.na(sd_post), NA_real_,
                         diff / sqrt((SD^2 + sd_post^2) / 2))
  )
}

# ----------------------------------------------------------------- demo plots --

#' Plot the SD-vs-mean curve and the three SMD curves for a given (M, SD, bounds).
plot_beta_smd <- function(M, SD, scale_min = 0, scale_max = 1) {
  fit <- beta_fit_mom(M, SD, scale_min, scale_max)
  curves <- smd_curves(M, SD, scale_min, scale_max, n = 400)
  range  <- scale_max - scale_min
  mu_grid <- seq(0.005, 0.995, length.out = 400)
  M_grid  <- scale_min + mu_grid * range
  sd_curve <- if (fit$impossible) rep(NA_real_, length(mu_grid)) else
    sqrt(mu_grid * (1 - mu_grid) / (1 + fit$phi)) * range
  bd_curve <- sqrt(mu_grid * (1 - mu_grid)) * range

  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), cex = 0.9)
  on.exit(par(op))

  # SD vs mean
  plot(M_grid, bd_curve, type = "l", lty = 2, col = "firebrick",
       xlab = "mean (raw scale)", ylab = "predicted SD (raw)",
       main = sprintf("SD vs mean  (phi = %.2f)",
                      ifelse(is.na(fit$phi), NA, fit$phi)),
       ylim = c(0, max(bd_curve)))
  if (!fit$impossible) lines(M_grid, sd_curve, col = "navy", lwd = 2)
  points(M, SD, pch = 19)
  legend("topright", bty = "n", cex = 0.85,
         legend = c("Beta-implied (constant phi)", "Bhatia-Davis ceiling", "observed (M, SD)"),
         col = c("navy", "firebrick", "black"),
         lty = c(1, 2, NA), pch = c(NA, NA, 19), lwd = c(2, 1, NA))

  # SMD curves
  ylim <- range(c(curves$smd_naive, curves$smd_pooled,
                  pmin(pmax(curves$smd_post, -10), 10)), na.rm = TRUE)
  plot(curves$target_mean, curves$smd_naive, type = "l", lty = 2, col = "grey30",
       xlab = "target (post) mean", ylab = "SMD", ylim = ylim,
       main = "SMD vs target mean")
  if (!fit$impossible) {
    lines(curves$target_mean, curves$smd_pooled, col = "navy", lwd = 2.5)
    lines(curves$target_mean, pmin(pmax(curves$smd_post, ylim[1]), ylim[2]),
          col = "firebrick", lwd = 2)
  }
  abline(h = 0, col = "grey80")
  abline(v = M, col = "grey80")
  legend("topleft", bty = "n", cex = 0.8,
         legend = c("Naive (control SD)", "Pooled SD (Beta)", "Post-only SD (Beta)"),
         col = c("grey30", "navy", "firebrick"),
         lty = c(2, 1, 1), lwd = c(1.5, 2.5, 2))

  invisible(list(fit = fit, curves = curves))
}

# --------------------------------------------------------------------- examples

if (interactive()) {
  # Worked example matching the back-of-envelope case
  res <- smd_bounds(M = 75, SD = 12, scale_min = 0, scale_max = 100)
  print(res$fit)
  print(rbind(naive = res$naive, pooled = res$pooled, post = res$post))

  plot_beta_smd(M = 75, SD = 12, scale_min = 0, scale_max = 100)

  # Bhatia-Davis violation: try M = 98, SD = 20 on [0, 100]
  try(smd_bounds(M = 98, SD = 20, scale_min = 0, scale_max = 100))

  # Vectorised application: flag suspect SMDs in a meta-analytic dataset
  # df <- data.frame(M_pre = ..., SD_pre = ..., M_post = ..., scale_min = ..., scale_max = ...,
  #                  reported_SMD = ...)
  # df$beta_pooled_ceiling <- mapply(function(Mpre, SDpre, smin, smax) {
  #   smd_bounds(Mpre, SDpre, smin, smax)$pooled["up"]
  # }, df$M_pre, df$SD_pre, df$scale_min, df$scale_max)
  # df$suspect <- df$reported_SMD > df$beta_pooled_ceiling
}
