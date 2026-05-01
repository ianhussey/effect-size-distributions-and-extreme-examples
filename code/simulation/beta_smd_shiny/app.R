# Beta-bound SMD explorer — Shiny app
#
# Mirrors the standalone HTML explorer. Run with:
#   shiny::runApp("app.R")
# from the directory that contains this file.
#
# Requires: shiny, ggplot2

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
})

# --------------------------------------------------------------- core math --

bhatia_davis_max_sd <- function(M, scale_min, scale_max) {
  sqrt((M - scale_min) * (scale_max - M))
}

beta_fit_mom <- function(M, SD, scale_min, scale_max) {
  range <- scale_max - scale_min
  mu <- (M - scale_min) / range
  sigma_tilde <- SD / range
  bd_ceiling <- bhatia_davis_max_sd(M, scale_min, scale_max)
  if (SD > bd_ceiling) {
    return(list(mu = mu, sigma_tilde = sigma_tilde, phi = NA_real_,
                alpha = NA_real_, beta = NA_real_,
                bd_ceiling = bd_ceiling, impossible = TRUE))
  }
  phi <- mu * (1 - mu) / sigma_tilde^2 - 1
  list(mu = mu, sigma_tilde = sigma_tilde, phi = phi,
       alpha = mu * phi, beta = (1 - mu) * phi,
       bd_ceiling = bd_ceiling, impossible = FALSE)
}

predict_sd_at_mean <- function(M_target, phi, scale_min, scale_max) {
  range <- scale_max - scale_min
  mu_t <- (M_target - scale_min) / range
  sqrt(mu_t * (1 - mu_t) / (1 + phi)) * range
}

# --------------------------------------------------------------------- UI --

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, BlinkMacSystemFont, system-ui, sans-serif; }
    h2.section { font-size: 15px; border-bottom: 1px solid #e0e0e0;
                 padding-bottom: 4px; margin-top: 22px; color: #333; }
    .readout { background: #f6f7f9; padding: 12px 14px; border-radius: 6px;
               font-size: 14px; line-height: 1.7; border: 1px solid #ececec;
               font-family: ui-monospace, SFMono-Regular, monospace; }
    .label { color: #666; font-family: -apple-system, system-ui, sans-serif; }
    .value { font-weight: 600; color: #1a3d6b; }
    .err  { color: #b00020; font-weight: 600; padding: 8px 10px;
            background: #fde8eb; border-radius: 4px; margin-bottom: 8px; }
    .warn { color: #8a5a00; padding: 8px 10px; background: #fdf3dd;
            border-radius: 4px; margin-bottom: 8px; font-size: 13px; }
    .lede { font-size: 14px; color: #555; }
    .note { font-size: 13px; color: #555; line-height: 1.55; }
  "))),

  titlePanel("Beta-bound SMD explorer"),
  p(class = "lede",
    "Given an observed mean and SD on a bounded scale, fit a Beta distribution",
    " by method of moments, predict how SD compresses as the mean moves toward",
    " the bounds (constant-precision assumption), and compute SMD bounds for",
    " between-groups Cohen's d or dependent Cohen's d_rm."),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      numericInput("smin", "Scale min", value = 0, step = 1),
      numericInput("smax", "Scale max", value = 100, step = 1),
      uiOutput("m_input"),
      uiOutput("sd_input"),
      hr(),
      selectInput("es_type", "Effect size type",
                  choices = c("Between-groups (Cohen's d)" = "between",
                              "Dependent (Cohen's d_rm)" = "dependent"),
                  selected = "between"),
      conditionalPanel(
        condition = "input.es_type == 'dependent'",
        numericInput("rho", "Pre-post correlation (r)",
                     value = 0.5, min = 0, max = 0.99, step = 0.05)
      )
    ),

    mainPanel(
      width = 9,
      h2("Implied Beta distribution", class = "section"),
      uiOutput("readout"),

      fluidRow(
        column(6,
          h2("Beta PDF (rescaled scale 0–1)", class = "section"),
          plotOutput("pdfPlot", height = "260px")),
        column(6,
          h2("Predicted SD as mean varies", class = "section"),
          plotOutput("sdPlot",  height = "260px"))
      ),

      h2("SMD vs post-intervention mean", class = "section"),
      p(class = "note",
        "Both curves assume the post-intervention distribution has the same",
        " precision phi as the observed one, so its SD is predicted rather",
        " than reported. The naive curve uses the observed SD as denominator;",
        " the Beta-implied curve uses the appropriate two-condition denominator",
        " for the selected effect size type."),
      plotOutput("smdPlot", height = "340px"),

      h2("SMD ceilings at the scale endpoints", class = "section"),
      uiOutput("bounds"),

      h2("Notes on interpretation", class = "section"),
      p(class = "note", strong("Bhatia–Davis ceiling. "),
        "For any distribution on [min, max] with mean M, ",
        "SD <= sqrt((M - min)(max - M)). If observed SD exceeds this, the ",
        "summary statistics are mathematically impossible regardless of ",
        "distributional assumptions."),
      p(class = "note", strong("Constant-phi assumption. "),
        "The Beta-implied SD curve is only as good as the assumption that ",
        "precision doesn't shift between conditions. Treat predicted SDs as ",
        "a reference distribution, not a hard cap."),
      p(class = "note", strong("Effect size conventions. "),
        "Between-groups Cohen's d uses sqrt((SD1^2 + SD2^2)/2) as the ",
        "denominator (equal-n simplification of d_s, equivalent to d_av in ",
        "within-subjects). Dependent Cohen's d_rm = (M2 - M1) * sqrt(2(1-r)) ",
        "/ sqrt(SD1^2 + SD2^2 - 2*r*SD1*SD2), which collapses to the ",
        "between-groups formula at r = 0.")
    )
  )
)

# ----------------------------------------------------------------- server --

server <- function(input, output, session) {

  # Slider inputs that rescale with the current bounds
  output$m_input <- renderUI({
    smin <- input$smin; smax <- input$smax
    if (is.null(smin) || is.null(smax) || smax <= smin) return(NULL)
    rng <- smax - smin
    numericInput("M", "Observed mean (M)",
                 value = isolate(input$M) %||% (smin + 0.75 * rng),
                 step = rng / 200)
  })

  output$sd_input <- renderUI({
    smin <- input$smin; smax <- input$smax
    if (is.null(smin) || is.null(smax) || smax <= smin) return(NULL)
    rng <- smax - smin
    numericInput("SD", "Observed SD",
                 value = isolate(input$SD) %||% (rng * 0.12),
                 step = rng / 500)
  })

  # Validate + fit
  fit_r <- reactive({
    req(input$smin, input$smax, input$M, input$SD)
    smin <- input$smin; smax <- input$smax
    M <- input$M; SD <- input$SD
    if (smax <= smin)        return(list(error = "Scale max must exceed scale min."))
    if (M <= smin || M >= smax) return(list(error = "Mean must lie strictly inside (min, max)."))
    if (SD <= 0)             return(list(error = "SD must be positive."))
    fit <- beta_fit_mom(M, SD, smin, smax)
    list(error = NULL, smin = smin, smax = smax, M = M, SD = SD, fit = fit)
  })

  # Implied-Beta readout
  output$readout <- renderUI({
    r <- fit_r()
    if (!is.null(r$error)) return(div(class = "err", r$error))
    f <- r$fit
    warn_html <- NULL
    if (f$impossible) {
      warn_html <- div(class = "err", sprintf(
        "SD = %.3f exceeds the Bhatia-Davis ceiling of %.3f. The reported summary statistics are mathematically impossible.",
        r$SD, f$bd_ceiling))
    } else if (f$phi < 1) {
      warn_html <- div(class = "warn", sprintf(
        "Implied phi = %.3f < 1: fitted Beta is U-shaped (mass concentrated at the bounds). Plausible for some scales but worth noting.",
        f$phi))
    }
    div(class = "readout",
        warn_html,
        HTML(sprintf(
          paste0("<span class='label'>Rescaled mean mu:</span> <span class='value'>%.4f</span> &nbsp;&middot;&nbsp;",
                 "<span class='label'>Rescaled SD sigma~:</span> <span class='value'>%.4f</span><br>",
                 "<span class='label'>Bhatia-Davis SD ceiling:</span> <span class='value'>%.3f</span> on raw scale",
                 " (sigma~ <= <span class='value'>%.4f</span>)<br>",
                 "<span class='label'>Implied precision phi:</span> <span class='value'>%s</span> &nbsp;&middot;&nbsp;",
                 "<span class='label'>Beta shape:</span> <span class='value'>alpha = %s, beta = %s</span>"),
          f$mu, f$sigma_tilde, f$bd_ceiling, sqrt(f$mu * (1 - f$mu)),
          if (is.na(f$phi)) "&mdash;" else sprintf("%.3f", f$phi),
          if (is.na(f$alpha)) "&mdash;" else sprintf("%.3f", f$alpha),
          if (is.na(f$beta))  "&mdash;" else sprintf("%.3f", f$beta))))
  })

  # Beta PDF plot
  output$pdfPlot <- renderPlot({
    r <- fit_r(); req(is.null(r$error)); f <- r$fit
    if (f$impossible) {
      return(ggplot() + theme_void() +
               annotate("text", x = 0, y = 0, label = "SD exceeds Bhatia-Davis ceiling.\nNo Beta fit possible.",
                        size = 4.5, colour = "#b00020"))
    }
    x <- seq(0.001, 0.999, length.out = 400)
    y <- dbeta(x, shape1 = f$alpha, shape2 = f$beta)
    ggplot(data.frame(x = x, y = y), aes(x, y)) +
      geom_area(fill = "#1a5fb4", alpha = 0.18) +
      geom_line(colour = "#1a5fb4", linewidth = 0.9) +
      labs(x = "rescaled value", y = "density") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
  })

  # SD-vs-mean plot
  output$sdPlot <- renderPlot({
    r <- fit_r(); req(is.null(r$error)); f <- r$fit
    rng <- r$smax - r$smin
    mu_grid <- seq(0.005, 0.995, length.out = 400)
    M_grid  <- r$smin + mu_grid * rng
    bd <- sqrt(mu_grid * (1 - mu_grid)) * rng
    df <- data.frame(M = M_grid, bd = bd)
    p <- ggplot(df, aes(M)) +
      geom_line(aes(y = bd, colour = "Bhatia-Davis ceiling"),
                linetype = "dashed", linewidth = 0.7) +
      labs(x = "mean (raw scale)", y = "predicted SD (raw)", colour = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", legend.text = element_text(size = 9),
            panel.grid.minor = element_blank())
    if (!f$impossible) {
      df$beta_sd <- sqrt(mu_grid * (1 - mu_grid) / (1 + f$phi)) * rng
      p <- p + geom_line(data = df, aes(y = beta_sd, colour = "Beta-implied (constant phi)"),
                         linewidth = 1.1) +
        geom_point(data = data.frame(M = r$M, SD = r$SD),
                   aes(x = M, y = SD), inherit.aes = FALSE,
                   colour = "black", size = 2.5)
    }
    p + scale_colour_manual(values = c("Beta-implied (constant phi)" = "#1a5fb4",
                                       "Bhatia-Davis ceiling" = "#c64600"))
  })

  # SMD plot
  output$smdPlot <- renderPlot({
    r <- fit_r(); req(is.null(r$error)); f <- r$fit
    rng <- r$smax - r$smin
    mu_grid <- seq(1/400, 1 - 1/400, length.out = 400)
    M_grid  <- r$smin + mu_grid * rng
    diff    <- M_grid - r$M
    smd_naive <- diff / r$SD

    es_type  <- input$es_type %||% "between"
    rho      <- if (es_type == "dependent") {
                  min(max(input$rho %||% 0.5, 0), 0.99)
                } else 0
    beta_lab <- if (es_type == "dependent") {
                  sprintf("Beta-implied (Cohen's d_rm, r = %.2f)", rho)
                } else "Beta-implied (between-groups Cohen's d)"

    if (f$impossible) {
      df <- data.frame(M = M_grid, smd = smd_naive,
                       conv = "Naive (observed SD)")
    } else {
      sd_post <- sqrt(mu_grid * (1 - mu_grid) / (1 + f$phi)) * rng
      if (es_type == "dependent") {
        sd_diff  <- sqrt(r$SD^2 + sd_post^2 - 2 * rho * r$SD * sd_post)
        smd_beta <- diff * sqrt(2 * (1 - rho)) / sd_diff
      } else {
        smd_beta <- diff / sqrt((r$SD^2 + sd_post^2) / 2)
      }
      df <- rbind(
        data.frame(M = M_grid, smd = smd_naive, conv = "Naive (observed SD)"),
        data.frame(M = M_grid, smd = smd_beta,  conv = beta_lab)
      )
    }
    df$conv <- factor(df$conv, levels = c("Naive (observed SD)", beta_lab))

    cols <- setNames(c("#444444", "#1a5fb4"),
                     c("Naive (observed SD)", beta_lab))
    ltys <- setNames(c("dashed", "solid"),
                     c("Naive (observed SD)", beta_lab))
    lwds <- setNames(c(0.7, 1.1),
                     c("Naive (observed SD)", beta_lab))

    ggplot(df, aes(M, smd, colour = conv, linetype = conv, linewidth = conv)) +
      geom_hline(yintercept = 0, colour = "grey80") +
      geom_vline(xintercept = r$M, colour = "grey80") +
      geom_line() +
      scale_colour_manual(values = cols) +
      scale_linetype_manual(values = ltys) +
      scale_linewidth_manual(values = lwds) +
      labs(x = "post-intervention mean (raw scale)", y = "SMD",
           colour = NULL, linetype = NULL, linewidth = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", legend.text = element_text(size = 9),
            panel.grid.minor = element_blank())
  })

  # Bounds readout
  output$bounds <- renderUI({
    r <- fit_r()
    if (!is.null(r$error)) return(NULL)
    f <- r$fit
    naive_up   <- (r$smax - r$M) / r$SD
    naive_down <- (r$smin - r$M) / r$SD

    es_type <- input$es_type %||% "between"
    rho     <- if (es_type == "dependent") {
                 min(max(input$rho %||% 0.5, 0), 0.99)
               } else 0
    beta_label <- if (es_type == "dependent") {
                    sprintf("Beta-implied bound (Cohen's d_rm, r = %.2f, mu' approaching bound)", rho)
                  } else "Beta-implied bound (between-groups Cohen's d, mu' approaching bound)"

    beta_up_str   <- "&mdash;"
    beta_down_str <- "&mdash;"
    if (!f$impossible) {
      eps <- 5e-4
      M_up   <- r$smin + (1 - eps) * (r$smax - r$smin)
      M_down <- r$smin + eps * (r$smax - r$smin)
      sd_up   <- predict_sd_at_mean(M_up,   f$phi, r$smin, r$smax)
      sd_down <- predict_sd_at_mean(M_down, f$phi, r$smin, r$smax)
      if (es_type == "dependent") {
        sd_diff_up   <- sqrt(r$SD^2 + sd_up^2   - 2 * rho * r$SD * sd_up)
        sd_diff_down <- sqrt(r$SD^2 + sd_down^2 - 2 * rho * r$SD * sd_down)
        beta_up_str   <- sprintf("%.3f", (M_up   - r$M) * sqrt(2 * (1 - rho)) / sd_diff_up)
        beta_down_str <- sprintf("%.3f", (M_down - r$M) * sqrt(2 * (1 - rho)) / sd_diff_down)
      } else {
        beta_up_str   <- sprintf("%.3f", (M_up   - r$M) / sqrt((r$SD^2 + sd_up^2)   / 2))
        beta_down_str <- sprintf("%.3f", (M_down - r$M) / sqrt((r$SD^2 + sd_down^2) / 2))
      }
    }
    div(class = "readout",
        HTML(sprintf(paste0(
          "<span class='label'>Naive linear bound (observed SD as denominator):</span><br>",
          "&nbsp;&nbsp;Max upward SMD &nbsp;= <span class='value'>%.3f</span> (at M' = %.3f)<br>",
          "&nbsp;&nbsp;Max downward SMD = <span class='value'>%.3f</span> (at M' = %.3f)<br><br>",
          "<span class='label'>%s:</span><br>",
          "&nbsp;&nbsp;Approaching upper bound: SMD &rarr; <span class='value'>%s</span><br>",
          "&nbsp;&nbsp;Approaching lower bound: SMD &rarr; <span class='value'>%s</span>"),
          naive_up, r$smax, naive_down, r$smin,
          beta_label, beta_up_str, beta_down_str)))
  })
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui = ui, server = server)
