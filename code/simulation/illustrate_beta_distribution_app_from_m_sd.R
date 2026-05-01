library(shiny)
library(ggplot2)
library(scales)

# ── helpers ──────────────────────────────────────────────────────────────────

# Method-of-moments: given mean & variance on [0,1], return alpha & beta
mom_beta <- function(mu, sigma2) {
  if (sigma2 <= 0 || sigma2 >= mu * (1 - mu)) return(NULL)
  kappa <- (mu * (1 - mu) / sigma2) - 1
  list(alpha = mu * kappa, beta = (1 - mu) * kappa)
}

shape_label_fn <- function(al, be) {
  if (al == 1 && be == 1) "Uniform"
  else if (al < 1 && be < 1) "U-shaped"
  else if (al < 1) "J-shaped (left)"
  else if (be < 1) "J-shaped (right)"
  else if (al == be) "Symmetric"
  else if (al > be) "Left-skewed"
  else "Right-skewed"
}

# Compute pointwise 95% CI band on the density curve via parametric bootstrap.
# Samples N obs from Beta(al,be), re-fits alpha/beta via MOM, recomputes density.
# Returns data frame with x, y_lo, y_hi in the display scale [lo, hi].
beta_density_ci <- function(al, be, lo, hi, n, n_boot = 1000) {
  x_unit <- seq(0.001, 0.999, length.out = 300)
  range  <- hi - lo

  mat <- matrix(NA_real_, nrow = n_boot, ncol = length(x_unit))
  mu0    <- al / (al + be)
  sigma2_0 <- (al * be) / ((al + be)^2 * (al + be + 1))

  for (i in seq_len(n_boot)) {
    # Draw N samples, re-estimate alpha/beta via MOM
    samp   <- rbeta(n, al, be)
    mu_s   <- mean(samp)
    s2_s   <- var(samp) * (n - 1) / n   # MLE variance
    # Guard against degenerate samples
    if (s2_s <= 0 || s2_s >= mu_s * (1 - mu_s)) next
    kappa  <- (mu_s * (1 - mu_s) / s2_s) - 1
    al_s   <- mu_s * kappa
    be_s   <- (1 - mu_s) * kappa
    if (al_s <= 0 || be_s <= 0) next
    mat[i, ] <- dbeta(x_unit, al_s, be_s) / range
  }

  y_lo <- apply(mat, 2, quantile, 0.025, na.rm = TRUE)
  y_hi <- apply(mat, 2, quantile, 0.975, na.rm = TRUE)

  data.frame(x = lo + x_unit * range, y_lo = y_lo, y_hi = y_hi)
}

beta_plot_fn <- function(al, be, lo, hi, ci_band = NULL) {
  x_unit <- seq(0.001, 0.999, length.out = 500)
  y      <- dbeta(x_unit, shape1 = al, shape2 = be)
  range  <- hi - lo
  x_plot <- lo + x_unit * range
  y_plot <- y / range
  mean_x <- lo + (al / (al + be)) * range
  mean_nudge <- range * 0.03

  df <- data.frame(x = x_plot, y = y_plot)

  p <- ggplot(df, aes(x = x, y = y))

  # CI ribbon behind the main curve
  if (!is.null(ci_band)) {
    p <- p +
      geom_ribbon(data = ci_band,
                  aes(x = x, ymin = y_lo, ymax = y_hi),
                  inherit.aes = FALSE,
                  fill = "#4dd8a0", alpha = 0.20) +
      geom_line(data = ci_band, aes(x = x, y = y_lo),
                inherit.aes = FALSE,
                colour = "#4dd8a0", linewidth = 0.6, linetype = "dashed", alpha = 0.8) +
      geom_line(data = ci_band, aes(x = x, y = y_hi),
                inherit.aes = FALSE,
                colour = "#4dd8a0", linewidth = 0.6, linetype = "dashed", alpha = 0.8)
  }

  p +
    geom_area(fill = "#7c6af7", alpha = 0.25) +
    geom_line(colour = "#7c6af7", linewidth = 1.4) +
    geom_vline(xintercept = mean_x, linetype = "dashed",
               colour = "#f0a868", linewidth = 0.8, alpha = 0.8) +
    annotate("text", x = mean_x + mean_nudge, y = max(y_plot) * 0.92,
             label = "mean", colour = "#f0a868", size = 5,
             family = "sans", hjust = 0) +
    scale_x_continuous(limits = c(lo, hi), expand = c(0.01, 0),
                       breaks = scales::breaks_pretty(n = 8)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    labs(x = "Score on native scale", y = "Density", title = NULL) +
    theme_minimal(base_size = 16) +
    theme(
      plot.background  = element_rect(fill = "#16151f", colour = NA),
      panel.background = element_rect(fill = "#16151f", colour = NA),
      panel.grid.major = element_line(colour = "#2a2838", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(colour = "#8a86a8", family = "sans", size = 15),
      axis.title       = element_text(colour = "#c8c4e0", family = "sans", size = 15),
      axis.ticks       = element_line(colour = "#2a2838"),
      plot.margin      = margin(12, 20, 12, 12)
    )
}

# ── CSS ───────────────────────────────────────────────────────────────────────

app_css <- "
  @import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500&display=swap');

  * { box-sizing: border-box; margin: 0; padding: 0; }

  body {
    background-color: #0e0e14;
    color: #e8e6f0;
    font-family: 'DM Sans', sans-serif;
    min-height: 100vh;
  }

  /* ── tabs ── */
  .nav-tabs {
    border-bottom: 1px solid #2a2838 !important;
    margin-bottom: 32px !important;
  }

  .nav-tabs > li > a {
    font-family: 'Space Mono', monospace !important;
    font-size: 11px !important;
    letter-spacing: 2px !important;
    text-transform: uppercase !important;
    color: #6b668a !important;
    background: transparent !important;
    border: none !important;
    border-bottom: 2px solid transparent !important;
    border-radius: 0 !important;
    padding: 10px 20px !important;
    transition: color 0.2s !important;
  }

  .nav-tabs > li > a:hover {
    color: #c8c4e0 !important;
    background: transparent !important;
    border-color: transparent !important;
  }

  .nav-tabs > li.active > a,
  .nav-tabs > li.active > a:focus,
  .nav-tabs > li.active > a:hover {
    color: #7c6af7 !important;
    background: transparent !important;
    border: none !important;
    border-bottom: 2px solid #7c6af7 !important;
  }

  .tab-content { background: transparent !important; }

  /* ── layout ── */
  .main-container {
    max-width: 1100px;
    margin: 0 auto;
    padding: 40px 24px;
  }

  .header {
    margin-bottom: 40px;
    border-left: 3px solid #7c6af7;
    padding-left: 16px;
  }

  .header h1 {
    font-family: 'Space Mono', monospace;
    font-size: 28px;
    font-weight: 700;
    color: #ffffff;
    letter-spacing: -0.5px;
    line-height: 1.2;
  }

  .header p {
    margin-top: 8px;
    color: #8a86a8;
    font-size: 14px;
    font-weight: 300;
    letter-spacing: 0.3px;
  }

  .layout {
    display: grid;
    grid-template-columns: 280px 1fr;
    gap: 24px;
    align-items: start;
  }

  /* ── panels ── */
  .controls-panel {
    background: #16151f;
    border: 1px solid #2a2838;
    border-radius: 12px;
    padding: 28px 24px;
  }

  .controls-panel h3 {
    font-family: 'Space Mono', monospace;
    font-size: 11px;
    letter-spacing: 2px;
    text-transform: uppercase;
    color: #7c6af7;
    margin-bottom: 28px;
  }

  .param-group { margin-bottom: 28px; }

  .param-label {
    font-family: 'Space Mono', monospace;
    font-size: 13px;
    color: #c8c4e0;
    margin-bottom: 10px;
    display: flex;
    justify-content: space-between;
    align-items: center;
  }

  .param-label .value-badge {
    background: #7c6af7;
    color: white;
    padding: 2px 8px;
    border-radius: 20px;
    font-size: 12px;
    font-weight: 700;
    min-width: 42px;
    text-align: center;
  }

  /* sliders */
  .irs--shiny .irs-bar {
    background: #7c6af7 !important;
    border-top: 1px solid #7c6af7 !important;
    border-bottom: 1px solid #7c6af7 !important;
  }
  .irs--shiny .irs-handle { background: #7c6af7 !important; border: 2px solid #fff !important; }
  .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single { background: #7c6af7 !important; }
  .irs--shiny .irs-line { background: #2a2838 !important; border: none !important; }

  /* stats grid */
  .stats-grid {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 10px;
    margin-top: 28px;
    padding-top: 24px;
    border-top: 1px solid #2a2838;
  }

  .stat-box {
    background: #0e0e14;
    border: 1px solid #2a2838;
    border-radius: 8px;
    padding: 12px;
    text-align: center;
  }

  .stat-box .stat-label {
    font-size: 10px;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    color: #6b668a;
    margin-bottom: 6px;
  }

  .stat-box .stat-value {
    font-family: 'Space Mono', monospace;
    font-size: 16px;
    font-weight: 700;
    color: #e8e6f0;
  }

  .stat-box .stat-value.derived {
    color: #7c6af7;
  }

  /* plot panel */
  .plot-panel {
    background: #16151f;
    border: 1px solid #2a2838;
    border-radius: 12px;
    padding: 24px;
  }

  .plot-panel .plot-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 20px;
  }

  .shape-tag {
    font-family: 'Space Mono', monospace;
    font-size: 11px;
    letter-spacing: 1px;
    color: #7c6af7;
    background: rgba(124, 106, 247, 0.1);
    border: 1px solid rgba(124, 106, 247, 0.3);
    border-radius: 20px;
    padding: 4px 12px;
  }

  .plot-title {
    font-family: 'Space Mono', monospace;
    font-size: 11px;
    letter-spacing: 2px;
    text-transform: uppercase;
    color: #6b668a;
  }

  /* native scale (tab 1) */
  .native-section {
    margin-top: 24px;
    padding-top: 24px;
    border-top: 1px solid #2a2838;
  }

  .native-toggle {
    display: flex;
    align-items: center;
    gap: 10px;
    cursor: pointer;
    margin-bottom: 0;
  }

  .native-toggle input[type='checkbox'] {
    appearance: none;
    -webkit-appearance: none;
    width: 36px;
    height: 20px;
    background: #2a2838;
    border-radius: 10px;
    position: relative;
    cursor: pointer;
    transition: background 0.2s;
    flex-shrink: 0;
  }

  .native-toggle input[type='checkbox']:checked { background: #7c6af7; }

  .native-toggle input[type='checkbox']::after {
    content: '';
    position: absolute;
    width: 14px;
    height: 14px;
    background: white;
    border-radius: 50%;
    top: 3px;
    left: 3px;
    transition: left 0.2s;
  }

  .native-toggle input[type='checkbox']:checked::after { left: 19px; }

  .native-toggle-label {
    font-family: 'Space Mono', monospace;
    font-size: 12px;
    color: #c8c4e0;
    letter-spacing: 0.3px;
  }

  .native-inputs {
    margin-top: 18px;
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 10px;
  }

  .native-input-group label {
    display: block;
    font-family: 'Space Mono', monospace;
    font-size: 10px;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    color: #6b668a;
    margin-bottom: 6px;
  }

  .native-input-group input[type='number'],
  .inv-input-group input[type='number'] {
    width: 100%;
    background: #0e0e14;
    border: 1px solid #2a2838;
    border-radius: 6px;
    color: #e8e6f0;
    font-family: 'Space Mono', monospace;
    font-size: 14px;
    padding: 8px 10px;
    outline: none;
    transition: border-color 0.2s;
  }

  .native-input-group input[type='number']:focus,
  .inv-input-group input[type='number']:focus { border-color: #7c6af7; }

  .native-input-group input[type='number']::-webkit-inner-spin-button,
  .native-input-group input[type='number']::-webkit-outer-spin-button,
  .inv-input-group input[type='number']::-webkit-inner-spin-button,
  .inv-input-group input[type='number']::-webkit-outer-spin-button { -webkit-appearance: none; }

  /* ── inverse tab inputs ── */
  .inv-grid {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 14px;
    margin-bottom: 0;
  }

  .inv-input-group label {
    display: block;
    font-family: 'Space Mono', monospace;
    font-size: 10px;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    color: #6b668a;
    margin-bottom: 6px;
  }

  .error-msg {
    margin-top: 16px;
    padding: 10px 14px;
    background: rgba(240, 100, 100, 0.1);
    border: 1px solid rgba(240, 100, 100, 0.3);
    border-radius: 8px;
    font-family: 'Space Mono', monospace;
    font-size: 11px;
    color: #f07070;
    line-height: 1.5;
    display: none;
  }

  /* Responsive */
  @media (max-width: 768px) {
    .layout { grid-template-columns: 1fr; }
  }
"

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  tags$head(tags$style(HTML(app_css))),

  div(class = "main-container",

    div(class = "header",
      tags$h1("Beta Distribution"),
      tags$p("Interactive visualiser — explore by parameters or recover parameters from summary statistics")
    ),

    tabsetPanel(id = "main_tabs",

      # ── Tab 1 : forward ────────────────────────────────────────────────────
      tabPanel("α / β → Distribution",

        br(),
        div(class = "layout",

          div(class = "controls-panel",
            tags$h3("Parameters"),

            div(class = "param-group",
              div(class = "param-label",
                span("α (alpha)"),
                span(class = "value-badge", textOutput("alpha_val", inline = TRUE))
              ),
              sliderInput("alpha", label = NULL, min = 0.01, max = 20, value = 8, step = 0.1, width = "100%")
            ),

            div(class = "param-group",
              div(class = "param-label",
                span("β (beta)"),
                span(class = "value-badge", textOutput("beta_val", inline = TRUE))
              ),
              sliderInput("beta_param", label = NULL, min = 0.01, max = 20, value = 8, step = 0.1, width = "100%")
            ),

            div(class = "native-section",
              tags$h3("Score on native scale"),
              div(class = "native-toggle",
                tags$input(type = "checkbox", id = "use_native",
                  onclick = "Shiny.setInputValue('use_native', this.checked)"),
                span(class = "native-toggle-label", "Rescale to native units")
              ),
              conditionalPanel(
                condition = "input.use_native == true",
                div(class = "native-inputs",
                  div(class = "native-input-group",
                    tags$label("Min"),
                    tags$input(type = "number", id = "native_min", value = "0",
                      oninput = "Shiny.setInputValue('native_min', parseFloat(this.value))")
                  ),
                  div(class = "native-input-group",
                    tags$label("Max"),
                    tags$input(type = "number", id = "native_max", value = "63",
                      oninput = "Shiny.setInputValue('native_max', parseFloat(this.value))")
                  )
                )
              )
            ),

            div(class = "native-section",
              tags$h3("Sample size"),
              div(class = "native-toggle",
                tags$input(type = "checkbox", id = "use_n_t1",
                  onclick = "Shiny.setInputValue('use_n_t1', this.checked)"),
                span(class = "native-toggle-label", "Show 95% CI on mean")
              ),
              conditionalPanel(
                condition = "input.use_n_t1 == true",
                div(class = "native-inputs",
                  div(class = "native-input-group",
                    tags$label("N"),
                    tags$input(type = "number", id = "n_t1", value = "100", min = "2",
                      oninput = "Shiny.setInputValue('n_t1', parseFloat(this.value))")
                  )
                )
              )
            ),

            div(class = "stats-grid",
              div(class = "stat-box",
                div(class = "stat-label", "Mean"),
                div(class = "stat-value", textOutput("mean_val", inline = TRUE))
              ),
              div(class = "stat-box",
                div(class = "stat-label", "SD"),
                div(class = "stat-value", textOutput("var_val", inline = TRUE))
              ),
              div(class = "stat-box",
                div(class = "stat-label", "Skewness"),
                div(class = "stat-value", textOutput("skew_val", inline = TRUE))
              )
            )
          ),

          div(class = "plot-panel",
            div(class = "plot-header",
              span(class = "plot-title", "Probability Density"),
              span(class = "shape-tag", textOutput("shape_label", inline = TRUE))
            ),
            plotOutput("beta_plot", height = "420px")
          )
        )
      ), # end tab 1

      # ── Tab 2 : inverse ────────────────────────────────────────────────────
      tabPanel("Summary stats → α / β",

        br(),
        div(class = "layout",

          div(class = "controls-panel",
            tags$h3("Summary Statistics"),

            div(class = "inv-grid",
              div(class = "inv-input-group",
                tags$label("Mean"),
                tags$input(type = "number", id = "inv_mean", value = "20",
                  oninput = "Shiny.setInputValue('inv_mean', parseFloat(this.value))")
              ),
              div(class = "inv-input-group",
                tags$label("SD"),
                tags$input(type = "number", id = "inv_sd", value = "10", min = "0.001",
                  oninput = "Shiny.setInputValue('inv_sd', parseFloat(this.value))")
              ),
              div(class = "inv-input-group",
                tags$label("Scale Min"),
                tags$input(type = "number", id = "inv_min", value = "0",
                  oninput = "Shiny.setInputValue('inv_min', parseFloat(this.value))")
              ),
              div(class = "inv-input-group",
                tags$label("Scale Max"),
                tags$input(type = "number", id = "inv_max", value = "63",
                  oninput = "Shiny.setInputValue('inv_max', parseFloat(this.value))")
              )
            ),

            # POMP toggle
            div(class = "native-section",
              tags$h3("Display scale"),
              div(class = "native-toggle",
                tags$input(type = "checkbox", id = "inv_use_pomp",
                  onclick = "Shiny.setInputValue('inv_use_pomp', this.checked)"),
                span(class = "native-toggle-label", "Rescale to standardized units (POMP)")
              )
            ),

            # N / CI toggle
            div(class = "native-section",
              tags$h3("Sample size"),
              div(class = "native-toggle",
                tags$input(type = "checkbox", id = "use_n_t2",
                  onclick = "Shiny.setInputValue('use_n_t2', this.checked)"),
                span(class = "native-toggle-label", "Show 95% CI on mean")
              ),
              conditionalPanel(
                condition = "input.use_n_t2 == true",
                div(class = "native-inputs",
                  div(class = "native-input-group",
                    tags$label("N"),
                    tags$input(type = "number", id = "n_t2", value = "100", min = "2",
                      oninput = "Shiny.setInputValue('n_t2', parseFloat(this.value))")
                  )
                )
              )
            ),

            # Derived parameters
            div(class = "stats-grid",
              div(class = "stat-box",
                div(class = "stat-label", "α (alpha)"),
                div(class = "stat-value derived", textOutput("inv_alpha", inline = TRUE))
              ),
              div(class = "stat-box",
                div(class = "stat-label", "β (beta)"),
                div(class = "stat-value derived", textOutput("inv_beta_out", inline = TRUE))
              ),
              div(class = "stat-box",
                div(class = "stat-label", "Skewness"),
                div(class = "stat-value", textOutput("inv_skew", inline = TRUE))
              )
            ),

            conditionalPanel(
              condition = "input.inv_use_pomp == true",
              div(class = "stats-grid",
                div(class = "stat-box",
                  div(class = "stat-label", "POMP Mean"),
                  div(class = "stat-value", textOutput("inv_pomp_mean", inline = TRUE))
                ),
                div(class = "stat-box",
                  div(class = "stat-label", "POMP SD"),
                  div(class = "stat-value", textOutput("inv_pomp_sd", inline = TRUE))
                )
              )
            ),

            uiOutput("inv_error")
          ),

          div(class = "plot-panel",
            div(class = "plot-header",
              span(class = "plot-title", "Probability Density"),
              span(class = "shape-tag", textOutput("inv_shape_tag", inline = TRUE))
            ),
            plotOutput("inv_plot", height = "420px")
          )
        )
      ) # end tab 2

    ) # end tabsetPanel
  ) # end main-container
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Tab 1 reactives ──────────────────────────────────────────────────────

  a <- reactive(input$alpha)
  b <- reactive(input$beta_param)

  use_native <- reactive(isTRUE(input$use_native))
  use_n_t1   <- reactive(isTRUE(input$use_n_t1))

  native_min <- reactive({ v <- input$native_min; if (is.null(v) || is.na(v)) 0 else v })
  native_max <- reactive({ v <- input$native_max; if (is.null(v) || is.na(v)) 63 else v })
  n_t1       <- reactive({ v <- input$n_t1; if (is.null(v) || is.na(v) || v < 2) 100 else v })

  to_native <- function(p) native_min() + p * (native_max() - native_min())

  beta_mean <- reactive(a() / (a() + b()))

  # CI helper for tab 1: returns bounds in current display scale
  t1_ci_band <- reactive({
    if (!use_n_t1()) return(NULL)
    lo <- if (use_native()) native_min() else 0
    hi <- if (use_native()) native_max() else 1
    beta_density_ci(a(), b(), lo, hi, n_t1())
  })

  output$alpha_val <- renderText(a())
  output$beta_val  <- renderText(b())

  output$mean_val <- renderText({
    m <- beta_mean()
    if (use_native()) round(to_native(m), 2) else round(m, 3)
  })

  output$var_val <- renderText({
    v <- (a() * b()) / ((a() + b())^2 * (a() + b() + 1))
    if (use_native()) round(sqrt(v) * (native_max() - native_min()), 2) else round(sqrt(v), 4)
  })

  output$skew_val <- renderText({
    s <- (2 * (b() - a()) * sqrt(a() + b() + 1)) / ((a() + b() + 2) * sqrt(a() * b()))
    round(s, 3)
  })

  output$shape_label <- renderText(shape_label_fn(a(), b()))

  output$beta_plot <- renderPlot({
    ci_band <- t1_ci_band()
    if (use_native()) {
      lo <- native_min(); hi <- native_max()
      if (is.na(lo) || is.na(hi) || hi <= lo) return(NULL)
      beta_plot_fn(a(), b(), lo, hi, ci_band)
    } else {
      beta_plot_fn(a(), b(), 0, 1, ci_band) +
        labs(x = "Score on proportion-of-maximum-possible scale")
    }
  }, bg = "#16151f")

  # ── Tab 2 reactives ──────────────────────────────────────────────────────

  inv_mean <- reactive({ v <- input$inv_mean; if (is.null(v) || is.na(v)) 20 else v })
  inv_sd   <- reactive({ v <- input$inv_sd;   if (is.null(v) || is.na(v)) 10 else v })
  inv_min  <- reactive({ v <- input$inv_min;  if (is.null(v) || is.na(v)) 0  else v })
  inv_max  <- reactive({ v <- input$inv_max;  if (is.null(v) || is.na(v)) 63 else v })

  inv_use_pomp <- reactive(isTRUE(input$inv_use_pomp))
  use_n_t2     <- reactive(isTRUE(input$use_n_t2))
  n_t2         <- reactive({ v <- input$n_t2; if (is.null(v) || is.na(v) || v < 2) 100 else v })

  inv_params <- reactive({
    lo <- inv_min(); hi <- inv_max(); mn <- inv_mean(); sd <- inv_sd()
    range <- hi - lo
    if (is.na(range) || range <= 0) return(list(error = "Scale max must be greater than scale min."))
    if (mn <= lo || mn >= hi)       return(list(error = "Mean must be strictly between scale min and max."))
    if (sd <= 0)                    return(list(error = "SD must be greater than 0."))

    # Convert to [0,1]
    mu     <- (mn - lo) / range
    sigma2 <- (sd / range)^2

    if (sigma2 >= mu * (1 - mu))
      return(list(error = sprintf(
        "SD is too large relative to the mean and scale range.\nMaximum possible SD at this mean: %.3f",
        sqrt(mu * (1 - mu)) * range
      )))

    p <- mom_beta(mu, sigma2)
    if (is.null(p)) return(list(error = "Could not compute parameters. Check your inputs."))
    list(alpha = p$alpha, beta = p$beta)
  })

  output$inv_error <- renderUI({
    p <- inv_params()
    if (!is.null(p$error)) {
      div(class = "error-msg", style = "display:block;", p$error)
    }
  })

  output$inv_alpha <- renderText({
    p <- inv_params()
    if (!is.null(p$error)) "—" else round(p$alpha, 3)
  })

  output$inv_beta_out <- renderText({
    p <- inv_params()
    if (!is.null(p$error)) "—" else round(p$beta, 3)
  })

  output$inv_skew <- renderText({
    p <- inv_params()
    if (!is.null(p$error)) return("—")
    al <- p$alpha; be <- p$beta
    s <- (2 * (be - al) * sqrt(al + be + 1)) / ((al + be + 2) * sqrt(al * be))
    round(s, 3)
  })

  output$inv_shape <- renderText({
    p <- inv_params()
    if (!is.null(p$error)) return("—")
    shape_label_fn(p$alpha, p$beta)
  })

  output$inv_shape_tag <- renderText({
    p <- inv_params()
    if (!is.null(p$error)) return("Invalid inputs")
    shape_label_fn(p$alpha, p$beta)
  })

  output$inv_pomp_mean <- renderText({
    p <- inv_params()
    if (!is.null(p$error)) return("—")
    lo <- inv_min(); hi <- inv_max()
    range <- hi - lo
    mu <- (inv_mean() - lo) / range
    round(mu * 100, 2)
  })

  output$inv_pomp_sd <- renderText({
    p <- inv_params()
    if (!is.null(p$error)) return("—")
    lo <- inv_min(); hi <- inv_max()
    range <- hi - lo
    round((inv_sd() / range) * 100, 2)
  })

  # CI band for tab 2 (in display scale: native or POMP)
  t2_ci_band <- reactive({
    if (!use_n_t2()) return(NULL)
    p <- inv_params()
    if (!is.null(p$error)) return(NULL)
    lo <- if (inv_use_pomp()) 0 else inv_min()
    hi <- if (inv_use_pomp()) 100 else inv_max()
    beta_density_ci(p$alpha, p$beta, lo, hi, n_t2())
  })

  output$t2_ci_lo <- renderText({
    ci <- t2_ci_band()
    if (is.null(ci)) "—" else round(min(ci$y_lo), 3)
  })

  output$t2_ci_hi <- renderText({
    ci <- t2_ci_band()
    if (is.null(ci)) "—" else round(max(ci$y_hi), 3)
  })

  output$inv_plot <- renderPlot({
    p <- inv_params()
    if (!is.null(p$error)) return(NULL)
    ci_band <- t2_ci_band()
    if (inv_use_pomp()) {
      beta_plot_fn(p$alpha, p$beta, 0, 100, ci_band) +
        labs(x = "POMP score (% of maximum possible)")
    } else {
      lo <- inv_min(); hi <- inv_max()
      beta_plot_fn(p$alpha, p$beta, lo, hi, ci_band)
    }
  }, bg = "#16151f")
}

shinyApp(ui, server)
