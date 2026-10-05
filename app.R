library(shiny)
library(MASS)
library(mvtnorm)
library(mclust)
library(teigen)

pt_col <- function(type) ifelse(type == 1, "blue", ifelse(type == 2, "red", "black"))
pt_pch <- function(type) ifelse(type == 1, 1, ifelse(type == 2, 4, 3))
cpal   <- c("blue", "red", "green", "orange", "purple",
            "brown", "cyan3", "magenta", "gray40")

# ---- 1. data generation ---------------------------------------------------
gen_clean <- function(n, pi1, d, rho, seed) {
  set.seed(seed)
  mus  <- list(c(0, -d), c(0, d))
  covs <- list(matrix(c(1, -rho, -rho, 1), 2), matrix(c(1, rho, rho, 1), 2))
  z <- sample(1:2, n, replace = TRUE, prob = c(pi1, 1 - pi1))
  X <- t(sapply(z, function(g) mvrnorm(1, mus[[g]], covs[[g]])))
  list(X = X, type = z)
}

add_noise <- function(X, type, n_noise, margin) {
  if (n_noise == 0) return(list(X = X, type = type))
  r1 <- range(X[, 1]) + c(-margin, margin)
  r2 <- range(X[, 2]) + c(-margin, margin)
  noise <- cbind(runif(n_noise, r1[1], r1[2]), runif(n_noise, r2[1], r2[2]))
  list(X = rbind(X, noise), type = c(type, rep(0, n_noise)))
}

# ---- 2. model fitting -----------------------------------------------------
fit_models <- function(X, Gmax, dfstart, maxit2) {
  list(
    g = Mclust(X, verbose = FALSE),   # default G = 1:9, as in the script
    t = teigen(X, Gs = 1:Gmax, models = "all", dfstart = dfstart,
               eps = c(0.001, 0.1), maxit = c(200, maxit2), verbose = FALSE)
  )
}

# ---- 3. parameter extraction (t params back-transformed to original scale) --
params_gauss <- function(f, p = 2) {
  list(G = f$G, pi = f$parameters$pro,
       mu = matrix(f$parameters$mean, nrow = p),
       S = array(f$parameters$variance$sigma, c(p, p, f$G)), nu = NULL)
}

params_t <- function(f, X) {
  p <- ncol(X); G <- f$G
  sx <- apply(X, 2, sd); mx <- colMeans(X)
  mu <- matrix(f$parameters$mean, nrow = G)
  mu <- t(sweep(sweep(mu, 2, sx, "*"), 2, mx, "+"))
  S <- array(f$parameters$sigma, c(p, p, G))
  for (g in seq_len(G)) S[, , g] <- S[, , g] * outer(sx, sx)
  list(G = G, pi = f$parameters$pig, mu = mu, S = S,
       nu = rep_len(f$parameters$df, G))
}

# ---- 4. plotting ----------------------------------------------------------
draw_points <- function(X, type, main) {
  plot(X[, 1], X[, 2], col = pt_col(type), pch = pt_pch(type), cex = 0.6,
       xlab = "x1", ylab = "x2", main = main)
}

draw_fit <- function(X, type, pr, main) {
  draw_points(X, type, main)
  xs <- seq(min(X[, 1]) - 1, max(X[, 1]) + 1, length.out = 200)
  ys <- seq(min(X[, 2]) - 1, max(X[, 2]) + 1, length.out = 200)
  grid <- expand.grid(x = xs, y = ys)
  for (g in seq_len(pr$G)) {
    d <- if (is.null(pr$nu)) {
      dmvnorm(grid, mean = pr$mu[, g], sigma = pr$S[, , g])
    } else {
      dmvt(grid, delta = pr$mu[, g], sigma = pr$S[, , g], df = pr$nu[g], log = FALSE)
    }
    contour(xs, ys, matrix(d, nrow = length(xs)), add = TRUE,
            col = cpal[(g - 1) %% length(cpal) + 1], lwd = 1.5, drawlabels = FALSE)
  }
  legend("topright", legend = paste("Comp", seq_len(pr$G)),
         col = cpal[(seq_len(pr$G) - 1) %% length(cpal) + 1], lwd = 2, bty = "n", cex = 0.7)
}

# ---- 5. tables ------------------------------------------------------------
est_table <- function(pr, model) {
  data.frame(Model = model, Comp = seq_len(pr$G), pi = round(pr$pi, 3),
             mu1 = round(pr$mu[1, ], 2), mu2 = round(pr$mu[2, ], 2),
             nu = if (is.null(pr$nu)) NA else round(pr$nu, 2),
             S11 = round(pr$S[1, 1, ], 2), S12 = round(pr$S[1, 2, ], 2),
             S22 = round(pr$S[2, 2, ], 2))
}

# ---- UI -------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Gaussian vs t mixture models under contamination"),
  sidebarLayout(
    sidebarPanel(
      h4("Data"),
      sliderInput("n", "Observations (n)", 50, 1000, 200, step = 50),
      sliderInput("pi1", "Mixing proportion of cluster 1", 0.1, 0.9, 0.3, step = 0.05),
      sliderInput("n_noise", "Uniform noise points", 0, 200, 30, step = 5),
      sliderInput("margin", "Noise margin beyond data range", 0, 5, 2, step = 0.5),
      numericInput("seed", "Random seed", 44),
      h4("Fitting"),
      numericInput("Gmax_clean", "teigen max G (clean data)", 3, min = 1, max = 9),
      numericInput("Gmax_contam", "teigen max G (contaminated data)", 9, min = 1, max = 9),
      numericInput("dfstart", "teigen starting df", 15, min = 3),
      actionButton("run", "Simulate & fit", class = "btn-primary")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data", plotOutput("p_data", height = 420)),
        tabPanel("Clean data fit", plotOutput("p_clean", height = 420)),
        tabPanel("Contaminated data fit", plotOutput("p_contam", height = 420)),
        tabPanel("Estimates",
                 h4("True parameters"), tableOutput("t_truth"),
                 h4("Clean data"), tableOutput("t_clean"),
                 h4("Contaminated data"), tableOutput("t_contam"),
                 p(tags$em("For the t mixture, S is the scale matrix; covariance = \u03bd/(\u03bd\u22122)\u00b7S (\u03bd > 2).")),
                 p(tags$em("Gaussian: S is the covariance matrix.")))
      )
    )
  )
)

# ---- server ---------------------------------------------------------------
server <- function(input, output, session) {
  
  sim <- eventReactive(input$run, {
    cl <- gen_clean(input$n, input$pi1, 3, 0.5, input$seed)   # means (0, -3), (0, 3); correlation -0.5, +0.5
    withProgress(message = "Fitting models...", value = 0, {
      incProgress(0.1, detail = "clean data")
      f1 <- fit_models(cl$X, input$Gmax_clean, input$dfstart, 500)
      # noise is drawn after the clean fit with no reseeding, as in the original script
      ct <- add_noise(cl$X, cl$type, input$n_noise, input$margin)
      incProgress(0.4, detail = "contaminated data")
      f2 <- fit_models(ct$X, input$Gmax_contam, input$dfstart, 700)
    })
    list(
      clean  = list(X = cl$X, type = cl$type, fits = f1),
      contam = list(X = ct$X, type = ct$type, fits = f2),
      truth  = data.frame(Comp = 1:2, pi = c(input$pi1, 1 - input$pi1),
                          mu1 = 0, mu2 = c(-3, 3),
                          S11 = 1, S12 = c(-0.5, 0.5), S22 = 1)
    )
  }, ignoreNULL = FALSE)
  
  fit_plot <- function(sc) {
    par(mfrow = c(1, 2))
    draw_fit(sc$X, sc$type, params_gauss(sc$fits$g),
             paste0("Gaussian mixture (G = ", sc$fits$g$G, ")"))
    draw_fit(sc$X, sc$type, params_t(sc$fits$t, sc$X),
             paste0("t mixture (G = ", sc$fits$t$G, ")"))
  }
  
  output$p_data <- renderPlot({
    s <- sim(); par(mfrow = c(1, 2))
    draw_points(s$clean$X, s$clean$type, "Clean data")
    draw_points(s$contam$X, s$contam$type, "With uniform noise (black +)")
  })
  output$p_clean  <- renderPlot(fit_plot(sim()$clean))
  output$p_contam <- renderPlot(fit_plot(sim()$contam))
  
  output$t_truth <- renderTable(sim()$truth)
  est <- function(sc) rbind(est_table(params_gauss(sc$fits$g), "Gaussian"),
                            est_table(params_t(sc$fits$t, sc$X), "t"))
  output$t_clean  <- renderTable(est(sim()$clean))
  output$t_contam <- renderTable(est(sim()$contam))
}

shinyApp(ui, server)