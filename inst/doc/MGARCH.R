## ----include = FALSE, warning=FALSE, message=FALSE----------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
start_time = Sys.time()
# Install locally
#  devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE )
# Build
#  setwd(R'(C:\Users\James.Thorson\Desktop\Git\dsem)'); devtools::build_rmd("vignettes/spatial_diffusion.Rmd")

## -----------------------------------------------------------------------------
library(dsem)
set.seed(123)

# Specify settings
n_times = 100
n_vars = 3

# SD over time
sigF_t = seq( 0.1, 0.3, length = n_times )

# Simulate and apply time-varying SD
eps_tc = matrix( rnorm(n_times*n_vars), ncol = n_vars )
eps_tc = sweep( eps_tc, MARGIN = 1, FUN = "*", STAT = sigF_t )

## -----------------------------------------------------------------------------
# Define data including latent factor for heteroskedasticity
dat = data.frame(
  setNames( data.frame(eps_tc),letters[seq_len(n_vars)]),
  F = NA
)

# Define SEM using F as latent moderating variable
sem = "
  a <-> a, 0, F
  b <-> b, 0, F
  c <-> c, 0, F
  F <-> F, 0, sdF, 0.1
  F -> F, 1, NA, 1
"

# exploratory fit
fit1 = dsem(
  tsdata = ts(dat),
  sem = sem,
  estimate_mu = colnames(dat),
  control = dsem_control(
    use_REML = FALSE,
    gmrf_parameterization = "full",
    logscale_moderating_variance = TRUE,
    quiet = TRUE
  )
)

# Inspect estimates
summary(fit1)

## -----------------------------------------------------------------------------
# Define data including latent factor for heteroskedasticity and covariate
dat = data.frame(
  setNames( data.frame(eps_tc),letters[seq_len(n_vars)]),
  F = NA,
  slope = scale( seq_len(n_times), center = TRUE, scale = TRUE )
)

# Randomly simulate 10% missing data for covariate
dat$slope[ sample(seq_len(n_times), n_times/2) ] = NA

# Define SEM using F as latent moderating variable
# and slope as covariate for F
sem = "
  a <-> a, 0, F
  b <-> b, 0, F
  c <-> c, 0, F
  F <-> F, 0, sdF, 0.1
  slope <-> slope, 0, sd_slope
  slope -> slope, 1, NA, 1
  slope -> F, 0, beta
"

# confirmatory MGARCH
fit2 = dsem(
  tsdata = ts(dat),
  sem = sem,
  estimate_mu = colnames(dat),
  control = dsem_control(
    use_REML = FALSE,
    gmrf_parameterization = "full",
    logscale_moderating_variance = TRUE,
    quiet = TRUE
  )
)

# Inspect estimates
summary(fit2)

## ----fig.width=4, fig.height=4------------------------------------------------
# Bundle true and estimated time-series
Y = cbind(
  True = sigF_t,
  exp(predict(fit1)[,4]),
  exp(predict(fit2)[,4])
)

#
matplot( 
  x = seq_len(n_times), y = Y, type = "l", lty = "solid",
  col = c("black","red","blue"), xlab = "Time", 
  ylab = "SD for heteroskedasticity"
)
legend( "topleft", fill = c("black","red","blue"), bty = "n",
        legend = c("True", "Exploratory", "Confirmatory"))

## ----include = FALSE, warning=FALSE, message=FALSE----------------------------
run_time = Sys.time() - start_time

