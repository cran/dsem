## ----setup, echo=FALSE, cache=FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = FALSE,
  error = FALSE,
  message = FALSE,
  warning = FALSE
)
# Install locally
#  devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE )
# Build and PDF
#  setwd(R'(C:\Users\James.Thorson\Desktop\Git\dsem)'); devtools::build_rmd("vignettes/model-description.Rmd"); rmarkdown::render( "vignettes/model-description.Rmd", rmarkdown::pdf_document())
#
# Recommended workflow:
#  * Open in Rstudio and knit using button there

## ----echo=TRUE, eval=FALSE----------------------------------------------------
# x -> x, 1, ar1
# x <-> x, 0, sd

## ----echo=TRUE, eval=FALSE----------------------------------------------------
# A -> B, 0, b_AB
# B -> C, 1, b_BC

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
dsem = " 
x -> x, 1, ar1, 0.8
x <-> x, 0, sd, 1
"

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
# Load package
library(dsem)

# call dsem without estimating parameters
out = dsem(
  tsdata = ts(data.frame( x = rep(1,10) )),
  sem = dsem,
  control = dsem_control(
    run_model = FALSE, 
    quiet = TRUE
  )
)

# Extract covariance
Sigma1 = solve(as.matrix(out$obj$report()$Q_kk))
plot( x=1:10, y = diag(Sigma1), xlab="time", 
      ylab="Marginal variance", type="l", 
      ylim = c(0,max(diag(Sigma1))))

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
# call dsem without estimating parameters
out = dsem(
  tsdata = ts(data.frame( x = rep(1,10) )),
  sem = dsem,
  control = dsem_control(
    run_model = FALSE, 
    quiet = TRUE, 
    constant_variance = "marginal"
  )
)

# Extract covariance
Sigma2 = solve(as.matrix(out$obj$report()$Q_kk))
plot( x=1:10, y = diag(Sigma2), xlab="time", 
      ylab="Marginal variance", type="l", 
      ylim = c(0,max(diag(Sigma1))))

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
#
dsem = "
  # Factor follows random walk with unit variance
  F <-> F, 0, NA, 1
  F -> F, 1, NA, 1
  # Loadings on two manifest variables
  F -> x, 0, b_x, 1
  F -> y, 0, b_y, 1
  # No residual variance for manifest variables
  x <-> x, 0, NA, 0
  y <-> y, 0, NA, 0
"
data = data.frame( 
  x = rnorm(10),
  y = rnorm(10),
  F = rep(NA,10)
)

# call dsem without estimating parameters
out = dsem(
  tsdata = ts(data),
  sem = dsem,
  family = c("normal","normal","fixed"),
  control = dsem_control(
    run_model = FALSE, 
    quiet = TRUE,
    gmrf_parameterization = "projection"
  )
)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
# Extract covariance
library(Matrix)
IminusRho_kk = out$obj$report()$IminusRho_kk
G_kk = out$obj$report()$Gamma_kk
Q_kk = t(IminusRho_kk) %*% t(G_kk) %*% G_kk %*% IminusRho_kk

# Display eigenvalues
eigen(Q_kk)$values

