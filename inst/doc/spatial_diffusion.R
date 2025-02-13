## ----include = FALSE, warning=FALSE, message=FALSE----------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# Install locally
#  devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE )
# Build
#  setwd(R'(C:\Users\James.Thorson\Desktop\Git\dsem)'); devtools::build_rmd("vignettes/spatial_diffusion.Rmd")

## ----simulate, echo=TRUE, message=FALSE---------------------------------------
library(igraph)
library(Matrix)

# Simulation
adjacency_graph = make_graph( ~ A - B - C - D - E )
A = as.matrix( adjacency_graph )

# Diffusion rate
Dprime = 1 * A
diag(Dprime) = -1 * rowSums(Dprime)

# Movement transition matrix
M = expm( Dprime )

# set seed for reproducibility
set.seed(101)

# Simulate densities
n_times = 100
n_burnin = 100
x_ti = matrix( NA, nrow=n_times+n_burnin, ncol = nrow(M) )
x_ti[1,] = rnorm(n=nrow(M), mean = 0, sd = 1 )
for( t in 2:nrow(x_ti) ){
  x_ti[t,] = (x_ti[t-1,] %*% M)[1,] + rnorm(n=nrow(M), mean = 0, sd = 0.1)
}

# Subset to times after burn-in
x_ti = x_ti[ n_burnin+seq_len(n_times), ]

## ----fit, echo=TRUE, message=FALSE--------------------------------------------
library(dsem)

# Specify SEM
sem = "
  # Spatial correlation
  A -> B, 0, d0
  B -> C, 0, d0
  C -> D, 0, d0
  D -> E, 0, d0
  E -> D, 0, d0
  D -> C, 0, d0
  C -> B, 0, d0
  B -> A, 0, d0

  # Spatio-temporal diffusion
  A -> B, 1, d
  B -> C, 1, d
  C -> D, 1, d
  D -> E, 1, d
  E -> D, 1, d
  D -> C, 1, d
  C -> B, 1, d
  B -> A, 1, d

  # Self-limitation
  A -> A, 1, rho
  B -> B, 1, rho
  C -> C, 1, rho
  D -> D, 1, rho
  E -> E, 1, rho
"

# Fit
colnames(x_ti) = c("A","B","C","D","E")
fit = dsem(
  tsdata = ts(x_ti),
  sem = sem
)

## ----predict, echo=TRUE, message=FALSE----------------------------------------
# Calculate total effect
effect = total_effect(fit, n_lags = 3)

# Calculate predicted movement
Mhat = array( subset(effect,lag==1)$total_effect, dim(M) )
dimnames(Mhat) = dimnames(M)

# Display predicted movement
knitr::kable( Mhat, digits=2)

## ----display, echo=TRUE, message=FALSE----------------------------------------
knitr::kable( as.matrix(M), digits=2)

