## ----include = FALSE, warning=FALSE, message=FALSE----------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# Install locally
#  devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE )
# Build
#  setwd(R'(C:\Users\James.Thorson\Desktop\Git\dsem)'); devtools::build_rmd("vignettes/dynamic_factor_analysis.Rmd")

## ----setup, echo=TRUE, message=FALSE------------------------------------------
library(dsem)
library(MARSS)
library(ggplot2)
data( harborSealWA, package="MARSS")

# Define helper function
grab = \(x,name) x[which(names(x)==name)] 

# Define number of factors
# n_factors >= 3 doesn't seem to converge using DSEM or MARSS without penalties
n_factors = 2 

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7----------------------
# Load data
dat <- t(scale(harborSealWA[,c("SJI","EBays","SJF","PSnd","HC")]))

# DFA with 3 states; used BFGS because it fits much faster for this model
fit_MARSS <- MARSS( dat, 
                    model = list(m=n_factors),
                    form="dfa",
                    method="BFGS",
                    silent = TRUE )

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=3.5--------------------
# Plots states using all data
plot(fit_MARSS, plot.type="xtT")

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7----------------------
# Plot expectation for data using all data
plot(fit_MARSS, plot.type="fitted.ytT")

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7, warning=FALSE-------
# Add factors to data
tsdata = ts( cbind(harborSealWA[,c("SJI","EBays","SJF","PSnd","HC")]), start=1978)

# Scale and center (matches with MARSS does when fitting a DFA)
tsdata = scale( tsdata, center=TRUE, scale=TRUE )

# Define SEM
sem = "
  # Random-walk process for variables 
  SJF -> SJF, 1, NA, 1
  SJI -> SJI, 1, NA, 1
  EBays -> EBays, 1, NA, 1
  PSnd -> PSnd, 1, NA, 1
  HC -> HC, 1, NA, 1
"

# Initial fit
mydsem0 = dsem( tsdata = tsdata,
               covs = c("SJF, SJI, EBays, PSnd, HC"),
               sem = sem,
               family = rep("normal", 5),
               control = dsem_control( quiet = TRUE,
                                       run_model = FALSE ) )   

# fix all measurement errors at diagonal and equal
map = mydsem0$tmb_inputs$map
map$lnsigma_j = factor( rep(1,ncol(tsdata)) )

#
mydsem_full = dsem( tsdata = tsdata,
               covs = c("SJF, SJI, EBays, PSnd, HC"),
               sem = sem,
               family = rep("normal", 5),
               control = dsem_control( quiet = TRUE,
                                       map = map ) )

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7----------------------
plot_states = function( out,
                        vars=1:ncol(out$tmb_inputs$data$y_tj) ){
  # 
  xhat_tj = as.list(out$sdrep,report=TRUE,what="Estimate")$z_tj[,vars,drop=FALSE]
  xse_tj = as.list(out$sdrep,report=TRUE,what="Std. Error")$z_tj[,vars,drop=FALSE]

  # 
  longform = expand.grid( Year=time(tsdata), Var=colnames(tsdata)[vars] )
  longform$est = as.vector(xhat_tj)
  longform$se = as.vector(xse_tj)
  longform$upper = longform$est + 1.96*longform$se
  longform$lower = longform$est - 1.96*longform$se
  longform$data = as.vector(tsdata[,vars,drop=FALSE])
  
  # 
  ggplot(data=longform) +  #, aes(x=interaction(var,eq), y=Estimate, color=method)) +
    geom_line( aes(x=Year,y=est) ) +
    geom_point( aes(x=Year,y=data), color="blue", na.rm=TRUE ) +
    geom_ribbon( aes(ymax=as.numeric(upper),ymin=as.numeric(lower), x=Year), color="grey", alpha=0.2 ) + 
    facet_wrap( facets=vars(Var), scales="free", ncol=2 )
}
plot_states( mydsem_full )

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7----------------------
# Add factors to data
tsdata = harborSealWA[,c("SJI","EBays","SJF","PSnd","HC")]
newcols = array( NA,
                 dim = c(nrow(tsdata),n_factors),
                 dimnames = list(NULL,paste0("F",seq_len(n_factors))) )
tsdata = ts( cbind(tsdata, newcols), start=1978)

# Scale and center (matches with MARSS does when fitting a DFA)
tsdata = scale( tsdata, center=TRUE, scale=TRUE )

# Automated version
#sem = make_dfa( variables = c("SJI","EBays","SJF","PSnd","HC"),
#                n_factors = n_factors )
# Manual specification to show structure, using equations-and-lags interface
equations = "
  # Loadings of variables onto factors
  SJI = L11(0.1) * F1
  EBays = L12(0.1) * F1 + L22(0.1) * F2
  SJF = L13(0.1) * F1 + L23(0.1) * F2
  PSnd = L14(0.1) * F1 + L24(0.1) * F2
  HC = L15(0.1) * F1 + L25(0.1) * F2

  # random walk for factors
  F1 = NA(1) * lag[F1,1]
  F2 = NA(1) * lag[F2,1]

  # Unit variance for factors
  V(F1) = NA(1)
  V(F2) = NA(1)

  # Zero residual variance for variables
  V(SJI) = NA(0)
  V(EBays) = NA(0)
  V(SJF) = NA(0)
  V(PSnd) = NA(0)
  V(HC) = NA(0)
"
sem = convert_equations(equations)

# Initial fit
mydsem0 = dsem( tsdata = tsdata,
               sem = sem,
               family = c( rep("normal",5), rep("fixed",n_factors) ),
               estimate_delta0 = TRUE,
               control = dsem_control( quiet = TRUE,
                                       run_model = FALSE,
                                       gmrf_parameterization = "projection" ) )

# fix all measurement errors at diagonal and equal
map = mydsem0$tmb_inputs$map
map$lnsigma_j = factor( rep(1,ncol(tsdata)) )

# Fix factors to have initial value, and variables to not
map$delta0_j = factor( c(rep(NA,ncol(harborSealWA)-1), 1:n_factors) )

# Fix variables to have no stationary mean except what's predicted by initial value
map$mu_j = factor( rep(NA,ncol(tsdata)) )

# profile "delta0_j" to match MARSS (which treats initial condition as unpenalized random effect)
mydfa = dsem( tsdata = tsdata,
               sem = sem,
               family = c( rep("normal",5), rep("fixed",n_factors) ),
               estimate_delta0 = TRUE,
               control = dsem_control( quiet = TRUE,
                                       map = map,
                                       use_REML = TRUE,
                                       #profile = "delta0_j",
                                       gmrf_parameterization = "projection" ) )

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=3.5--------------------
# Plot estimated factors
plot_states( mydfa, vars=5+seq_len(n_factors) )

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7----------------------
# Plot estimated variables
plot_states( mydfa, vars=1:5 )

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7----------------------
# Extract internal settings
modmats <-  summary(fit_MARSS$model, silent=TRUE)

# Redefine defaults
modmats$V0 <- matrix(0, n_factors, n_factors )
modmats$x0 <- "unequal" 

# Refit
fit_MARSS2 = MARSS( dat, 
                    model = modmats,
                    silent = TRUE,
                    control = list( abstol = 0.001,
                                    conv.test.slope.tol = 0.01, 
                                    maxit = 1000 ))

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=3.5--------------------
# Plots states using all data
plot(fit_MARSS2, plot.type="xtT")

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7----------------------
# Compare likelihood for MARSS and DSEM
Table = c( "MARSS" = logLik(fit_MARSS), 
           "DSEM" = logLik(mydfa), 
           "MARSS_no_pen" = logLik(fit_MARSS2) )
knitr::kable( Table, digits=3)       

## ----echo=TRUE, message=FALSE, fig.width=7, fig.height=7----------------------
Table = cbind( "MARSS" = as.vector(fit_MARSS$par$Z),
       "DSEM" = grab(mydfa$opt$par,"beta_z"),
       "MARSS_no_pen" = as.vector(fit_MARSS2$par$Z) )
rownames(Table) = names(fit_MARSS$coef)[1:nrow(Table)]
knitr::kable( Table, digits=3)       

