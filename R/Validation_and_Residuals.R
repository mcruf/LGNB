library(TMB)



dyn.load(dynlib("model"))


## Get model 1 results (survey model)
load("results_WBScod_m1_A2_survey_No_One_noprofile_.RData")
stopifnot(levels(datatot$Data)[1] == "survey")
obj1 <- obj; sdr1 <- env1$sdr


## Get model 2 results (integrated model)
load("results_WBScod_m1_A2_both_No_One_noprofile_.RData")
stopifnot(levels(datatot$Data)[1] == "survey")
obj2 <- obj; sdr2 <- env1$sdr


obj1$fn()
obj1$gr()
obj2$fn()
obj2$gr()

pl1 <- as.list(sdr1,"Est")
pl2 <- as.list(sdr2,"Est")

## Parameters from first model (survey alone)
par1 <- sdr1$par.fixed

## Input:
## - obj1: Model object of survey alone
## - pl1:  Estimate list of L1(theta1)
## - pl2:  Estimate list of L2(theta1, theta2)
## Output:
## Estimate of theta1 based on L2
extractBaseParameters <- function(obj1, pl1, pl2, all=FALSE) {
  npar1 <- sapply(pl1,length) ## Number of survey parameters by component
  plnew <- Map(head, pl2, npar1) ## Extract survey parameters from second fit
  applyMap <- function(parameters, map) {
    param.map <- lapply(names(map),
                        function(nam)
                        {
                          TMB:::updateMap(parameters[[nam]], map[[nam]])
                        })
    parameters[names(map)] <- param.map
    parameters
  }
  par2 <- unlist(applyMap(plnew, obj1$env$map))
  if (!all) par2 <- par2[-obj1$env$random]
  par2
}

## Parameters comparable to 'par1' obtained from expanded model + data
par2 <- extractBaseParameters(obj1, pl1, pl2)

## Test that results of second fit doesn't contradict first fit
f1 <- as.numeric( obj1$fn(par1) )
f2 <- as.numeric( obj1$fn(par2) )

## Get degrees of freedom - recall we use '0' for betas that are not estimated by obj1
df <- sum(par1 != 0)
1 - pchisq( 2 * (f2-f1), df=df )

## Similar test now including random effects
par2.all <- extractBaseParameters(obj1, pl1, pl2, all=TRUE)
f1.all <- obj1$env$f(obj1$env$last.par.best)
f2.all <- obj1$env$f(par2.all)
df.all <- df + length(obj1$env$random)
1 - pchisq( 2 * (f2.all-f1.all), df=df.all )

## Compare estimates visually
i <- (par1 != 0) ## consider non missing parameters
plot(par1[i], par2[i]); abline(0,1, col="red")

## Compare spatial random effects
plot(as.vector(pl1$eta_density), as.vector(pl2$eta_density)); abline(0,1, col="red")



########################################################################
## Residuals
## In this case based on model 'obj2' with survey and commercial data
########################################################################

## Get names of random effects in the model
random <- unique(names(obj2$env$par[obj2$env$random]))
## Get (logical) non random effects indices
fixed <- !(names(pl2) %in% random)

## Fix non-random parameters to their estimated values
map <- lapply(pl2[fixed], function(x) factor(rep(NA,length(x))) )

## Construct corresponding new function object
obj <- MakeADFun(obj2$env$data, pl2, map=map, DLL="model")

## Run MCMC to get posterior sample of random effects given data
library(tmbstan)
set.seed(123)
qw <- tmbstan(obj, chains=1, iter=100, warmup=99) ## FIXME: Increase ???
s <- extract(qw)

## Plugin simulation in parameter list (copy and override)
pl.sim <- pl2
for (nm in random) {
  pl.sim[[nm]][] <- s[[nm]]
}

## Once again construct a function object, this time using simulated
## random effects
obj.sim <- MakeADFun(obj2$env$data, pl.sim, map=map, DLL="model")

## Luckily the model template has 'REPORT' statements for what we need (mu and var)
rep <- obj.sim$report()
mu <- exp(rep$log_mu)                 ## Mean
var <- exp(rep$log_var_minus_mu) + mu ## Variance
obs <- obj.sim$env$data$response      ## Observed

## Transform nbinom parameters (mu, var) to R's parameterization (mu, size)
## var = mu + mu^2/size
size <- mu^2 / (var - mu)

## Transform observations with CDF and randomize point probabilities
Fx <- pnbinom(obs, mu=mu, size=size)
px <- dnbinom(obs, mu=mu, size=size)
u <- runif(length(Fx))
residual <- qnorm(Fx - u * px)

## Finally:
qqnorm(residual)
abline(0, 1)
ks.test(residual, "pnorm")
