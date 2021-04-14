##########################################################################################
#                                                                                        #
##              LGNB model: A practical statistical framework to combine                ##
##               commercial & survey data to model the spatio-temporal                  ##
##                          dynamics of marine harvested species                        ##
##                                    (Rufener et al.)                                  ##
#                                                                                        #
##########################################################################################


## The following script computes the results related to the consistency checks
## (i.e., validation of the integrated model). The model consistency consits of two steps:
## 1) Consistency of random and fixed-effect parameters (statistically and visually)
## 2) Evaluate model's goodness-of-fit based on residual normality (statistically and visually)


# last update: April 2021


#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><


# Load libraries
#~~~~~~~~~~~~~~~~~~
library(TMB)
library(tmbstan)
library(ggplot2)


# Load dynamic library of the LGNB model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("~/LGNB/src")
dyn.load(dynlib("LGNB"))


#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   1) Consistency of fixed and random effect parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1.1) Get results
#~~~~~~~~~~~~~~~~~~~~

#setwd("~/LGNB/src")

## Survey model
load("res_m1_A3_survey_No_alpha_YearQuarter_.RData")
stopifnot(levels(datatot$Data)[1] == "survey")
obj1 <- obj; sdr1 <- env1$sdr #change "env1" by the model number that was ran; default is M1 in LGNB_Rmodel.R script (hence, env1)


## Integrated model
load("res_m1_A3_both_No_alpha_YearQuarter_.RData")
stopifnot(levels(datatot$Data)[1] == "survey")
obj2 <- obj; sdr2 <- env1$sdr #change "env1" by the model number that was ran; default is M1 in LGNB_Rmodel.R script (hence, env1)


# 1.2) Evaluate likelihood function and gradients
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Likelihood function
obj1$fn() #survey model
obj2$fn() #integrated model

## Gradients
obj1$gr() #survey model
obj2$gr() #integrated model



# 1.3) Extract shared parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 1.3.1) First, extract the parameters from each mdoel

pl1 <- as.list(sdr1,"Est") #survey model parameters
pl2 <- as.list(sdr2,"Est") #integrated model parameters
#names(pl1); names(pl2)

par1 <- sdr1$par.fixed # Parameters from first model (survey alone)


## 1.3.2) Then, extract the parameters that are shared between the two models

# Input:
# - obj1: Model object of survey alone
# - pl1:  Estimate list of L1(theta1)
# - pl2:  Estimate list of L2(theta1, theta2)
# Output:
# Estimate of theta1 based on L2
extractBaseParameters <- function(obj1, pl1, pl2, all=FALSE) {
  npar1 <- sapply(pl1,length) #Number of survey parameters by component
  plnew <- Map(head, pl2, npar1) #Extract survey parameters from second fit
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



# 1.4) Re-evaluate the survey likelihood function based on extracted parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test that results of second fit doesn't contradict first fit
f1 <- as.numeric(obj1$fn(par1)) #Evaluate survey likelihood from survey model parameters
f2 <- as.numeric(obj1$fn(par2)) #Evaluate survey likelihood from integrated model parameters



# 1.5) Test statistics
#~~~~~~~~~~~~~~~~~~~~~~~~
# Checks that the estimates from the integrated model are within the multivariate
# estimates based on the parameters estimated by the survey model
# p > 0.05 implies that the parameters of the integrated model are within the 95% confidence
# region of the survey model.

## For the fixed-effect parameters
df <- sum(par1 != 0) ## Get degrees of freedom - recall we use '0' for betas that are not estimated by obj1
fixed <- 1 - pchisq( 2 * (f2-f1), df=df )
fixed


## Similar test now including random effects
par2.all <- extractBaseParameters(obj1, pl1, pl2, all=TRUE)
f1.all <- obj1$env$f(obj1$env$last.par.best) #Best evaluated parameters
f2.all <- obj1$env$f(par2.all)
df.all <- df + length(obj1$env$random)
random <- 1 - pchisq( 2 * (f2.all-f1.all), df=df.all )
random


## Compare estimates visually
i <- (par1 != 0) ## consider non missing parameters

par(mfrow=c(1,2))
plot(par1[i], par2[i], main="fixed effects"); abline(0,1, col="red")
legend("topleft", 
       legend = paste("p = ",round(fixed,6)),
       cex = 1.2)


plot(as.vector(pl1$eta_density), as.vector(pl2$eta_density), main="random effects"); abline(0,1, col="red")
legend("topleft", 
       legend = paste("p = ", round(random,6)),
       cex = 1.2)


# # Plotting in ggplot....
# dfrandom <- data.frame(pl1 = as.vector(pl1$eta_density), pl2 = as.vector(pl2$eta_density)) #Random effects
# dfixed <- data.frame(par1 = par1[i], par2 = par2[i]) #Fixed effects
# 
# ## random-effects
# ggplot(dfrandom, aes(x=pl1,y=pl2)) +
#   geom_point(alpha = 0.1, size=2.5) +
#   geom_abline(col="#0073C2FF",size=1.4, linetype = "dashed") +
#   theme_bw() +
#   ggtitle("Age 2") +
#   ylab("Random Effects (survey data)") +
#   xlab("Random Effects (integrated data)") +
#   theme(panel.border = element_blank(),
#         axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
#         axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"),
#         
#         plot.title = element_text(hjust = 0.5, margin=margin(b=15),size=18,face="bold"),
#         
#         axis.text.x = element_text(face="bold",size=11),
#         axis.text.y = element_text(size=11,face="bold"),
#         axis.title.y = element_text(margin=margin(t=0,r=20,b=0,l=0),size=12,face="bold"),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=0,l=0),size=12,face="bold"),
#         #axis.line = element_line(size=1, colour = "black"),
#         
#         legend.position = "right",
#         
#         plot.margin = unit(c(1,1,1,1),"cm"))
# 
# 
# ## Fixed-effects
# ggplot(dfixed, aes(x=par1,y=par2)) +
#   geom_point(alpha = 0.3, size=3.5) +
#   geom_abline(col="#0073C2FF",size=1.4, linetype = "dashed") +
#   theme_bw() +
#   ggtitle("Age 2") +
#   ylab("Fixed Effects (survey data)") +
#   xlab("Fixed Effects (integrated data)") +
#   theme(panel.border = element_blank(),
#         axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
#         axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"),
#         
#         plot.title = element_text(hjust = 0.5, margin=margin(b=15),size=18,face="bold"),
#         
#         axis.text.x = element_text(face="bold",size=11),
#         axis.text.y = element_text(size=11,face="bold"),
#         axis.title.y = element_text(margin=margin(t=0,r=20,b=0,l=0),size=12,face="bold"),
#         axis.title.x = element_text(margin=margin(t=20,r=0,b=0,l=0),size=12,face="bold"),
#         #axis.line = element_line(size=1, colour = "black"),
#         
#         legend.position = "right",
#         
#         plot.margin = unit(c(1,1,1,1),"cm"))




#~~~~~~~~~~~~~~~~~~~~~~~~~~
#   2) Goodness-of-fit 
#~~~~~~~~~~~~~~~~~~~~~~~~~~

# In this case based on model 'obj2' with survey and commercial data (integrated model)
# Computing residuals from midex-effect models is not a trivial task.
# Here, we compute simulation-based residuals, which are based on a sample from the
# posterior distribution of the latent variables given the data and with the parameters 
# replaced by their MLE.


## Get names of random effects in the model
random <- unique(names(obj2$env$par[obj2$env$random]))


## Get (logical) non random effects indices
fixed <- !(names(pl2) %in% random)


## Fix non-random parameters to their estimated values
map <- lapply(pl2[fixed], function(x) factor(rep(NA,length(x))) )


## Construct corresponding new function object
obj <- MakeADFun(obj2$env$data, pl2, map=map, DLL="LGNB")


## Run MCMC to get posterior sample of random effects given data
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
obj.sim <- MakeADFun(obj2$env$data, pl.sim, map=map, DLL="LGNB")

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



# Plot data-specific residuals
dfres <- data.frame(res = residual, Data = datatot$Data)


ggplot(dfres, aes(x=res, color=Data,fill=Data)) +
  geom_histogram(alpha=0.2, position="identity",size=1.2) +
  scale_colour_manual(values=c("#EFC000FF", "#0073C2FF")) +
  scale_fill_manual(values=c("#EFC000FF", "#0073C2FF")) +
  geom_vline(xintercept=0, linetype="dashed", color = "black",size=1.2) +
  theme_bw() +
  xlab("Standardized residuals") +
  ylab("Frequency") +
  ggtitle(AGE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(panel.border = element_blank(),
        
        axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"),
        
        plot.title = element_text(hjust = 0.5, margin=margin(b=15),size=18,face="bold"),
        
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(size=11,face="bold"),
        axis.title.y = element_text(margin=margin(t=0,r=20,b=0,l=0),size=12,face="bold"),
        axis.title.x = element_text(margin=margin(t=20,r=0,b=0,l=0),size=12,face="bold"),
        #axis.line = element_line(size=1, colour = "black"),
        
        legend.position = "right",
        #legend.position = "none",
        
        plot.margin = unit(c(1,1,1,1),"cm"))
