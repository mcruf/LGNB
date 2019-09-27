##########################################################################################
#                                                                                        #
##              LGNB model: A practical statistical framework to combine                ##
##               commercial & survey data to model the spatio-temporal                  ##
##                          dynamics of marine harvested species                        ##
##                                    (Rufener et al.)                                  ##
#                                                                                        #
##########################################################################################

# last update: September 2019



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 1: Default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
STOCK <- c("WBS","KAT")[1] # Specify to which stock the data should be subsetted (Western Baltic Sea or Kattegat)
INCLUDE  <- c("commercial", "survey", "both") [3] # Specify the desired input data for the model; default is commercial data-
RESPONSE <- c("SizeGroup","AgeGroup","Cohort")[2] #Choose whether the model is applied for each SizeGroup, AgeGroup, or on a Cohort basis (when Nage).Default is set to SizeGroup.
ALPHA <- c("No","Single","Multi")[3] # Define how many alpha parameters should be computed
SUPPORT_AREA <- c("One","Several")[2] #Choose whether to use one or several support areas to describe the commercial fisheries data; Survey data is described by only single support area

# @ ALPHA = "No" -> Models without alpha parameter; can be applied either to model with one or several support areas.
# @ ALPHA = "Single I" -> Models with one alpha parameter for each data soruce; Used when commercial data are described by only one support area for the entire time series; Also applied when using survey data;
# @ ALPHA = "Multi" -> Models with several alpha parameters to describe commercial data (when using several support areas); Does NOT apply to survey data;


# We go from the assumption that survey data has ALWAYS only one support area; 
# one can choose whether to estimate the alpha-parameter or not. 
if(INCLUDE == "survey"){
  SUPPORT_AREA <- "One" 
  ALPHA <- c("No","Single")[1] 
}


# Specify the model structure; default model is m2 (see lines ... to ...) for either size groups, age groups or cohort
# Default size group is 5 (S5), age group is 3 (A3) and cohort from 2005.
if(RESPONSE == "SizeGroup"){
  MODEL_CONFIG <- "m2_S5" #Default model and SizeGroup 5
}else if(RESPONSE == "AgeGroup"){
  MODEL_CONFIG <- "m2_A3" #Default model and AgeGroup 3
} else if(RESPONSE == "Cohort")
  MODEL_CONFIG <- "m1_2005" #Default model when applied on cohort-basis (2005 cohort) 

MODEL_CONFIG <- strsplit(MODEL_CONFIG, "_")[[1]]
MODEL_FORMULA <- MODEL_CONFIG[1]

if(RESPONSE == "SizeGroup"){
  SIZE <- MODEL_CONFIG[2]
} else if(RESPONSE == "AgeGroup"){
  AGE <- MODEL_CONFIG[2]
} else if(RESPONSE == "Cohort")
  YEARCLASS <- MODEL_CONFIG[2]



# For scripting (Useful when running on a HPC)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)
stopifnot(INCLUDE %in% c("commercial", "survey", "both"))



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 2: Load data files & R packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2.1) Load helper functions and R libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("C:/Users/mruf/Documents/LGNB/R/utilities.R")

#devtools::install_github("kaskr/gridConstruct",subdir="gridConstruct") # To install the gridConstruct package
mLoad(raster,rgeos,maptools,maps,data.table,dplyr,TMB,sp,
      DATRAS,gridConstruct,rgdal,geosphere,devtools,plyr,fields,forcats)



# 2.2.1) Load commercial fisheries data (fishery-depedent)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comFULL <- readRDS("C:/Users/mruf/Documents/LGNB/Data/commercial.rds")
comFULL$stock <- ifelse(comFULL$Area=="21","KAT","WBS")


# Subset data for a particular stock and/or time frame
commercial <- subset(comFULL, stock == STOCK) #Setting stock based on ICES area (KAT=21, WBS=22-24)


# Drop unused factor levels
commercial[,c("Month","Year","Quarter","Area","Sediment","Metiers","Data","Haul_ID","TimeYear","VE_LENcat")] <- lapply(commercial[,c("Month","Year","Quarter","Area","Sediment","Metiers","Data","Haul_ID","TimeYear","VE_LENcat")], factor)



# 2.2.2) Load scientific survey data (Fishery~indepdendent)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
survey <- readRDS("C:/Users/mruf/Documents/LGNB/Data/survey.rds")
survey$stock <- ifelse(survey$Area=="21","KAT","WBS") #Setting stock based on ICES area (KAT=21, WBS=22-24)


# Subset data for a particular time frame
survey <- filter(survey, stock == STOCK)


# Drop unused factor levels
survey[,c("Haul_ID","Month","Year","Quarter","Area","Data","TimeYear","Sediment")] <- lapply(survey[,c("Haul_ID","Month","Year","Quarter","Area","Data","TimeYear","Sediment")], factor)



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 3: Binding survey and commercial data files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Binding depends whether response variable is on a cohort basis or AgeGroup/SizeGroup basis


# 3.1.1) Cohort-basis
#~~~~~~~~~~~~~~~~~~~~~~~
# First we need to create a dataframe in such way that both datasets have the same timelevels;
# This is VERY important, as uneven timelevels will cause problems for the AR1 process.

if(RESPONSE == "Cohort"){
  df_cohort  <- extractCohortLevels(commercial,survey,yearclass = as.numeric(YEARCLASS)) #Df with equal time-steps
  NageGroup  <- length(grep("Age_", names(commercial), value = TRUE))
  tim        <- (as.numeric(YEARCLASS) + NageGroup)-1 #
  df_cohort  <- subset(df_cohort, Year %in% c(YEARCLASS:tim))
  
  cohort_com <- extract_cohort_quarter(commercial,df_cohort) #Extract cohort for commercial data
  cohort_sur <- extract_cohort_quarter(survey,df_cohort) #Extract cohort for survey data
  
  cohort_com <- transform(cohort_com, HLID=Haul_ID, HaulDur=HaulDuration_hours)
  cohort_sur <- transform(cohort_sur, latStart=lat, lonStart=lon, latEnd=lat, lonEnd=lon, HLID=Haul_ID)
  
  
  # 3.1.2) AgeGroup-basis
  #~~~~~~~~~~~~~~~~~~~~~~~~
} else if (RESPONSE == "AgeGroup"){
  age_com   <- transform(commercial, HLID=Haul_ID, HaulDur=HaulDuration_hours)
  age_sur   <- transform(survey, latStart=lat, lonStart=lon, latEnd=lat, lonEnd=lon,HLID=Haul_ID)
  
  
  # 3.1.3) SizeGroup-basis
  #~~~~~~~~~~~~~~~~~~~~~~~~
} else if(RESPONSE == "SizeGroup"){
  size_com <- transform(commercial, HLID=Haul_ID, HaulDur=HaulDuration_hours)
  size_sur <- transform(survey, latStart=lat, lonStart=lon, latEnd=lat, lonEnd=lon,HLID=Haul_ID)
  
}



# 3.2) Bind both datasets into a single data frame
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(RESPONSE == "Cohort" & INCLUDE =="both"){
  datatot <- mybind(cohort_com, cohort_sur)
} else if(RESPONSE == "Cohort" & INCLUDE == "survey"){
  datatot <- cohort_sur
} else if(RESPONSE == "Cohort" & INCLUDE == "commercial"){
  datatot <- cohort_com
} else if (RESPONSE == "AgeGroup" & INCLUDE == "both"){
  datatot <- mybind(age_com, age_sur)
} else if (RESPONSE == "AgeGroup" & INCLUDE == "survey") {
  datatot <- age_sur
} else if (RESPONSE == "AgeGroup" & INCLUDE == "commercial"){
  datatot <- age_com
} else if (RESPONSE == "SizeGroup" & INCLUDE == "both"){
  datatot <- mybind(size_com, size_sur)
} else if (RESPONSE == "SizeGroup" & INCLUDE == "survey"){
  datatot <- size_sur
} else if(RESPONSE == "SizeGroup" & INCLUDE == "commercial")
  datatot <- size_com




# 3.3) Define response variable for model applied on AgeGroup our SizeGroup basis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 3.3.1) For AgeGroups
#~~~~~~~~~~~~~~~~~~~~~~
if(RESPONSE == "AgeGroup"){
  if(AGE == "A0"){
    datatot$Response <- as.numeric(paste(datatot$Age_0))
  } else if(AGE == "A1"){
    datatot$Response <- as.numeric(paste(datatot$Age_1))
  } else if (AGE == "A2") {
    datatot$Response <- as.numeric(paste(datatot$Age_2))
  } else if (AGE == "A3"){
    datatot$Response <- as.numeric(paste(datatot$Age_3))
  } else if (AGE == "A4"){
    datatot$Response <- as.numeric(paste(datatot$Age_4))
  } else if (AGE == "A5"){
    datatot$Response <- as.numeric(paste(datatot$Age_5))  
  } else if (AGE == "A6"){
    datatot$Response <- as.numeric(paste(datatot$Age_6))
  } else if (AGE == "A7"){
    datatot$Response <- as.numeric(paste(datatot$Age_7))
  } else if (AGE == "A8"){
    datatot$Response <- as.numeric(paste(datatot$Age_8))
  }
}


# 3.3.2) For SizeGroups
#~~~~~~~~~~~~~~~~~~~~~~~
if(RESPONSE == "SizeGroup"){
  if(SIZE == "S1"){
    datatot$Response <- as.numeric(paste(datatot$SG_1))
  } else if (SIZE == "S2") {
    datatot$Response <- as.numeric(paste(datatot$SG_2))
  } else if (SIZE == "S3"){
    datatot$Response <- as.numeric(paste(datatot$SG_3))
  } else if (SIZE == "S4"){
    datatot$Response <- as.numeric(paste(datatot$SG_4))
  } else if (SIZE == "S5"){
    datatot$Response <- as.numeric(paste(datatot$SG_5))  
  } else if (SIZE == "S6"){
    datatot$Response <- as.numeric(paste(datatot$SG_6))
  } else if (SIZE == "S7"){
    datatot$Response <- as.numeric(paste(datatot$SG_7))
  } else if (SIZE == "S8"){
    datatot$Response <- as.numeric(paste(datatot$SG_8))
  } else if (SIZE == "S9"){
    datatot$Response <- as.numeric(paste(datatot$SG_9))
  } else if (SIZE == "S10"){
    datatot$Response <- as.numeric(paste(datatot$SG_10))
  } else if (SIZE == "S11"){
    datatot$Response <- as.numeric(paste(datatot$SG_11))
  } else if (SIZE == "S12"){
    datatot$Response <- as.numeric(paste(datatot$SG_12))
  } else if (SIZE == "S13"){
    datatot$Response <- as.numeric(paste(datatot$SG_13))
  } else if (SIZE =="S14"){
    datatot$Response <- as.numeric(paste(datatot$SG_14))
  }
}




# 3.4) Create equally time spaced intervals - VERY important for the AR1 process
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datatot$Year <- as.numeric(as.character(datatot$Year))
timeLevels <- as.vector(t(outer(min(datatot$Year):max(datatot$Year), 1:4, paste)))
datatot$YearQuarter <- factor(paste(datatot$Year, datatot$Quarter), levels=timeLevels)



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 4: Building grid for the study area
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grid for both commercial and survey data should be the same.
# It doesn't matter wheter to construct grid based on commercial or survey data, and one can arbitrary choose which dataset to use.
# Here we take the commercial data to do this.


# 4.1) Creating a dataframe with mean values of long and lat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comFULL$lon_mean <- rowMeans(comFULL[,c("lonStart", "lonEnd")])
comFULL$lat_mean <- rowMeans(comFULL[,c("latStart", "latEnd")])

df <- data.frame(lon=comFULL$lon_mean, lat=comFULL$lat_mean) #temporary df


# 4.2) Building the grid
#~~~~~~~~~~~~~~~~~~~~~~~~~
if(.Platform$OS.type == "windows") setwd("C:/Users/mruf/Documents/LGNB/Shapefiles")
grid <- gridConstruct2(df,km=5,scale=1.2)
gr <- gridFilter2(grid,df,icesSquare = T,connected=T) # filter out unnecessary spatial extensions
# plot(gr)



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 5: Discretize and associate hauls along grid cells 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.1) Setting a data frame containing the haul ID, and start and end long/lat of the haul
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- data.frame(sampleID=datatot$HLID, start_long=datatot$lonStart,
                  start_lat=datatot$latStart, end_long=datatot$lonEnd, end_lat=datatot$latEnd)

# Plotting trawl paths
#segments(dat$start_long,dat$start_lat, dat$end_long, dat$end_lat, col="red",lwd=1.2)



# 5.2) Define a matrix for the haul´s start and end position
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mStart <- matrix(c(dat$start_long,dat$start_lat), ncol=2)
mEnd <- matrix(c(dat$end_long,dat$end_lat), ncol=2)



# 5.3) Interpolate points at regular distance (default is 1km)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nbpts <- floor(distance(dat$start_long,dat$start_lat, dat$end_long, dat$end_lat) / 1) # one point every 1 km ~ reasonable for a 5x5 km grid
inter_pts <- gcIntermediate(mStart,mEnd, n=nbpts, addStartEnd=FALSE) 



# 5.4) Associate the discretized hauls to the grid ID 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note that the haul Id MUST be a factor, where each level is the frequency of a particular haul crossing a specific grid ID
tmp <- lapply(1:length(inter_pts), function(i) {
  print(i)
  x <- inter_pts[[i]]
  colnames(x) <- c("lon", "lat") #Needs to be the same names as those in inter_pts
  x <- as.data.frame(x)
  haul.id <- datatot$HLID[i] #Pick the specific haul id
  ind <- gridLocate(gr, x) #Locate the grid ID
  data.frame(haul.id=haul.id, ind=ind, rowID = i) 
})
tmp2 <- do.call("rbind", tmp)
tmp2$haulid <- factor(tmp2$haul.id)
tmp2$gf <- factor(tmp2$ind, 1:nrow(gr))
tmp2 <- tmp2[c("haulid","gf","rowID")] 



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 6: Defining the support areas 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# To be used later in association with the alpha-parameter.


# 6.1) For single support area
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here the haul positions of the aggregated time-series are considered to assigne
# a unique support area for the whole time-series.
# This is the MSA (model-single-alpha) case for the preferential sampling correction method.

datatot$split_area <- ifelse(as.character(datatot$Data)=="commercial", "commercial", "survey") #Takes the haul positions of the aggregated time-series 
datatot$split_area <- as.factor(datatot$split_area) #IMPORTANT - needs to be a factor!!
tmpOne <- tmp2; tmpOne$split <- datatot$split_area[tmp2$rowID]
SupportAreaMatrix    <- table(tmpOne$gf, tmpOne$split)
SupportAreaMatrix[]  <- SupportAreaMatrix>0
SupportAreaMatrix    <- ifelse(SupportAreaMatrix==0,FALSE,TRUE)
# levels(datatot$split_area); levels(datatot$Data) #Check


# 6.2) For multiple support areas
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Similar approach as above, but the amount of support areas will
# reflect the time-frame one is analyzing (note the difference between datatot$split_area and datatot$split_area2).
# This is the MMA (model-multiple-alphas) case for the preferential sampling correction method.

if(SUPPORT_AREA=="Several"){ # For multiple support areas
  
  datatot$split_area2 <- ifelse(as.character(datatot$Data)=="commercial", as.character(datatot$YearQuarter), "survey") #Takes the haul positions of the disaggregated time-series (e.g., monthly, quarterly,etc.)
  datatot$split_area2 <- as.factor(datatot$split_area2)
  tmpMulti <- tmp2; tmpMulti$split <- datatot$split_area2[tmp2$rowID]
  SupportAreaMatrix2 <- table(tmpMulti$gf, tmpMulti$split)
  SupportAreaMatrix2[] <- SupportAreaMatrix2>0
  SupportAreaMatrix2 <- ifelse(SupportAreaMatrix2==0,FALSE,TRUE)
}


# 6.3) Setting support areas based on chosen input
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(INCLUDE == "commercial" & SUPPORT_AREA == "Several"){
  SupportAreaMatrix2[] <- SupportAreaMatrix[,1]
  SupportAreaMatrix <- SupportAreaMatrix2 
} else if(INCLUDE == "both" & SUPPORT_AREA == "Several"){
  SupportAreaMatrix2[,1:(ncol(SupportAreaMatrix2)-1)] <- SupportAreaMatrix[,1]
  SupportAreaMatrix <- SupportAreaMatrix2 
} else if(INCLUDE == "commercial" & SUPPORT_AREA == "One"){
  SupportAreaMatrix <- SupportAreaMatrix
} else if (INCLUDE == "both" & SUPPORT_AREA =="One"){
  SupportAreaMatrix <- SupportAreaMatrix
}

#image(gr, SupportAreaMatrix[,1]) #To see the progress..


# 6.4) Map haulid to data.frame rows 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VERY IMPORTANT: haulid must match with the dataframe's row number (in increasing order)
rowid <- match(as.character(tmp2$haulid),as.character(datatot$HLID))
stopifnot(all(is.finite(rowid)))
rowid <- factor(rowid)



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 7: TMB processing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 7.1) Sparse matrices for GMRF: Q = Q0+delta*I
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Q0 <- -attr(gr,"pattern")
diag(Q0) <- 0
diag(Q0) <- -rowSums(Q0)
I <- .symDiagonal(nrow(Q0))



# 7.2) Compile LGNB model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(.Platform$OS.type == "windows") setwd("C:/Users/mruf/Documents/LGNB/src") #Set your own directory where the C++ is stored
compile("LGNB.cpp")
dyn.load(dynlib("LGNB"))



# 7.3) Prepare TMB data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TMB data are set in such way that it automatically
# recognizes the data-specific inputs. 



# 7.3.1) Linking R data to TMB data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- list(
  time = datatot$YearQuarter[rowid], # Temporal correlation 
  gf = tmp2$gf, # Satial correlation
  Q0 = Q0, # Spatial correlation
  I = I, # Spatial correlation
  Xpredict = matrix(0,0,0), # Covariate estimation matrix
  Apredict = factor(numeric(0)), # Covariate prediction matrix (FIXME: disabled in this current version)
  rowid = rowid ,
  response = datatot$Response,
  SupportAreaMatrix = SupportAreaMatrix,
  SupportAreaGroup = if(SUPPORT_AREA == "Several"){ # Links the support area matrix to TMB
    as.factor(datatot$split_area2)
  } else if(SUPPORT_AREA == "One"){
    as.factor(datatot$split_area)
  },
  Data=datatot$Data,
  h = mean(summary(as.polygons(gr))$side.length) # Optional; Can be used  to plot the spatial decorrelation as a function of distance
)



# 7.3.2) Optimize and minimize the objective function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_model <- function(data, model_struct=NULL, with_static_field=FALSE) {
  time_levels <- levels(data$time) # There must be at least one time levels; They MUST match between the two data types!
  grid_nlevels <- nlevels(data$gf)
  data$doPredict <- 0
  data$offset <- model_struct$offset 
  
  if(FALSE) { ## FIXME: prediction disabled
    ## Stuff for prediction: Dummy dataset that matches the space time grid:
    ## Xpredict: design matrix for prediction
    ## Apredict: Area factor for prediction
    ##DFpredict <- expand.grid(gf=levels(data$gf), time=levels(data$time))
    DFpredict <- expand.grid(gf=levels(data$gf), time=levels(data$time))
    ## FIXME: We should include depth and covariates here !
    ##        But that requires depth on the entire grid...
    Xpredict <- model.matrix(~time, data=DFpredict) ## <-- gear removed
    stopifnot( all( colnames(Xpredict) %in% colnames(data$X) ) ) ## Validity check
    tmp <- matrix(0, nrow(Xpredict), ncol(data$X))
    tmp[,match(colnames(Xpredict), colnames(data$X))] <- Xpredict
    data$Xpredict <- tmp
    data$Apredict <- factor(icesSquare(gr))
  }
  ## Perhaps we want measure all indices relative to a fixed reference square:
  ##   plot(cod, plot.response = FALSE)
  ## "42G1" seems appropriate
  ## data$refindex <- which(levels(data$Apredict) == "42G1")
  data$refindex <- 0 ## <-- Disable
  
  parameters <- list(
    eta_density = matrix(0,nrow(Q0),length(time_levels)),
    eta_nugget = numeric(0),
    logdelta = -4,       # Check values
    logscale = 0,         # Check values
    logsd_nugget = 0,    # Check values
    time_corr = 2,       # Check values
    beta = rep(0, ncol(data$X)),
    logphi = rep(0, nlevels(data$Data)),
    alpha = rep(0, nlevels(data$SupportAreaGroup))
  )
  parameters$eta_static <- rep(0, grid_nlevels * with_static_field )
  parameters$logdelta_static <- rep(0, 1 * with_static_field )
  parameters$logscale_static <- rep(0, 1 * with_static_field )
  
  ## Plugin model specification
  if(!is.null(model_struct)) {
    data$X                  <- model_struct$Xf
    data$beta_r_fac         <- model_struct$beta_r_fac
    parameters$beta         <- model_struct$beta
    parameters$beta_r       <- model_struct$beta_r
    parameters$beta_r_logsd <- model_struct$beta_r_logsd
  }
  
  ## Prior std dev on fixed effects (for robustness only)
  data$huge_sd <- 100
  
  map <- list()
  if(TRUE) map$logsd_nugget <- factor(NA)
  
  if(ALPHA == "No" & INCLUDE == "both"){
    map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
  } else if(ALPHA == "No" & INCLUDE == "commercial") {
    map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
  } else if(ALPHA == "No" & INCLUDE == "survey"){
    map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
  } else if(ALPHA == "Single I"){
    # map$alpha = factor(map$alpha)
  } else if(ALPHA == "Single II"){
    map$alpha <- factor(levels(data$SupportAreaGroup) != "survey") 
    map$alpha = factor(map$alpha)
  } else if(ALPHA == "Multi"){
    # map$alpha = factor(map$alpha)
  } 
  
  
  if(length(parameters$beta)>0 || length(parameters$beta)>0) {
    profile <- c("beta")
  } else {
    profile <- NULL
  }
  
  obj <- MakeADFun(data, parameters, random=c("eta_density","eta_nugget","eta_static","beta_r"),
                   profile = profile,
                   map=map, DLL="model")
  
  obj$env$tracepar <-TRUE
  
  .GlobalEnv$obj <- obj ## Copy for debugging
  runSymbolicAnalysis(obj)
  fit <- nlminb(obj$par,obj$fn,obj$gr)
  if(FALSE) { ## FIXME: disabled
    rep <- obj$report(obj$env$last.par.best)
    rownames(rep$logindex) <- levels(data$Apredict)
    colnames(rep$logindex) <- levels(data$time)
  }
  
  sdr <- sdreport(obj)
  s <- summary(sdr)
  est <- s[rownames(s) != "eta_density" & rownames(s) != "eta_nugget",]
  s <- summary(sdr,p.value = TRUE)
  s1 <- s[rownames(s) == "beta",]
  rownames(s1) <- head(colnames(data$X), nrow(s1)) # extracting names of the fixed effects
  s.fixed <- s1
  s1 <- s[rownames(s) == "beta_r",]
  rownames(s1) <- tail(colnames(data$X), nrow(s1)) # extracting names of random effects
  s.random <- s1
  
  s1
  s2 <- s[rownames(s) == "eta_density",]
  return(environment())
}




# 7.4) Model matrix
#~~~~~~~~~~~~~~~~~~~
buildModelMatrices <- function(fixed, random=NULL, ..., offset=NULL, data) {
  mm <- function(formula, data) {
    ##myna<-function(object,...){object[is.na(object)]<-0; object}
    mf <- model.frame(formula, data, na.action=na.pass)
    ans <- model.matrix(formula, mf)
    ans[is.na(ans)] <- 0
    ans
  }
  Xf <- mm(fixed, data=data)
  if(!is.null(random))
    Xr <- lapply(list(random, ...), mm, data=data)
  else
    Xr <- list(matrix(NA, nrow(Xf), 0))
  nr <- sapply(Xr, ncol)
  nf <- ncol(Xf)
  beta <- rep(0, nf)
  beta_r <- rep(0, sum(nr))
  beta_r_fac <- factor(rep(1:length(nr), nr))
  beta_r_logsd <- rep(0, nlevels(beta_r_fac))
  offset <- eval(offset, data)
  if(!is.null(offset)) stopifnot(is.numeric(offset)) else offset <- numeric(0)
  list(Xf=cbind(Xf, do.call("cbind", Xr)), nr=nr, nf=nf, beta=beta, beta_r=beta_r, beta_r_fac=beta_r_fac, beta_r_logsd=beta_r_logsd, offset=offset)
}



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 8: Fitting LGNB model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Model matrices were now built in such way that it allows to include both fixed and random effects;
# If random effects are included, then a separated formula needs to be specified, apart from the fixed effect formula;
# E.g.: buildModelMatrices (~ Fix_1 + Fix_2, ~ Random_1 -1);
# In this example, two variables are included as fixed effect, whereas one variable was considered as random;
 
# Obs.: When including random effect, the intercept needs to be substracted (-1) from it.

# Note that when a model contains one or more factor covariate (either as fixed or random effect)
# that is data-specific (e.g. metier), one needs to fit the model in such way that it runs with
# the whole formula structure (~ time, ~ metier - 1) for the data in which the factor is present
# (in this case commercial and combined model), and with a subsetted formula structure (~ time)
# for the data in which the factor does not appear (in this case the survey model)
# See examples below



## M1: Fixed effects: time-period; Random effects: vessel length & metier 
if (MODEL_FORMULA == "m2") {
  if (INCLUDE %in% c("commercial", "both")) {
    m1 <- buildModelMatrices(~ TimeYear + offset(log(HaulDur)), ~Metiers-1 + VE_LENcat-1, data=datatot)
  } else {
    m1 <- buildModelMatrices(~ TimeYear + offset(log(HaulDur)), data=datatot)
  }
  env1 <- fit_model(data, m1, with_static_field = F)
}


## M2: Fixed effects: time-period and sea bottom Depth; Random effects: vessel length & metier
if (MODEL_FORMULA == "m3") {
  if (INCLUDE %in% c("commercial", "both")) {
    m2 <- buildModelMatrices(~ timeyear2 + Depth + offset(log(HaulDur)), ~metiers-1 + VE_LENcat-1, data=datatot)
  } else {
    m2 <- buildModelMatrices(~ timeyear2 + Depth + offset(log(HaulDur)), data=datatot)
  }
  env2 <- fit_model(data, m2, with_static_field = F)
}


## M3: Fixed effects: time-period & Depth^2; Random effects: vessel length & metier
if (MODEL_FORMULA == "m4") {
  if (INCLUDE %in% c("commercial", "both")) {
    m3 <- buildModelMatrices(~ timeyear2 + poly(Depth,2) + offset(log(HaulDur)), ~metiers-1 + VE_LENcat-1, data=datatot)
  } else {
    m3 <- buildModelMatrices(~ timeyear2 + poly(Depth,2) + offset(log(HaulDur)), data=datatot)
  }
  env3 <- fit_model(data, m3, with_static_field = F)
}


## M5: Fixed effects: time levels & Sediment; Random effects: vessel length & metier as random effect
if (MODEL_FORMULA == "m5") {
  if (INCLUDE %in% c("commercial", "both")) {
    m5 <- buildModelMatrices(~ timeyear2 + Sediment + offset(log(HaulDur)), ~metiers-1 + VE_LENcat-1, data=datatot)
  } else {
    m5 <- buildModelMatrices(~ timeyear2 + Sediment + offset(log(HaulDur)), data=datatot)
  }
  env5 <- fit_model(data, m5, with_static_field = F)
}


## M6: Fixed effects: time levels, Depth & Sediment; Random effects: vessel length & metier as random effect
if (MODEL_FORMULA == "m6") {
  if (INCLUDE %in% c("commercial", "both")) {
    m6 <- buildModelMatrices(~ timeyear2 + Depth + Sediment + offset(log(HaulDur)), ~metiers-1 + VE_LENcat-1, data=datatot)
  } else {
    m6 <- buildModelMatrices(~ timeyear2 + Depth + Sediment + offset(log(HaulDur)), data=datatot)
  }
  env6 <- fit_model(data, m6, with_static_field = F)
}



## M7: Fixed effects: time levels, Depth^2 & Sediment; Random effects: vessel length & metier as random effect
if (MODEL_FORMULA == "m7") {
  if (INCLUDE %in% c("commercial", "both")) {
    m7 <- buildModelMatrices(~ timeyear2 + poly(Depth,2) + Sediment + offset(log(HaulDur)), ~metiers-1 + VE_LENcat-1, data=datatot)
  } else {
    m7 <- buildModelMatrices(~ timeyear2 + poly(Depth,2) + Sediment + offset(log(HaulDur)), data=datatot)
  }
  env7 <- fit_model(data, m7, with_static_field = F)
}



## M8: Fixed effects: time levels, Depth, Sediment & Depth:Sediment; Random effects: vessel length & metier as random effect
if (MODEL_FORMULA == "m8") {
  if (INCLUDE %in% c("commercial", "both")) {
    m8 <- buildModelMatrices(~ timeyear2 + Depth*Sediment + offset(log(HaulDur)), ~metiers-1 + VE_LENcat-1, data=datatot)
  } else {
    m8 <- buildModelMatrices(~ timeyear2 + Depth*Sediment + offset(log(HaulDur)), data=datatot)
  }
  env8 <- fit_model(data, m8, with_static_field = F)
}



#~~~~~~~~~~~~~
# Save results
#~~~~~~~~~~~~~
#rm(list=setdiff(ls(), ls(pattern="datatot|gr"))) #CHECK IF IT WORKS AFTER MODELS HAVE BEEN RUN.

rm(list=c("cohort_com","cohort_sur","dat","df","df_cohort","envpred_commercial","envpred_survey",
          "fd_data","fd_data2","fid_data","mStart","mEnd","colsel_com","colsel_sur","inter_pts","nbpts","WBS","tmp")) #Drop reduntand objects to decrease memory usage


#obj$env$L.created.by.newton <- NULL ## Trim off very large object before saving
#local({obj<-NULL}, get(grep("^env",ls(),value=T))) 

OUTFILE  <- paste0("results_WBScod", INCLUDE, ".RData")
save.image(file=OUTFILE)
