##########################################################################################
#                                                                                        #
##              LGNB model: A practical statistical framework to combine                ##
##               commercial & survey data to model the spatio-temporal                  ##
##                          dynamics of marine harvested species                        ##
##                                    (Rufener et al.)                                  ##
#                                                                                        #
##########################################################################################

# last update: April 2021


# To keep track of all the processes in a logical and concise manner,
# the following script was divided into 8 different sections:

## Section 1: Default inputs
## Section 2: Load data files & packages
## Section 3: Binding data files 
## Section 4: Grid building
## Section 5: Discretize commercial hauls 
## Section 6: Defining the preferential sampling, if present 
## Section 7: TMB processing 
## Section 8: Fitting the LGNB model 





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 1: Default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
STOCK <- c("WBS","KAT")[1] # Specify to which stock the data should be subsetted (Western Baltic Sea or Kattegat)
DATA  <- c("commercial", "survey", "both") [3] # Specify the desired input data for the model; default is commercial data-
RESPONSE <- c("SizeGroup","AgeGroup","Cohort")[2] #Choose whether the model is applied for each SizeGroup, AgeGroup, or on a Cohort basis (when Nage).Default is set to SizeGroup.
PS <- c("No","One","Two")[1] #Define how the sampling nature should be accounted for
TIME <- c("YearMonth","YearQuarter","Year")[2] #Define the temporal resolution; default is set on a monthly basis


# @ PS = "No" -> For both datasets (commercial, survey), no preferential sampling is accounted for.
# @ PS = "One" -> Preferential sampling is accounted only for the commercial data.
# @ PS = "Two" -> Preferential sampling is accounted in both survey and commercial data; applies only to the integrated model (DATA="both")




# Specify the model structure; default model is m2for either size groups, age groups or cohort
# Default size group is 5 (S5), age group is 3 (A3) and cohort from 2005.
if(RESPONSE == "SizeGroup"){
  MODEL_CONFIG <- "m1_S5" #Default model and SizeGroup 5
}else if(RESPONSE == "AgeGroup"){
  MODEL_CONFIG <- "m1_A3" #Default model and AgeGroup 3
} else if(RESPONSE == "Cohort")
  MODEL_CONFIG <- "m1_2011" #Default model when applied on cohort-basis (2005 cohort) 


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
stopifnot(DATA %in% c("commercial", "survey", "both"))



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 2: Load data files & R packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2.1) Load helper functions and R libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("C:/Users/mruf/Documents/LGNB/R/utilities.R")

#devtools::install_github("kaskr/gridConstruct",subdir="gridConstruct") # To install the gridConstruct package
mLoad(raster,rgeos,maptools,maps,data.table,dplyr,TMB,sp,
      DATRAS,gridConstruct,rgdal,geosphere,devtools,plyr,fields)



# 2.2.1) Load commercial fisheries data (fishery-depedent)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comFULL <- readRDS("C:/Users/mruf/Documents/LGNB/Data/commercial.rds")
comFULL$stock <- ifelse(comFULL$Area=="21","KAT","WBS")


# Subset data for a particular stock
commercial <- if(RESPONSE == "Cohort"){ #Subsetting when response is on a cohort level also subsets the data automatically for the time-period related to the yearclass being analyzed
  t <- as.numeric(as.character(comFULL$Year))
  tmax <- as.factor(max(t))
  
  yearclass <- as.numeric(as.character(YEARCLASS))
  tperiod <- as.factor(c(yearclass:as.numeric(as.character(tmax))))
  
  subset(comFULL, stock == STOCK & Year %in% tperiod) #Setting stock based on ICES area (KAT=21, WBS=22-24)
} else{
  subset(comFULL, stock == STOCK)
}



# Drop unused factor levels
commercial[,c("Month","Year","Quarter","Area","metiers","Data","HAULID","VESSELID","TimeYear")] <- lapply(commercial[,c("Month","Year","Quarter","Area","metiers","Data","HAULID","VESSELID","TimeYear")], factor)


# Set haulduration to 1 (ensures that no offset is accounted in the commercial data)
commercial$haulduration_hours <- 1 #Note: log(1)=0!


# Create a fake "Age_0" column - for consistency when applying the model on a cohort-basis
commercial$Age_0 <- NA



# 2.2.2) Load scientific survey data (Fishery~indepdendent)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
survey <- readRDS("C:/Users/mruf/Documents/LGNB/Data/survey.rds")
survey$stock <- ifelse(survey$Area=="21","KAT","WBS") #Setting stock based on ICES area (KAT=21, WBS=22-24)



# Subset data for a particular time frame
survey <- if(RESPONSE == "Cohort"){
  t <- as.numeric(as.character(survey$Year))
  tmax <- max(t)
  
  yearclass <- as.numeric(as.character(YEARCLASS))
  tperiod <- as.factor(c(yearclass:tmax))
  
  subset(survey, stock == STOCK & Year %in% tperiod) #Setting stock based on ICES area (KAT=21, WBS=22-24)
} else{
  
  subset(survey, stock == STOCK)
}



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
  
  cohort_com <- transform(cohort_com, HLID=HAULID, HaulDur=haulduration_hours)
  cohort_sur <- transform(cohort_sur, latStart=lat, lonStart=lon, latEnd=lat, lonEnd=lon, HLID=Haul_ID)
  
  
  # 3.1.2) AgeGroup-basis
  #~~~~~~~~~~~~~~~~~~~~~~~~
} else if (RESPONSE == "AgeGroup"){
  age_com   <- transform(commercial, HLID=HAULID, HaulDur=haulduration_hours)
  age_sur   <- transform(survey, latStart=lat, lonStart=lon, latEnd=lat, lonEnd=lon,HLID=Haul_ID)
  
  
  # 3.1.3) SizeGroup-basis
  #~~~~~~~~~~~~~~~~~~~~~~~~
} else if(RESPONSE == "SizeGroup"){
  size_com <- transform(commercial, HLID=HAULID, HaulDur=haulduration_hours)
  size_sur <- transform(survey, latStart=lat, lonStart=lon, latEnd=lat, lonEnd=lon,HLID=Haul_ID)
  
}



# 3.2) Bind both datasets into a single data frame
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(RESPONSE == "Cohort" & DATA =="both"){
  datatot <- mybind(cohort_com, cohort_sur)
} else if(RESPONSE == "Cohort" & DATA == "survey"){
  datatot <- cohort_sur
} else if(RESPONSE == "Cohort" & DATA == "commercial"){
  datatot <- cohort_com
} else if (RESPONSE == "AgeGroup" & DATA == "both"){
  datatot <- mybind(age_com, age_sur)
} else if (RESPONSE == "AgeGroup" & DATA == "survey") {
  datatot <- age_sur
} else if (RESPONSE == "AgeGroup" & DATA == "commercial"){
  datatot <- age_com
} else if (RESPONSE == "SizeGroup" & DATA == "both"){
  datatot <- mybind(size_com, size_sur)
} else if (RESPONSE == "SizeGroup" & DATA == "survey"){
  datatot <- size_sur
} else if(RESPONSE == "SizeGroup" & DATA == "commercial")
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

datatot$Year <- as.numeric(as.character(datatot$Year)) #Set as numeric to run the lines below


if(TIME=="YearMonth"){
  
  timeLevels <- as.vector(t(outer(min(datatot$Year):max(datatot$Year), 1:12, paste)))
  datatot$YearMonth <- factor(paste(datatot$Year, datatot$Month), levels=timeLevels)
  
} else if(TIME=="YearQuarter"){
  timeLevels <- as.vector(t(outer(min(datatot$Year):max(datatot$Year), 1:4, paste)))
  datatot$YearQuarter <- factor(paste(datatot$Year, datatot$Quarter), levels=timeLevels)
  
} 


datatot$Year <- as.factor(datatot$Year) #Set back to factor





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
grid <- gridConstruct2(df,km=20,scale=1.2)
gr <- gridFilter2(grid,df,icesSquare = T,connected=T) # filter out unnecessary spatial extensions
# shape <- readOGR(".", "CDK2_cutted"); plot(gr); plot(shape,add=T,col="grey70")



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
nbpts <- floor(distance(dat$start_long,dat$start_lat, dat$end_long, dat$end_lat) / 1) #one point every 1 km
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 6: Defining the preferential sampling 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If a preferential sampling behaviour in the commercial and/or survey data
# is accounted for, we have to define the so-called sampling support area.


# 6.1) For single support area
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here the haul positions of the aggregated time-series are considered to assigne
# a unique support area for the whole time-series.

datatot$split_area <- ifelse(as.character(datatot$Data)=="commercial", "commercial", "survey") #Takes the haul positions of the aggregated time-series 
datatot$split_area <- as.factor(datatot$split_area) #IMPORTANT - needs to be a factor!!
tmpOne <- tmp2; tmpOne$split <- datatot$split_area[tmp2$rowID]
SupportAreaMatrix    <- table(tmpOne$gf, tmpOne$split)
SupportAreaMatrix[]  <- SupportAreaMatrix>0
SupportAreaMatrix    <- ifelse(SupportAreaMatrix==0,FALSE,TRUE)
# levels(datatot$split_area); levels(datatot$Data) #Check



# 6.2) Setting support areas based on chosen input
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if(DATA == "commercial" || DATA == "both"){
#   SupportAreaMatrix <- SupportAreaMatrix
# }
#image(gr, SupportAreaMatrix[,1]) #To see the progress..


# 6.3) Map haulid to data.frame rows 
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
#Sys.setenv(PATH="%PATH%;C:/Rtools/gcc-4.6.3/bin;c:/Rtools/bin") #Run only when on windows
compile("LGNB.cpp")
dyn.load(dynlib("LGNB"))



# 7.3) Prepare TMB data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TMB data are set in such way that it automatically
# recognizes the data-specific inputs. 



# 7.3.1) Linking R data to TMB data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- list(
  time = if(TIME=="Year") {datatot$Year[rowid]} #Yearly resolution
  else if(TIME=="YearQuarter") {datatot$YearQuarter[rowid]} #Quarterly resolution
  else if(TIME=="YearMonth"){datatot$YearMonth[rowid]}, #Monthly resolution
  gf = tmp2$gf,              
  rowid = rowid,
  response = datatot$Response,
  Q0 = Q0,
  I = I,
  Xpredict = matrix(0,0,0),
  Apredict = factor(numeric(0)),
  SupportAreaMatrix = SupportAreaMatrix,
  SupportAreaGroup = as.factor(datatot$split_area),
  Data=datatot$Data,
  h = mean(summary(as.polygons(gr))$side.length) #To be used later to plot the spatial decorrelation as a function of distance
)


# 7.3.2) Optimize and minimize the objective function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Fitting the model
#~~~~~~~~~~~~~~~~~~~~

## 'profile' speeds up fitting but may be inconvenient
fit_model <- function(data, model_struct=NULL, with_static_field=FALSE, profile=TRUE) {
  time_levels <- levels(data$time) # There must be at least one time levels; They MUST match between the two data types!
  grid_nlevels <- nlevels(data$gf)
  data$doPredict <- 0
  data$offset <- model_struct$offset 
  
  if(FALSE) { ## FIXME: prediction disabled
    ## Stuff for prediction: Dummy dataset that matches the space time grid:
    ## Xpredict: design matrix for prediction
    ## Apredict: Area factor for prediction
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
  data$refindex <- 0 ## <-- Disable
  
  parameters <- list(
    eta_density = matrix(0,nrow(Q0),length(time_levels)),
    eta_nugget = numeric(0),
    logdelta = -4,       # Check values
    logscale = 0,         # Check values
    logsd_nugget = 0,    # Check values
    time_corr = 2,       # Check values
    beta = rep(0, ncol(data$X)),
    logphi = rep(0, nlevels(data$Data) ),
    alpha = rep(0, nlevels(data$SupportAreaGroup))
  )
  parameters$eta_static <- rep(0, grid_nlevels * with_static_field )
  parameters$logdelta_static <- rep(0, 1 * with_static_field )
  parameters$logscale_static <- rep(0, 1 * with_static_field )
  
  ## Plugin model specification
  if(!is.null(model_struct)) {
    data$Xs                 <- model_struct$Xs
    data$X                  <- model_struct$Xf
    data$beta_r_fac         <- model_struct$beta_r_fac
    parameters$beta         <- model_struct$beta
    parameters$beta_r       <- model_struct$beta_r
    parameters$beta_s       <- model_struct$beta_s
    parameters$beta_r_logsd <- model_struct$beta_r_logsd
    data$which_beta_time    <- guess_time_effects_in_beta(data, time_levels)
  }
  
  ## Prior std dev on fixed effects (for robustness only)
  data$huge_sd <- 100
  
  map <- list()
  if(TRUE) map$logsd_nugget <- factor(NA)
  
  if(PS == "No" & DATA == "both"){
    map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
  } else if(PS == "No" & DATA == "commercial") {
    map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
  } else if(PS == "No" & DATA == "survey"){
    map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
  } else if(PS == "One"){
    lev <- levels(data$SupportAreaGroup)
    lev[lev=="survey"] <- NA
    map$alpha <- factor(lev)
  } else if(PS == "Two"){
    lev <- levels(data$SupportAreaGroup)
    #lev[lev=="survey"] <- NA
    map$alpha <- factor(lev)
    # map$alpha = factor(map$alpha)
  } 
  
  
  if(profile && (length(parameters$beta)>0 || length(parameters$beta)>0) ) {
    profile <- c("beta")
  } else {
    profile <- NULL
  }
  
  obj <- MakeADFun(data, parameters, random=c("eta_density","eta_nugget","eta_static","beta_r"),
                   profile = profile,
                   map=map, DLL="LGNB")
  
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





# 7.3.3) Include environmental covariates for prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NOTE: Ignore this step if you want to predict the abundance only 
# as a function of the spatio-temporal correlation paramaters


source("C:/Users/mruf/Documents/LGNB/R/Extract_Covariates_Prediction.R")
covariates <- readRDS("~/pred_covariates.rds")

attr(datatot, "static") <- covariates




#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Section 8: Fitting LGNB model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The model matrices were now built in such way that it allows to include both fixed and random effects;
# If random effects are included, then a separated formula needs to be specified for each random effect,apart from the fixed effect formula;
# E.g.: buildModelMatrices (~ Fix_A + Fix_B, ~ Random_C -1, ~ Random_D -1);
# In this example, two variables (A and B) are included as fixed effect, and two as random (C and D);

# Note that when a model contains one or more factor covariate (either as fixed or random effect)
# that is data-specific (e.g. metier), one needs to fit the model in such way that it runs with
# the whole formula structure (~ time, ~ metier - 1) for the data in which the factor is present
# (in this case commercial and combined model), and with a subsetted formula structure (~ time)
# for the data in which the factor does not appear (in this case the survey model)

# Note, also, that the buildMoldeMatrices function allows you to specify
# two ways of including the fixed effect of the covariates (if present).
# The first way includes the covariates in the observation process, while the
# second way does it in the latent process. 
# Including covariates in the observation process allows to account for factors that are
# thought to affect the catchability, whereas in the second case it accounts for factors
# that are believed to affect the underlying abundance density field.


# But first, we need to revers the level order of the data, such that
# the survey data comes always first. THIS is very important, as the survey
# data will always be treated as the reference level; hence, everything is being
# parametrized upon the survey data


## Reverse levels (to re-parameterize)
data$Data <- factor(data$Data, rev(levels(data$Data))) #Make sure that survey data comes always first!
datatot$Data <- factor(datatot$Data, rev(levels(datatot$Data))) #Make sure that survey data comes always first!

# See some examples below:


#~~~~~~~
## M1 ##
#~~~~~~~

# Fixed effects on observation process: time-period (in both data sources), metiers (commercial data), Ship (survey data)
# Fixed effects on latent process: None
# Random effects: vessel ID (whenever commercial data is considered)
# Offset: log(haul duration)

if (MODEL_FORMULA == "m1") {
  if(DATA == "commercial"){
    if(TIME == "YearMonth"){
      
      m1 <- buildModelMatrices(fixed = ~  -1 + YearMonth + metiers, random = ~ -1 + VESSELID, offset = quote(log(HaulDur)), data=datatot)
    } else if(TIME == "YearQuarter"){
      m1 <- buildModelMatrices(fixed = ~  -1 + YearQuarter + metiers, random = ~ -1 + VESSELID, offset = quote(log(HaulDur)), data=datatot)
    } else if(TIME == "Year"){
      m1 <- buildModelMatrices(fixed = ~  -1 + Year + metiers, random = ~ -1 + VESSELID, offset = quote(log(HaulDur)), data=datatot)
    }
    
  } else if(DATA == "survey"){
    if(TIME == "YearMonth"){
      m1 <- buildModelMatrices(fixed = ~ -1 + YearMonth + Ship, offset = quote(log(HaulDur)), data=datatot)
    } else if(TIME == "YearQuarter"){
      m1 <- buildModelMatrices(fixed = ~ -1 + YearQuarter + Ship, offset = quote(log(HaulDur)), data=datatot)
    } else if(TIME == "Year"){
      m1 <- buildModelMatrices(fixed = ~ -1 + Year + Ship, offset = quote(log(HaulDur)), data=datatot)
    }
    
  } else {
    if(TIME == "YearMonth"){
      m1 <- buildModelMatrices(fixed = ~ -1 + YearMonth + Ship + Data + metiers, random = ~ -1 + VESSELID, offset= quote(log(HaulDur)), data=datatot)
    } else if (TIME == "YearQuarter"){
      m1 <- buildModelMatrices(fixed = ~ -1 + YearQuarter + Ship + Data + metiers, random = ~ -1 + VESSELID, offset= quote(log(HaulDur)), data=datatot)
    } else if(TIME == "Year"){
      m1 <- buildModelMatrices(fixed = ~ -1 + Year + Ship + Data + metiers, random = ~ -1 + VESSELID, offset= quote(log(HaulDur)), data=datatot)
    }
  }
  system.time( env1 <- fit_model(data, m1, with_static_field = F, profile=F) )
}



#~~~~~~~
## M2 ##
#~~~~~~~

# Fixed effects on observation process: time-period (in both data sources), metiers (commercial data), Ship (survey data), depth
# Fixed effects on latent process: None
# Random effects: vessel ID (whenever commercial data is considered)
# Offset: log(haul duration)
if (MODEL_FORMULA == "m2") {
  if(DATA == "commercial"){
    m2 <- buildModelMatrices(~  -1 + YearQuarter + metiers + Depth, random = ~ -1 + VESSELID, offset= quote(log(HaulDur)), data=datatot)
  } else if(DATA == "survey"){
    m2 <- buildModelMatrices(~ -1 + YearQuarter + Ship + Depth, offset = quote(log(HaulDur)), data=datatot)
  } else {
    m2 <- buildModelMatrices(~ -1 + YearQuarter + Ship + Data + metiers + Depth, random = ~ -1 + VESSELID, offset= quote(log(HaulDur)), data=datatot)
  }
  system.time( env2 <- fit_model(data, m2, with_static_field = F, profile=F) )
}




#~~~~~~~
## M3 ##
#~~~~~~~

# Fixed effects on observation process: time-period (in both data sources), metiers (commercial data), Ship (survey data)
# Fixed effects on latent process: Depth
# Random effects: vessel ID (whenever commercial data is considered)
# Offset: log(haul duration)
if (MODEL_FORMULA == "m2") {
  if(DATA == "commercial"){
    m3 <- buildModelMatrices(~  -1 + YearQuarter + metiers + Depth, random = ~ -1 + VESSELID, offset = quote(log(HaulDur)), data = datatot)
  } else if(DATA == "survey"){
    m3 <- buildModelMatrices(~ -1 + YearQuarter + Ship + Depth, offset = quote(log(HaulDur)), data = datatot)
  } else {
    m3 <- buildModelMatrices(~ -1 + YearQuarter + Ship + Data + metiers + Depth, static = ~ -1 + depth, ~ random = -1 + VESSELID, offset = quote(log(HaulDur)), data = datatot)
  }
  system.time( env3 <- fit_model(data, m3, with_static_field = F, profile=F) )
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Extract the results
#~~~~~~~~~~~~~~~~~~~~~~


env1$sdr #get TMB sdreport
env1$s.fixed #get fixed-effect values
env1$s.random

# Replace env2 by the appropriate model ID that was run


# Plotting
pl <- as.list(env1$sdr, "Estimate") #for the estimated abundances
plsd <- as.list(env1$sdr, "Std. Error") #for the SE of the estimated abundances


nr <- nlevels(factor(datatot$Year))
nc <- 4

par(mfrow=c(nr,nc))



# Plot estimated abundance fields
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in 1:ncol(pl$eta_density)){
  image(gr, concTransform(pl$eta_density[,i]),col=tim.colors(99))
  title(levels(datatot$YearQuarter)[i])
}



# For the SE of the estimated abundance field
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in 1:ncol(pl$eta_density)){
  image(gr, concTransform(plsd$eta_density[,i]),col=tim.colors(99))
  title(levels(datatot$YearQuarter)[i])
}

