#=====================
# Some extra functions
#=====================


# These functions will be used along the LGNB model



####################################################################################################


# Function to read all packages in a single line
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mLoad <- function(...) {
  sapply(sapply(match.call(), as.character)[-1], require, character.only = TRUE)
}



# Function to extract cohorts and create new variable response variable based on that (Response)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Base function to ensure that both datasets have the same timelevels
extractCohortLevels <- function(d1, d2 = NULL, yearclass=2005, quarterclass=1) {
  firstYear1 <- min(as.numeric(levels(d1$Year)))
  firstYear2 <- min(as.numeric(levels(d2$Year)))
  firstYear <- min(firstYear1, firstYear2)
  lastYear1 <- max(as.numeric(levels(d1$Year)))
  lastYear2 <- max(as.numeric(levels(d2$Year)))
  lastYear <- max(lastYear1, lastYear2)
  numYears <- lastYear - firstYear + 1
  #ny <- 1:numYears # to consider only from age_1 group onwards
  ny <- (1:numYears)-1 # to include age_0 group
  age0 <- (paste("Age", ny, sep="_"))
  #year0 <- as.character(yearclass + ny - 1) # to consider only from age_1 group onwards
  year0 <- as.character(yearclass + ny) # to include age_0 group
  Quarter <- as.character(rep(1:4, length(age0)))
  Age <- rep(age0,each=4)
  Year <- rep(year0,each=4)
  df <- data.frame(Year,Age,Quarter,stringsAsFactors=FALSE)
  if (quarterclass > 1) {
    rem <- (1:(quarterclass-1))
    df <- df[-rem,]
  }
  df$levels <- local(paste("Yearclass_", yearclass, ":", Year,":",Quarter,":", Age, sep=""), df)
  df
}



# Base function to ensure that both datasets have the same timelevels, and which should be specific to each data type
extractCohortLevels2 <- function(d1, d2 = NULL, yearclass=2005, quarterclass=1) {
  firstYear1 <- min(as.numeric(levels(d1$Year)))
  firstYear2 <- min(as.numeric(levels(d2$Year)))
  firstYear <- min(firstYear1, firstYear2)
  lastYear1 <- max(as.numeric(levels(d1$Year)))
  lastYear2 <- max(as.numeric(levels(d2$Year)))
  lastYear <- max(lastYear1, lastYear2)
  numYears <- lastYear - firstYear + 1
  ny <- 1:numYears
  age0 <- (paste("age", ny, sep="_"))
  year0 <- as.character(yearclass + ny - 1)
  datType1 <- rep(levels(d1$Data),length(age0))
  datType2 <- rep(levels(d2$Data),length(age0))
  datType  <- rep(c(datType1,datType2), each=4)
  quarter <- as.character(rep(1:4, length(age0), len=8*length(age0)))
  age <- rep(age0,each=4,times=2)
  year <- rep(year0,each=4,times=2)
  df <- data.frame(year,age,quarter,datType,stringsAsFactors=FALSE)
  if (quarterclass > 1) {
    rem <- (1:(quarterclass-1))
    df <- df[-rem,]
  }
  df$levels <- local(paste("Yearclass_", yearclass, ":", year,":",quarter,":", age,"_", datType, sep=""), df)
  df
}


####################################################################################################


# Function to extract cohort on a yearly basis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_cohort_year <- function(d, yearclass=2005) {
  ##yearclass <- 2005 #Year when data starts
  n <- 1:nlevels(d$Year) #no. of Year levels in the dataset - column name should be the same as in the dataset
  #n<-1:15
  age <- paste("age", n, sep="_")
  year <- as.character(yearclass + n - 1)
  dnew <- NULL
  for(i in n) {
    tmp <- subset(d, Year==year[i])
    tmp$numberTotal_observer <- tmp[[age[i]]]
    dnew <- rbind(dnew, tmp)
  }
  levels(dnew$Year) <- paste("Yearclass_", yearclass, ":", year,":", age, sep="")
  dnew
}


#########################################################################################


# Function to extract cohort on a quarter basis within each year
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_cohort_quarter <- function(d, df) {
    n <- nrow(df)
    dnew <- NULL
    for(i in 1:n) {
        tmp <- subset(d, Year==df$Year[i] & Quarter==df$Quarter[i])
        if (nrow(tmp) > 0) {
          tmp$Response <- tmp[[df$Age[i]]]
          tmp$year <- df$levels[i]
          dnew <- rbind(dnew, tmp)
        }
    }
    dnew$Cohort <- factor(dnew$year, levels=df$levels)
    dnew$year <- NULL
    dnew
}

##df <- extractCohortLevels(d, d)
##ec <- extract_cohort_quarter(d, df)

####################################################################################################


# Function for binding commercial and survey datasets into one single dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mybind <- function(x, y) {
  i <- intersect(names(x), names(y))
  ix <- setdiff(names(x), names(y))
  iy <- setdiff(names(y), names(x))
  z <- rbind(x[i], y[i])
  xna <- data.frame(lapply(iy, function(...)rep(NA, nrow(x))))
  names(xna) <- iy
  yna <- data.frame(lapply(ix, function(...)rep(NA, nrow(y))))
  names(yna) <- ix
  zx <- rbind(x[ix], yna)
  zy <- rbind(xna, y[iy])
  cbind(z, zx, zy)
}


####################################################################################################

# Constructing grid based on existing shapefile
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Use this function if using the contracted versions of the DK shapefile;
gridConstruct2 <- function(d, km,scale,nearestObs){
  ## Grid
  gr <- gridConstruct::gridConstruct(d=d,km=km,scale=scale,nearestObs=nearestObs,filter=FALSE)
  ## Read shape file data
  shape <- readOGR(".", "CDK2_cutted") #Importing contracted DK shapefile with width=-0.010
  proj4string(shape) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  map <- spTransform(shape,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  gr2 <- as.data.frame(gr)
  coordinates(gr2) <- ~lon + lat
  proj4string(gr2) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  xtra <- over(gr2, map)
  gr3 <- gr[is.na(xtra$id),]
  cc <- connectedComponents(gr3)
  gr3 <- gr3[cc[[which.max(sapply(cc,length))]],]
  while(any(rowSums(attr(gr3,"pattern")) <=2 ) ){
    gr3 <- gr3[rowSums(attr(gr3,"pattern"))>2,]
  }
  ##plot(map,add=TRUE)
  f <- function(i){
    x <- map@polygons[[i]]
    do.call("rbind",lapply(x@Polygons,function(x)rbind(x@coords,NA)))
  }
  qw <- do.call("rbind",lapply(1:length(map@polygons),f))
  plot(gr3)
  #points(qw[,1],qw[,2],type="l")
  polygon(qw[,1],qw[,2],col="grey70")
  ##gr2 <- spTransform(gr2,CRS(proj4))
  ##xtra <- over(gr2, shape)
  ##image(gr,xtra$Lower_,map=NULL)
  #gr <<- gr3
  map.pol <<- qw
  return(gr3)
  
}




# Function to filterout unwanted grid cells
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is a function applied after having constructed de grid.
# It is used to reduce the ammount of gridcells, and hence keep the grid only for the area of interest.
# In Kasper's original function, one can fitler out gridcells based on the ICES squares. To do so, ICESsquares are assigned
# based on the data points. Because we want to extend the grid a bit further to the actual ICESsquares that are based on the datapoints,
# I did a quick and dirty manual adaptation, where additional ICESsquares were added.

Kattegat <- c("44G0","44G1","43G0","43G1","43G2","42G0","42G1","42G2","41G0","41G1","41G2") #subdivision 21
Belt <- c("40F9","40G0","40G1","39F9","39G0","39G1","38F9","38G0","38G1","37G0","37G1") #subdivision 22
Sound <- c("40G2") #subdivision 23
Baltic_Sea <- c("39G2","39G3","39G4","38G2","38G3","38G4","37G2","37G3","37G4") #subdivision 24
#Extrasquares <- c("36G4","36G3","40G4","37G5","38G5","39G5","40G5","41G5","41G4","45G1","45G0")
Extrasquares <- c("36G4","36G3","40G4","37G5","38G5","39G5","40G5","41G5","41G4","40G3")

studyarea <- c(Kattegat,Belt,Sound,Baltic_Sea,Extrasquares)


##' @param icesSquare Remove grid points outside ICES squares in the data?
##' @param nearestObs Remove grid points with closest data point greater than \code{nearestObs}.
##' @param wet Remove grid points on land.
##' @param wetEdges Alternative: Keep edges passing through water. 
##' @param ordertol Require at least \code{ordertol} neighbors to every grid point.
##' @param connected Keep only largest connected component.
##' @return Filtered grid object
##' @rdname grid
gridFilter2 <- function(grid,data,
                       icesSquare=FALSE,
                       nearestObs=Inf,
                       wet=FALSE,
                       connected=FALSE,
                       ordertol=0,
                       wetEdges=FALSE,
                       ...){
  keep <- rep(TRUE,nrow(grid))
  ## icesSquare filtering
  if(icesSquare){
    require(DATRAS)
    #sq <- unique(DATRAS::icesSquare(data)) #Kasper's original code
    sq <- studyarea #My own contribution
    ##pol <- icesSquare2coord(sq,format="polygons")
    ##keep <- keep & !is.na(cut(grid,pol))
    sqg <- DATRAS::icesSquare(grid)
    keep <- keep & (sqg %in% sq)
  }
  ## nearest observation filtering
  if(is.finite(nearestObs)){
    ## i <- gridLocate(data,grid)
    ## d <- dist.km(grid,data[i,],outer=FALSE)
    ## keep <- keep & (d<=nearestObs)
    i <- gridLocate(data,grid[keep,])
    d <- dist.km(grid[keep,],data[i,],outer=FALSE)
    keep[keep] <- (d<=nearestObs)
  }
  ## map filtering
  if(wet){
    ##keep <- keep & is.na(cut(grid))
    keep[keep] <- is.na(cut(grid[keep,]))
  }
  if(wetEdges){
    ## Reduce pattern matrix first
    Q <- attr(grid,"pattern")
    Q[!keep,] <- 0
    Q[,!keep] <- 0
    attr(grid,"pattern") <- Q
    ## Now "wetEdgeFilter" has less work
    grid <- wetEdgeFilter(grid,...)
  }
  grid <- grid[keep,]
  ## Remove gridpoints with too few neighbors
  if(ordertol>0){
    rem <- function(i0) {
      p0 <- attr(grid, "pattern")[i0, i0]
      cs <- colSums(p0) - 1
      keep <- cs > ordertol
      i0 <- i0[keep]
      if (any(!keep)) 
        i0 <- rem(i0)
      i0
    }
    i0 <- rem(seq(length=nrow(grid)))
    grid <- grid[i0,]
  }
  ## Keep largest connected component
  if(connected){
    c <- connectedComponents(grid)
    if(length(c)>1){
      c <- c[[which.max(sapply(c,length))]]
      grid <- grid[c,]
    }
  }
  grid
}


####################################################################################################


# Function to calculate distance between two points
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
distance <- function (lon, lat, lonRef, latRef)
{
  pd <- pi/180
  a1 <- sin(((latRef - lat) * pd)/2)
  a2 <- cos(lat * pd)
  a3 <- cos(latRef * pd)
  a4 <- sin(((lonRef - lon) * pd)/2)
  a <- a1 * a1 + a2 * a3 * a4 * a4
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  return(6371 * c)
}


####################################################################################################


# Function to locate points on nearest grid vertex
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gridLocateBin2 <- function (grid, points) 
{
  recursiveLocate <- function(x, y) {
    if (nrow(x) == 0 | nrow(y) == 0) 
      return(NULL)
    if (nrow(x) == 1) 
      return(rep(x$key, nrow(y)))
    x[1:2] <- x[2:1]
    y[1:2] <- y[2:1]
    med <- median(x[, 1])
    i <- factor(x[, 1] < med, levels = c(TRUE, FALSE))
    spx <- split(x, i)
    j <- factor(y[, 1] < med, levels = c(TRUE, FALSE))
    spy <- split(y, j)
    rm(x, y)
    ans <- mapply(recursiveLocate, spx, spy, SIMPLIFY = FALSE)
    if (any(sapply(ans, length) == 0)) 
      return(unlist(ans))
    else return(unsplit(ans, j))
  }
  grid <- as.data.frame(grid[c("lon", "lat")])
  points <- as.data.frame(points[c("lon", "lat")])
  grid$key <- 1:nrow(grid)
  recursiveLocate(grid, points)
}


# Modified gridFactor function from Kasper's original function
gridFac <- function (data, grid, ...) {
            #data <- data[c("lon", "lat")] #Kasper's code
            data <- data[c("lonMean", "latMean")]
            #ch <- paste(data$lon, data$lat)
            ch <- paste(data$lonMean, data$latMean)
            nd <- which(!duplicated(ch))
            fac <- unclass(factor(ch, levels = ch[nd]))
            i <- gridLocate(grid, data[nd, ])
            fac <- factor(i[fac], 1:nrow(grid))
            levelAttrib(fac)$grid <- grid
            class(fac) <- c("gridFactor", class(fac))
            return(fac)
}


####################################################################################################


# Guess time effects (needed by index calculation)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
guess_time_effects_in_beta <- function(data, time_levels) {
  cn <- colnames(data$X)
  cn <- gsub(" ",":",time_levels)
  tl <- gsub(" ",":",time_levels)
  li <- lapply(tl, function(nm) grep(nm, cn) )
  if (any( sapply(li, length) > 1 )) {
    print(lapply(li, function(i)cn[i]))
    warning("Time effect is not unique. Index calc will be omitted")
    return (integer(0))
  }
  found <- sapply(li, function(x)x[1])
  found[is.na(found)] <- 0
  ## NOTE: Output length = length(time_levels)
  ## NOTE: NAs coded as "-1" ===> SKIP index calc or crash !!!
  as.integer(found) - 1L
}



####################################################################################################


# Whenever you are running a model (e.g., GLM within R basic functions), all your information is
# stored in matrices. So if you, for example, specify a model as response ~ intercept + covariate A + covariate B, 
# these information will be stored in a so-called design matrix, where you would have typically two
# vectors (one for the response, another for the intercept) and one covariate matrix.
# To do so, you could usually rely on the model.matrix R basic function. However, this function
# provides some serious limitation when some factor levels of a particular covariate present NAs. 
# The function normally removes empty factor levels, but when dealing with both datasets together 
# (combined model), keeping all factor levels is very important to fit the model. For instance, 
# the métier effect is only present in the commercial data. Thus, when binding both datasets, 
# the final dataframe (herein called as datatot) will have empty factor levels for the survey 
# data whenever the métier column is regarded.

# Thus, to circumvent this problem, we modified the model.matrix R function within our
# own created function (buildModelMatrices). The function provides a user-friendly interface, 
# where you can specify your model in a similar way as done within a GLM(M)/GAM(M) formulation. ¨
# The nice aspect of this function is that the fixed effect(s), random effect(s) and offset(s)
# can be defined as a standard model formula (see section 8 for more details).

# Here it should be noted that the function also ensures that the model matrix recognizes
# immediately which data type is used as an input (commercial, survey or commercial+survey),
# and thereby adapting itself according to the data considered.



# Build Matrices for the model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
buildModelMatrices <- function(fixed, random=NULL, ..., offset=NULL, static=NULL, data) {
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
  if(!is.null(static))
    Xs <- mm(static, data=attr(data,"static"))
  else
    #Xs <- matrix(NA, 0, 0) # Kasper's version; works on my local machine but not in hpc
    Xs <- matrix(NA, nrow(gr), 0) #works on both local machine and hpc
  nr <- sapply(Xr, ncol)
  nf <- ncol(Xf)
  beta <- rep(0, nf)
  beta_r <- rep(0, sum(nr))
  beta_r_fac <- factor(rep(1:length(nr), nr))
  beta_r_logsd <- rep(0, nlevels(beta_r_fac))
  offset <- eval(offset, data)
  ns <- ncol(Xs)
  beta_s <- rep(0, ns)
  if(!is.null(offset)) stopifnot(is.numeric(offset)) else offset <- numeric(0)
  list(Xf=cbind(Xf, do.call("cbind", Xr)), Xs=Xs, nr=nr, nf=nf, ns=ns, beta=beta, beta_s=beta_s, beta_r=beta_r, beta_r_fac=beta_r_fac, beta_r_logsd=beta_r_logsd, offset=offset)
}
