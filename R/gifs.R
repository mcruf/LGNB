
##################################################################
#   Generate GIS from the model output - maps based on ggplot    #
##################################################################



library(ggplot2)
library(maptools)
library(gridConstruct)
library(raster)
library(ggsn)
library(mapdata)
library(fields)
library(rgdal)
library(pals)
library(tidyr)
library(gganimate)
library(magick)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directoy and load the data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/PhD/Manuscript_01/NEW submission/Results/WBScod")



#~~~~~~~~~~~~~~~~~~
# Survey data
#~~~~~~~~~~~~~~~~~~
load("results_WBScod_m1_A2_survey_No_One_noprofile_.RData")

sur <- as.list(env1$sdr,"Estimate"); sur  <- as.data.frame(sur$eta_density)
#sur <- as.list(env1$sdr,"Std. Error"); sur  <- as.data.frame(sur$eta_density)

rm(list=setdiff(ls(), ls(pattern=c("com|sur|both"))))



#~~~~~~~~~~~~~~~~~~
# combined data
#~~~~~~~~~~~~~~~~~~
load("results_WBScod_m1_A2_both_No_One_noprofile_.RData")

both <- as.list(env1$sdr,"Estimate"); both  <- as.data.frame(both$eta_density)
#both <- as.list(env1$sdr,"Std. Error"); both  <- as.data.frame(both$eta_density)

rm(list=setdiff(ls(), ls(pattern=c("com|sur|both"))))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract and transform the abundance estimated values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Transform abundance values to 0-1 interval for better visualization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

both_ct <- apply(both,2,concTransform)
sur_ct <- apply(sur,2,concTransform)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate raster files from the estimated abundances
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We need to do this step if we want to produce the maps on ggplot.
# To do so, we first neet do create a dataframe with the grids long/lat and the abundance associated to this locations.
# Afterwards we need to create an "empty raster" because we have an unregularly spaced grid

load("results_WBScod_m1_A2_both_No_One_noprofile_.RData")


# Create a dataframe with the grid lon/lat and the associated abundance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_bothCT <- cbind(gr,both_ct)

df_surCT <- cbind(gr,sur_ct)

# spdf <- SpatialPixelsDataFrame(points=df[c("lon", "lat")], data = df,tolerance = 0.959642)
# test_df <- as.data.frame(spdf)
# colnames(test_df) <- c("value", "x", "y")

# Create an empty rasterfile based on the dataframe from the previous step
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e <- extent(as.matrix(gr))
r <- raster(e,ncol=55,nrow=55) #WBS cod
#r <- raster(e,ncol=45,nrow=45) #KAT cod

abulist <- list(df_bothCT,df_surCT)
names(abulist) <- c("df_bothCT","df_surCT")


# Rasterize the previously created raster - FIXME: OTIMIZE THAT IN A LOOP!!!!!!!!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convert df from wide to long format
df_bothCT2 <- gather(df_bothCT,key=YearQuarter,value=Abundance,-c(lon,lat))

df_surCT2 <- gather(df_surCT,key=YearQuarter,value=Abundance,-c(lon,lat))


df_bothCT2$YearQuarter <- as.factor(df_bothCT2$YearQuarter)
df_surCT2$YearQuarter <- as.factor(df_surCT2$YearQuarter)


# # WBS cod
df_bothCT2$YearQuarter <- factor(df_bothCT2$YearQuarter, levels = c(paste("V",1:48,sep="")))
df_surCT2$YearQuarter <- factor(df_surCT2$YearQuarter, levels = c(paste("V",1:48,sep="")))




dflist_bothCT2 <- split(df_bothCT2,list(df_bothCT2$YearQuarter))
dflist_surCT2 <- split(df_surCT2,list(df_surCT2$YearQuarter))





raster_bothCT <- list() 
for(i in 1:length(dflist_bothCT2)){
  raster_bothCT [[i]] <- rasterize(dflist_bothCT2[[i]][,c("lon","lat")], r, dflist_bothCT2[[i]][,"Abundance"], fun=mean)
  raster_bothCT [[i]] <- disaggregate(raster_bothCT[[i]],2,method="bilinear")
}
names(raster_bothCT) <- names(dflist_bothCT2) #name new list accordingly



raster_surCT <- list() 
for(i in 1:length(dflist_surCT2)){
  raster_surCT [[i]] <- rasterize(dflist_surCT2[[i]][,c("lon","lat")], r, dflist_surCT2[[i]][,"Abundance"], fun=mean)
  raster_surCT [[i]] <- disaggregate(raster_surCT[[i]],2,method="bilinear")
}
names(raster_surCT) <- names(dflist_surCT2) #name new list accordingly





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate GIFs of the predicted abundance with ggplot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Plotting raster on ggplot isn't very straightforward.
#In order to plot them, one needs first to convert the raster file into a dataframe.
#But before doing that, one needs to convert the raster to a spatialpixeldataframe (SPDF) and afterwards,
#convert the SPDF to a conventional dataframe (DF).
#I also created manually a vector of titles for all the 48 maps, and which will be included within the ggplot script.



# Convert rasterlayer to a DF (raster -> SPDF -> DF)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#r.spdf2_bothCT <- as(raster_bothCT, "SpatialPixelsDataFrame")
#r.df2 <- as.data.frame(r.spdf2)

# Combined
r.spdf2_bothCT <- list()
r.df2_bothCT <- list()
for(i in 1:length(raster_bothCT)){
  r.spdf2_bothCT[[i]] <- as(raster_bothCT[[i]],"SpatialPixelsDataFrame")
  r.df2_bothCT[[i]] <- as.data.frame(r.spdf2_bothCT[[i]])
}


# Survey
r.spdf2_surCT <- list()
r.df2_surCT <- list()
for(i in 1:length(raster_surCT)){
  r.spdf2_surCT[[i]] <- as(raster_surCT[[i]],"SpatialPixelsDataFrame")
  r.df2_surCT[[i]] <- as.data.frame(r.spdf2_surCT[[i]])
}




# WBS cod
titles <- as.factor(c("2005:Q1","2005:Q2","2005:Q3","2005:Q4",
                      "2006:Q1","2006:Q2","2006:Q3","2006:Q4",
                      "2007:Q1","2007:Q2","2007:Q3","2007:Q4",
                      "2008:Q1","2008:Q2","2008:Q3","2008:Q4",
                      "2009:Q1","2009:Q2","2009:Q3","2009:Q4",
                      "2010:Q1","2010:Q2","2010:Q3","2010:Q4",
                      "2011:Q1","2011:Q2","2011:Q3","2011:Q4",
                      "2012:Q1","2012:Q2","2012:Q3","2012:Q4",
                      "2013:Q1","2013:Q2","2013:Q3","2013:Q4",
                      "2014:Q1","2014:Q2","2014:Q3","2014:Q4",
                      "2015:Q1","2015:Q2","2015:Q3","2015:Q4",
                      "2016:Q1","2016:Q2","2016:Q3","2016:Q4"))



# Create shapefile for study area
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data("worldHiresMapEnv")
DK_coast_poly <- map("worldHires",  fill=TRUE, col="transparent",
                     plot=FALSE, xlim=c(9,15.5), ylim=c(54.5,58))
DK_coast_poly$names
IDs <- sapply(strsplit(DK_coast_poly$names, ":"), function(x) x[1])
DK_poly <- map2SpatialPolygons(DK_coast_poly, IDs=IDs,
                               proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

#####################################################################################################################



# Bind all single data-frames
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# But first we need to include the time-specific info on the columns

# Include TimeYear column in each df
names(r.df2_bothCT) <- titles
names(r.df2_surCT) <- titles

for(i in 1:length(r.df2_bothCT)){
  r.df2_bothCT[[i]]$TimeYear <- as.factor(paste(titles[i]))
}

for(i in 1:length(r.df2_surCT)){
  r.df2_surCT[[i]]$TimeYear <- as.factor(paste(titles[i]))
}


# Bind all dataframes
dboth <- do.call("rbind",r.df2_bothCT)
dsur <- do.call("rbind",r.df2_surCT)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For all plots in one page
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#dboth2 <- subset(dboth, TimeYear%in%c("2010:Q1","2010:Q2","2010:Q3","2010:Q4"))
#dboth2$TimeYear <- factor(dboth2$TimeYear)

BOTH <- ggplot() +
              geom_tile(data=dboth, aes(x=x,y=y,fill = layer)) +
              ggtitle(paste("Integrated data -",titles[1])) +
              scale_fill_gradientn(colours = rev(parula(99)[99:1]),limits=c(0,1)) +
              labs(x = "Longitude (°)", y="Latitude (°)") +
              guides(fill = guide_colorbar(barwidth = 1, barheight = 16,ticks=F, title=NULL,frame.colour="black",frame.linewidth = 1.5)) +
              coord_map(xlim = c(9,15.5), ylim = c(53.3,56.5)) +
              geom_polygon(data=DK_poly, aes(x=long, y=lat, group=group),fill="grey70", colour="black") + 
              theme_bw() +
                  theme(#legend.text = element_blank(),
                  axis.text = element_text(size=13),
                  legend.text = element_text(size=13,face="bold"),
                  axis.text.x = element_text(size=14),
                  axis.text.y = element_text(size=14),
                  axis.title.y = element_text(margin=margin(t = 0, r = 18, b = 0, l = 0),size=16),
                  axis.title.x = element_text(margin=margin(t=18,r=0,b=0,l=0),size=16),
                  plot.title = element_text(hjust = 0.5,size=16,face="bold"),
                  plot.subtitle =  element_text(hjust = 0.5,size=14,face="bold"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
                  labs(title="Integrated data", subtitle = "{closest_state}") +
                  transition_states(TimeYear,transition_length=48,state_length=1) #tansition_lenght=nlevels(TimeYear)
                  #transition_time(Year)

gifBOTH <- animate(BOTH,fps=10,duration = 30)







SUR <- ggplot() +
  geom_tile(data=dsur, aes(x=x,y=y,fill = layer)) +
  ggtitle(paste("Survey data -",titles[1])) +
  scale_fill_gradientn(colours = rev(parula(99)[99:1]),limits=c(0,1)) +
  labs(x = "Longitude (°)", y="Latitude (°)") +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 16,ticks=F, title=NULL,frame.colour="black",frame.linewidth = 1.5)) +
  coord_map(xlim = c(9,15.5), ylim = c(53.3,56.5)) +
  geom_polygon(data=DK_poly, aes(x=long, y=lat, group=group),fill="grey70", colour="black") + 
  theme_bw() +
  theme(#legend.text = element_blank(),
    axis.text = element_text(size=13),
    legend.text = element_text(size=13,face="bold"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(margin=margin(t = 0, r = 18, b = 0, l = 0),size=16),
    axis.title.x = element_text(margin=margin(t=18,r=0,b=0,l=0),size=16),
    plot.title = element_text(hjust = 0.5,size=16,face="bold"),
    plot.subtitle =  element_text(hjust = 0.5,size=14,face="bold"),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  labs(title="Survey data", subtitle = "{closest_state}") +
  transition_states(TimeYear,transition_length=48,state_length=1) #tansition_lenght=nlevels(TimeYear)
#transition_time(Year)

gifSUR <- animate(SUR,fps=10,duration = 30)


s_gif <- image_read(gifSUR)
b_gif <- image_read(gifBOTH)


#Save gis sepparately
magick::image_write(s_gif, path="C:/Users/mruf/Documents/LGNB/Demo/SUR_animation.gif")
magick::image_write(b_gif, path="C:/Users/mruf/Documents/LGNB/Demo/BOTH_animation.gif")


new_gif <- image_append(c(s_gif[1], b_gif[1]))


for(i in 2:48){
  combined <- image_append(c(s_gif[i], b_gif[i]))
  new_gif <- c(new_gif, combined)
}

#
#new_gif #FIXME: not displaying all the time-steps.
#magick::image_write(new_gif2, path="C:/Users/mruf/Documents/LGNB/Demo/full_animation2.gif")


## Acess https://ezgif.com/combine to combine manually the indivudal-saved gifs.