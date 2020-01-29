#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Covariates for prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# Extract covariates for the spatial grid 
# based on grid defined in the LGNB_Rmodel script


# 1) Define spatial extent and coords
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ext <- extent(8,16,53.5,58.5) #  ~extent of the study area
gr_coords <- data.frame(lon=gr$lon, lat=gr$lat) # Take the grid object from the "model_complete" script


# 2) Depth
#~~~~~~~~~~~~
bathy <- raster("C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/PhD/Environmental predictors/bathymetry.tif")
depth <- crop(bathy, ext) ; res(depth)
depth <- depth*-1

depth_gr <- extract(depth,gr_coords) # extract depht for the grid

summary(depth_gr) # some NA's - replaced them by the median value
depth_gr[is.na(depth_gr)] <- median(depth_gr,na.rm=T)



# 3) Sediment type
#~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/PhD/Environmental predictors/EUNIS")
eunis <- readOGR(".", "Export_Output")
eunis2 <- eunis[,"EUNIS_name"]; rm(eunis)
map <- spTransform(eunis2,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")); rm(eunis2)


coordinates(gr_coords) <- ~ lon + lat
proj4string(gr_coords) <- proj4string(map)
sediments_gr <- over(gr_coords, map[,"EUNIS_name"]) #This tackes a long time (~15-20min)

sediments_gr$EUNIS_name <- factor(sediments_gr$EUNIS_name)
colnames(sediments_gr) <- "sediment"


levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Atlantic and Mediterranean high energy circalittoral rock"] <- "rock"
#levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Atlantic and Mediterranean high energy infralittoral rock"] <- "rock"
#levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Atlantic and Mediterranean low energy infralittoral rock"] <- "rock"
#levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Atlantic and Mediterranean moderate energy circalittoral rock"] <- "rock"
#levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Atlantic and Mediterranean moderate energy infralittoral rock"] <- "rock"
levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Baltic moderately exposed circalittoral rock"] <- "rock"
levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Baltic moderately exposed infralittoral rock"] <- "rock"
levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Data not available"] <- NA
#levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Deep-sea mud"] <- "mud"
levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Sublittoral coarse sediment"] <- "mixed_sediment"
levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Sublittoral mixed sediments"] <- "mixed_sediment"
levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Sublittoral mud"] <- "mud"
levels(sediments_gr$sediment)[levels(sediments_gr$sediment)=="Sublittoral sand"] <- "sand"

sedna <- which(is.na(sediments_gr$sediment))

sediments_gr$sediment <- as.character(sediments_gr$sediment)
sediments_gr$sediment[is.na(sediments_gr$sediment)] <- "rock"
sediments_gr$sediment <- as.factor(sediments_gr$sediment)



## locate grid ID
gr_gf <- gridLocate(gr, gr) ## CHECK!!!!!
covariates <- data.frame(depth = depth_gr, sediment = sediments_gr$sediment, gf = gr_gf)


saveRDS(covariates,"C:/Users/mruf/Documents/LGNB/Data/pred_covariates.rds")

