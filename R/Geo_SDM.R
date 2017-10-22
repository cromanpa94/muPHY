Geo_SDM<-function(data=data, res=2.5){

  rasters_min<-list()
  rasters_total<-list()
  for(i in 1:length(unique(ma$species))){

    ma<-as.data.frame(data[,c("species","lon", "lat")])
    sma<-ma[ma$species==unique(ma$species)[i],]

    sma$species<-rep(1, dim(sma)[1])
    names(sma)[1]<-gsub(" ", ".", unique(ma$species)[i], fixed=TRUE)

    print("Looking for Bioclim raster files...")

    preds <- raster::getData("worldclim",var="bio",res=res)

    print("Cropping rasters and removing collinearity...")

    occurrence<-sma
    coord_cols <- match(c('lon', 'lat'), colnames(occurrence))
    occurrence <- SpatialPointsDataFrame(occurrence[, coord_cols],
                                         occurrence,
                                         coords.nrs  = coord_cols)

    ##Crop rasters
    rc <- crop(preds, extent(gBuffer(occurrence, width = 2)))

    ##
    print("Generating SDM")
    list_unco<-removeCollinearity(rc, select.variables = TRUE)

    kuu <- sdmData(as.formula(paste(names(sma)[1],"~", paste(list_unco, collapse="+")  )),
                   train=occurrence,predictors = rc, bg=list(n=1000,method='gRandom',remove=TRUE))
    m1 <- sdm(as.formula(paste(names(sma)[1],"~", paste(list_unco, collapse="+")  )),data=kuu,methods=c('glm','brt',"gam"))
    e1 <- ensemble(m1,newdata=rc,setting=list(method='weighted',stat='AUC'))
    plot(e1, main=paste(names(sma)[1], "All"))
    rasters_total[[i]]<-e1

    e1[e1 < cellStats(e1, max)*0.1] <- NA
    plot(e1, main=paste(names(sma)[1], ">30%") )

    rasters_min[[i]]<-e1

    cat("Done: Species",i,"of", length(unique(ma$species)), "\n")

  }

  return(rasters_min)


}
