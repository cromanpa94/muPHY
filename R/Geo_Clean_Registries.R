Geo_Clean_Registries<-function(name= file,ncores=6){


  matrixpr <- fread(file, header=T)
  names(matrixpr)[18]<-"lon"
  names(matrixpr)[17]<-"lat"

  matrixpr<-matrixpr[!institutioncode=="iNaturalist",]
  matrixpr<-matrixpr[!basisofrecord=="UNKNOWN",]


  matrixpr_2<-matrixpr[!with(matrixpr,is.na(lon)& is.na(lat)),]
  matrixpr_2 <- matrixpr_2[(matrixpr_2$genus %in% row.names(newdata_4)),]
  matrixpr_2<-as.data.frame(matrixpr_2)

  # install.packages("rangeBuilder");


  print("Filtering registries on land")

  ncores <- ncores

  onLand <- filterByLand(matrixpr_2[,c('lon','lat')], returnGood = TRUE)
  crotalus <- matrixpr_2[onLand,]
  crotalus <- cbind(crotalus, match = crotalus$species)
  crotalus <- crotalus[which(!is.na(crotalus$match)),]
  crotalus$match <- as.character(crotalus$match)
  crotalus <- cbind(crotalus,
                    matchedCountries=standardizeCountry(crotalus$country, nthreads = ncores))
  crotalus$matchedCountries <- as.character(crotalus$matchedCountries)

  #for each record, do the coordinates actually fall in the listed country?
  #for records with no country, fill that country in
  #if returned country does not match country, check for sign errors with flipSign
  pb <- txtProgressBar(min = 0, max = nrow(crotalus), style = 3)
  for (i in 1:nrow(crotalus)) {
    setTxtProgressBar(pb, i)
    x <- closestCountry(crotalus[i, c('lon','lat')])
    if (!crotalus$matchedCountries[i] %in% x & crotalus$matchedCountries[i] != '') {
      sign <- flipSign(crotalus[i, c('lon','lat')],
                       country = crotalus$matchedCountries[i])
      if (sign$matched == TRUE) {
        crotalus[i, c('lon', 'lat')] 	<- as.numeric(sign$newcoords)
      }
    } else if (crotalus$matchedCountries[i] == '') {
      crotalus$matchedCountries[i] <- x[1]
    }
  }

  ##remove registries from where the species is invasive
  #install.packages("originr")
  library(originr)
  print("\n Checking for invasive species")

  corrected<-NULL
  dataa<-NULL
  countries<-NULL
  dars<-NULL
  newspw<-list()

  pb <- txtProgressBar(min = 0, max = length(unique(crotalus$species)), style = 3)
  for (i in 1:length(unique(crotalus$species))){
    species<-if(unique(crotalus$species)[i] =="") "Silverstoneia nubicola" else unique(crotalus$species)[i]
    dataa<-crotalus[(crotalus$species == species),]
    countries<-unique(dataa$matchedCountries)
    dars<-griis(name = species)
    newspw[[i]] <- if(unique(crotalus$species)[i] =="") NA else {dataa[!(dataa$matchedCountries %in% dars$Country),]}
    setTxtProgressBar(pb, i)
  }

  na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

  newspw<-na.omit.list(newspw)


  corrected <- ldply(newspw, data.frame)


  crotalus<-corrected


  crotalusn<-crotalus %>% group_by(species) %>% filter (! duplicated(lon))
  crotalusn<-crotalus %>% group_by(species) %>% filter (! duplicated(lat))
  tt <- table(crotalusn$species)

  print("\n Filtering species with at least 10 occurrences")

  head(crotalusn)
  write.csv(crotalusn,"Curated.registries.csv")
  Curated.DB <- subset(crotalusn, species %in% names(tt[tt > 5]))

  return(Curated.DB)
}
