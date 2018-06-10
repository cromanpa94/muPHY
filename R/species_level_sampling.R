species_level_sampling<-function(pattern="Cluster"){

  temp = list.files(pattern="*.fasta")
  temp<-temp[grep(pattern, temp)]
  temp<-temp[order(nchar(temp), temp)]

  fst<- lapply(temp, read.dna, format = "f")
  dfs<-  do.call("rbind",  lapply(1:length(fst), function(x){
    d<-as.data.frame(cbind(x, t(sapply(strsplit(labels(fst[[x]]), " "), head, 3))))
    colnames(d)<-c("Cluster", "AN", "Genus", "Species")
    d
  }))
  ntotal<-unique(paste(dfs$Genus, dfs$Species))
  ts<-unlist(lapply(1:max(as.numeric(dfs$Cluster)), function(x){
    lspp<-unique(paste(dfs[dfs$Cluster == unique(dfs$Cluster)[x],c(3)],
                       dfs[dfs$Cluster == unique(dfs$Cluster)[x],c(4)]))

    c(length(lspp)/length(ntotal))

  }))

  dir.create("Species_level")

  spp_lev<-do.call("rbind",lapply(1:max(as.numeric(dfs$Cluster)), function(x){
    lspp<-duplicated(paste(dfs[dfs$Cluster == unique(dfs$Cluster)[x],c(3)],
                           dfs[dfs$Cluster == unique(dfs$Cluster)[x],c(4)]))

    tan<-dfs[dfs$Cluster == unique(dfs$Cluster)[x],]

    jhn<-fst[[x]][lspp==F,]

    dimnames(jhn)[[1]]<-paste0(tan[lspp==F,3], "_",tan[lspp==F,4])

    write.dna(jhn, paste0("Species_level/Cluster",x,"_sf=",round(ts[x],4),".fasta"), format = "fasta")
    return(tan[lspp==F,])

  }))
  write.csv(spp_lev,"Species_level/Sampling.csv")

}
species_level_sampling()
