Phylo_Summary<-function(tree){

  genera<-  unique(sub("_.*", "", tree$tip.label))

  spp_vec<-c()
  spp_vec_con<-c()
  for (i in 1:length(genera)){
    spp_vec<-tree$tip.label[sub("_.*", "", tree$tip.label) %in% genera[i]]
    spp_vec_con<-rbind(spp_vec_con,paste(spp_vec, collapse = ', '))
  }


  spp_vec_csin<-cSplit(spp_vec_con, 'V1', sep=", ", type.convert=FALSE)
  spp_vec_csin<-as.matrix(spp_vec_csin)


  ##Obtain node number

  noden<-c()
  node<-c()
  for (i in 1:length(genera)){
    node<-getMRCA(tree, na.omit(spp_vec_csin[i,]))
    node_1<-replace(node,which(is.null(node)),0)
    noden<-c(noden,node_1)
  }

  ##Organize genus and node number
  organiza<-cbind(genera, noden)

  ###extract and obtain Nnumber and Ntips
  tree$node.label<-(length(tree$tip.label)+1: c(length(tree$tip.label)+tree$Nnode))

  nnodes<-c()
  nntips<-c()
  formanodes<-c()
  formatips<-c()

  for (i in 1:dim(organiza)[1]){
    stree<-extract.clade(tree, replace(organiza[,2][i], organiza[,2][i]==0, tree$node.label[1]))
    nnodes<-stree$Nnode
    nntips<-length(stree$tip.label)
    formanodes<-c(formanodes,nnodes)
    formatips<-c(formatips,nntips)
  }

  newdata<-cbind(organiza, formanodes, formatips)
  newdata<-newdata[,-1]
  newdata[formanodes ==(tree$node.label[1]-2)]<-0
  newdata[formatips ==(tree$node.label[1]+1)] <-0
  rownames(newdata)<-genera



  ##Number of species in each genus

  sppnum<-c()
  for (i in 1:dim(spp_vec_csin)[1]){
    ri<-length(na.omit(spp_vec_csin[i,]))
    sppnum<-c(sppnum,ri)
  }

  newdata<-cbind(newdata, sppnum)
  newdata<-as.data.frame(newdata)
  nsppintree<- replace(as.vector(newdata$formatips), as.vector(newdata$formatips)==0, 1)
  newdata$formatips<-nsppintree


  ##IF ELSE
  mon<-c()
  monoph<-c()
  for (i in 1:dim(newdata)[1]){
    mon<-ifelse(newdata[i,3] == newdata[i,4], "Monophyletic", ifelse( as.numeric(newdata[i,3]) < c(2*as.numeric(as.character(newdata[i,4]))), "Paraphyletic", "Poly"))
    monoph<-rbind(monoph,mon)
  }

  newdata<-cbind(newdata, monoph)

  newdataf<-as.data.frame(newdata)
  newdata_2<-newdataf
  #newdata_2 <-newdataf[as.numeric(newdataf[,3]) < c(2*as.numeric(as.vector(newdataf[,4]))), ] ##Include paraphyletic
  #newdata_2 <- subset(newdata , monoph == "Mono") ##Without paraphyletic
  colnames(newdata_2)<-c("node number", "# nodes", "# tips (clade)","# tips (species)", "Monophyly")

  newdata_2[newdata_2==0] <- NA

  return(newdata_2)

}
