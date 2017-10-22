Phylo_Get_Div_Single<-function(tree=tree,data=newdata_2, e.frac=0.5){

  newtree<-tree

  newdata_3<-data
  newdata_3[] <- lapply(newdata_3, as.character)
  newdata_3[] <- lapply(newdata_3, as.numeric)

  diversification<-c()
  div<-c()

  for (i in 1:dim(newdata_3)[1]){

    div<-if(is.na(newdata_3[i,2]) == T) bd.ms(phy=NULL, time= as.numeric(inf[grep(grep(row.names(newdata_3)[i], tree$tip.label,value=TRUE), inf[,1]),][2]), e=e.frac
                                              , n= newdata_3[i,6], missing=(newdata_3[i,6]-1), crown = F) else bd.ms( e=e.frac,phy=extract.clade(newtree, newdata_3[i,1], root.edge=4), missing=(newdata_3[i,6]-newdata_3[i,4]), crown=F)
    diversification<-rbind(diversification, div)
    cat("\r", i, "of", dim(newdata_3)[1])
    flush.console()
  }

  diversification_2<-diversification
  diversification_2<-replace(diversification_2,which(is.nan(diversification_2)),0)
  newdata_4<-cbind(newdata_3, diversification_2)
  colnames(newdata_4)[7]<-"diversification"

  print(newdata_4)

}
