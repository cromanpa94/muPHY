Phylo_Get_Div<-function(tree2=NULL,richness_tur=NULL){

  ifelse(is.null(tree2$node.label)==T,
  tree2$node.label<-seq(length(tree2$tip.label)+1, (length(tree2$tip.label)+tree2$Nnode)),
  tree2$node.label)

  library(phytools)
  library(geiger)
  #####Div_2
  div_crown<-list()
  div_stem<-list()
  for( i in 1:length(richness_tur[,1]) ){
    tryCatch({

      if(length(tree2$tip.label[grep(as.character(richness_tur[i,1]),tree2$tip.label)])==1){
        miss<-richness_tur[i,2]-1
        n<-richness_tur[i,2]
        nu<-length(tree2$tip.label)
        ee<-setNames(tree2$edge.length[sapply(1:nu,function(x,y)   which(y==x),y=tree2$edge[,2])],tree2$tip.label)
        age_mono<-as.numeric(ee[grep(as.character(richness_tur[i,1]),names(ee))])
        div_stem[[i]]<-data.frame(e0=bd.ms(phy=NULL,time=age_mono ,n, missing = miss, crown=FALSE, epsilon = 0),
                                  e0.5=bd.ms(phy=NULL,time=age_mono ,n, missing = miss, crown=F, epsilon = 0.5),
                                  e0.9=bd.ms(phy=NULL,time=age_mono ,n, missing = miss, crown=F, epsilon = 0.9),genus=richness_tur[i,1])

      } else{


        n=length(extract.clade(tree2, findMRCA(tree2,  tree2$tip.label[grep(as.character(richness_tur[i,1]),tree2$tip.label)]))$tip.label)
        miss<-richness_tur[i,2]-n
        sb_tree<-extract.clade(tree2, findMRCA(tree2,  tree2$tip.label[grep(as.character(richness_tur[i,1]),tree2$tip.label)]))
        node<-findMRCA(tree2,  tree2$tip.label[grep(as.character(richness_tur[i,1]),tree2$tip.label)])
        age_crown<-as.numeric(branching.times(tree2)[which(names(branching.times(tree2))==node)])
        H<-nodeHeights(tree2)
        age_stem<-H[tree2$edge==phytools:::getAncestors(tree2,node,"parent")][1]

        div_crown[[i]]<-data.frame(e0=bd.ms(phy=NULL,time=age_crown, n, missing = miss, crown=TRUE, epsilon = 0),
                                   e0.5=bd.ms(phy=NULL,time=age_crown, n, missing = miss, crown=TRUE, epsilon = 0.5),
                                   e0.9=bd.ms(phy=NULL,time=age_crown,n,  missing = miss, crown=TRUE, epsilon = 0.9),genus=richness_tur[i,1])

        div_stem[[i]]<-data.frame(e0=bd.ms(phy=NULL,time=age_stem ,n, missing = miss, crown=FALSE, epsilon = 0),
                                  e0.5=bd.ms(phy=NULL,time=age_stem ,n, missing = miss, crown=F, epsilon = 0.5),
                                  e0.9=bd.ms(phy=NULL,time=age_stem ,n, missing = miss, crown=F, epsilon = 0.9),genus=richness_tur[i,1])

      }

    }, error=function(e){})

    print(i)

  }

  null_rows_S<-div_stem[(sapply(div_stem, is.null)==F)]
  null_rows_C<-div_crown[(sapply(div_crown, is.null)==F)]
  div.rates <- cbind.data.frame(do.call("rbind", null_rows_S),do.call("rbind", null_rows_C))
  return(div.rates)

}
