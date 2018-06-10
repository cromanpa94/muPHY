get_ortho_clusters<-function(ingroup, outgroup, MSA){

  ##Delete files
  if(file.exists("sequences.fasta")){unlink("sequences.fasta", recursive = T)}else{}
  if (length(grep("cluster_*", list.files("."), value = T)) == 0) {
  } else {
    unlink(grep("cluster_*", list.files("."), value = T))
  }


  taxa<-c(ingroup, outgroup)


  lapply(1:length(taxa), function(x){
    tr<-entrez_search(db="nuccore", term= paste0(taxa[x], "[ORGN]"), use_history=TRUE)
    for( seq_start in seq(1,tr$count,tr$count/10)){
      recs <- entrez_fetch(db="nuccore", web_history=tr$web_history,
                           rettype="fasta", retstart=seq_start)
      cat(recs, file="sequences.fasta", append=TRUE)
      cat(round(seq_start+49), "sequences downloaded\r")
    }

  })

  seqs<- read.dna("sequences.fasta", "fasta")

  seqs<- seqs[! duplicated(names(seqs)) ]
  seqs<- seqs[- grep("sp.|cf.|aff.", names(seqs) ) ]

  write.dna(seqs, "sequences.fasta", format = "fasta")

  args <-
    paste0("-in ", getwd(), "/sequences.fasta", " -input_type fasta -dbtype nucl")
  system2(
    command = 'makeblastdb',
    args = args,
    wait = FALSE,
    stdout = TRUE
  )


  que <-
    paste0("-outfmt 6 -query sequences.fasta -out test.txt -db ",
           getwd(),
           "/sequences.fasta", " -num_threads 4")
  system2(
    command = 'blastn',
    args = que,
    wait = FALSE,
    stdout = TRUE
  )


  BLAST1 <- read.blast(file = "test.txt")
  sim_mat<-simMatrix(BLAST1)

  res.d <- sim2dist(sim_mat)
  hc<-hclust(res.d)
  hc2 <- cutree(hc, h=0.999999)

  ##Write clusters

  seqs<- read.dna("sequences.fasta", "fasta")
  lapply(1:length(unique(hc2)), function(x){
    names<-sapply(strsplit(names(seqs), " "), head, 1)
    seqs2<-seqs[ which(names %in% names(which(hc2==x)))]

    cat(file=paste0("cluster_",x,".fasta"), paste(paste0(">",names(seqs2)),
                                                  sapply(as.character(seqs2), paste, collapse=""), sep="\n"), sep="\n")

  })

  ##Do MSA--------
  if(MSA==TRUE){

    cat("\nalignment in process. Please be patient\n")

    temp = list.files(pattern="*.fasta")
    temp<-temp[grep("cluster", temp)]


    ##Chose method: "Muscle", "ClustalOmega", "ClustalW"

    align<-function(file.names){
      mySequences <- readDNAStringSet(file.names)
      if(length(mySequences)<2){}else{
        myFirstAlignment <- msa(mySequences, "ClustalOmega")
        aln <- as.DNAbin(msaConvert(myFirstAlignment, type="phangorn::phyDat"))}
    }


    alignments<-pblapply(temp,align)

    alignments<-Filter(Negate(is.null), alignments)

    is.null(alignments)
    options(warn=-1)
    dir.create( "Aligned")
    options(warn=-0)


    for (i in 1: length(alignments)){
      write.dna(alignments[[i]], file=paste0("Aligned/Alignment_","Cluster" ,i, ".fasta"), format="fasta")
    }

    temp = list.files(path = "Aligned/",pattern="*.fasta", full.names = T)

    fst<- lapply(temp, read.dna, format = "f")

    dfs<-  do.call("rbind",  lapply(1:length(fst), function(x){
      d<-as.data.frame(cbind(x, t(sapply(strsplit(labels(fst[[x]]), " "), head, 3))))
      colnames(d)<-c("Cluster", "AN", "Genus", "Species")
      d
    }))

    write.csv(dfs, "Aligned/Sampling.csv")



  } else{}
  ##Sampling df
  temp = list.files(pattern="*.fasta")
  temp<-temp[grep("cluster", temp)]

  fst<- lapply(temp, read.dna, format = "f")

  dfs<-  do.call("rbind",  lapply(1:length(fst), function(x){
    d<-as.data.frame(cbind(x, t(sapply(strsplit(labels(fst[[x]]), " "), head, 3))))
    colnames(d)<-c("Cluster", "AN", "Genus", "Species")
    d
  }))

  write.csv(dfs, "Sampling.csv")
}
