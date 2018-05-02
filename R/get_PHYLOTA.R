get_PHYLOTA<-function(clade, MSA = FALSE, ALI =FALSE, c_directory=FALSE){

  cat("Starting...")

  ##GetWD
  mainDir <- ifelse(c_directory == FALSE, tempdir() , getwd())
  setwd(mainDir)
  ti<-get_uid(sciname = clade)[1]

  ##Get Main website
  url <- paste0("http://phylota.net/cgi-bin/sql_getclusterset.cgi?ti=",ti, "&ntype=1&piflag=1&dflag=0&db=194")
  doc <- htmlParse(url)
  links <- xpathSApply(doc, "//a/@href")
  clusters<-grep("getcluster",links, value=T )

  clusters<-clusters[grep(ti,clusters)]

  ##Sub-websites containing the different sequences

  getwebsites<- function( websites ){
    url_2 <- websites
    doc_2 <- htmlParse(url_2)
    links_2 <- xpathSApply(doc_2, "//a/@href")
    todo<-grep("getalign",links_2, value=T )
    matches <- regmatches(todo, gregexpr("[[:digit:]]+", todo))
    num<-as.numeric(unlist(matches))
    newweb<-paste0("http://phylota.net/cgi-bin/sql_getcluster_fasta.cgi?format=gi&db=194&ti=",num[1], "&cl=",num[2],"&ntype=1" )
  }

  cat("Retrieving clusters")
  seqs<-pblapply(clusters,getwebsites)

  ##Get sequences per cluster

  ReadFasta<-function(file) {
    options(warn=-1)
    fasta<-readLines(file)
    options(warn=0)
    fasta[1]<-gsub("<html><pre>", "", fasta[1])
    fasan<-fasta
    fasan_2<-fasan[-length(fasta)[1]]
    # Identify header lines
    ind<-grep(">", fasan_2)
    # Identify the sequence lines
    s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasan_2)))
    # Process sequence lines
    seqs<-rep(NA, length(ind))
    for(i in 1:length(ind)) {
      seqs[i]<-paste(fasan_2[s$from[i]:s$to[i]], collapse="")
    }
    # Create a data frame
    DF<-data.frame(name=gsub(">", "", fasan_2[ind]), sequence=seqs)
    # Return the data frame as a result object from the function
    return(DF)
  }
  options(warn=-1)
  dir.create(file.path(mainDir, "Unaligned"))
  options(warn=-0)

  setwd(file.path(mainDir, "Unaligned"))

  plot.progress <- function(...)	{
    vectOfBar <- c(...)*100
    numOfBar <- length(vectOfBar)
    plot(c(0,100), c(0,numOfBar), type='n', xlab='', ylab='', yaxt='n', mar=c(3,3,3,3))
    for(i in 1:numOfBar) {
      rect(0, 0.1+i-1, vectOfBar[i], 0.9+i-1, col=rainbow(numOfBar)[i])
      text(0.5, 0.5+i-1, paste('Running!!) ', i, ': ', round(vectOfBar[i],2), '%', sep=''), adj=0)
    }
    title('Progress...')
  }

  cat("Downloading clusters. Look at your working directory")


  for (i in 1:length(seqs)){  #
    sequences<-ReadFasta(gsub("format=gi","format=all",seqs[i]$href))
    sequences$gi<- word(sequences$name, 1, 1, fixed("_"))
    sequences$name<- gsub(" ", "_" ,word(sequences$name, 3, 3, fixed("_")))

    summary_unique_sequences<-sequences[!duplicated(sequences$name), ]
    summary_non_sp<-summary_unique_sequences[- grep("_sp.", summary_unique_sequences$name) , ]

    if(length(grep("_sp.", summary_unique_sequences$name))== 0) {
      seqRFLP::write.fasta(dataframe2fas(summary_unique_sequences[c("name","sequence")]), file=paste0(clade, "_","Cluster" ,i, ".fasta"))
      write.csv(summary_unique_sequences[c("name","gi")],
                file=paste0(clade, "_","Cluster" ,i, ".csv") )

    } else{
      seqRFLP::write.fasta(dataframe2fas(summary_non_sp[c("name","sequence")]), file=paste0(clade, "_","Cluster" ,i, ".fasta"))
      write.csv(summary_non_sp[c("name","gi")],
                file=paste0(clade, "_","Cluster" ,i, ".csv"))
    }
    summary_non_sp<-NULL
    summary_unique_sequences<-NULL
    plot.progress(i/length(seqs))

  }

  cat("\nDone!!")
  setwd(mainDir)



  ##Make sampling matrix------
  #library(plyr)
  #sampling_matrix<- join_all(for_Sup[c(23:24)], "species", match="first")

  ##Do MSA--------
  ##Read fasta files
  if(MSA==TRUE){

    cat("\nalignment in process. Please be patient\n")

    setwd(file.path(mainDir, "unaligned"))

    temp = list.files(pattern="*.fasta")

    ##Chose method: "Muscle", "ClustalOmega", "ClustalW"

    align_PHYLOTA<-function(file.names){
      mySequences <- readDNAStringSet(file.names)
      myFirstAlignment <- msa(mySequences, "ClustalOmega")
      aln <- as.DNAbin(msaConvert(myFirstAlignment, type="phangorn::phyDat"))
    }

    alignments<-pblapply(temp,align_PHYLOTA)

    options(warn=-1)
    dir.create(file.path(mainDir, "Aligned"))
    options(warn=-0)

    setwd(file.path(mainDir, "Aligned"))

    for (i in 1: length(alignments)){
      write.dna(alignments[[i]], paste0(clade, "_", "Alignment","Cluster" ,i, ".fasta"), format="fasta")
    }
    cat("\nAligned and unaligned sequences are in separated folders within your working directory")

    if (ALI==TRUE){
      ####For Aliscore------

      setwd(mainDir)

      if( "Aliscore_v.2.0" %in% list.files() == FALSE){
      download.file("http://software.zfmk.de/ALISCORE_v2.0.zip",'Aliscore_v.2.0.zip')
        unzip("Aliscore_v.2.0.zip")
        unzip("ALISCORE_v2.0/Aliscore_v.2.0.zip")
        } else
      options(warn=-1)
      dir.create(file.path(mainDir, "Aliscore"))
      options(warn=0)

      setwd(file.path(mainDir, "Aliscore"))

      aln_aliscored<-list()
      for (i in 1: length(alignments)){
        aln_aliscored[[i]]<-aliscore(alignments[[i]], gaps = "ambiguous", w = 3, path = paste0(mainDir, "/Aliscore_v.2.0"))
        write.dna(aln_aliscored[[i]], paste0(clade, "_", "AliScore", "_Alignment","_Cluster_" ,i, ".fasta"), format="fasta")
        plot.progress(i/length(alignments))
        cat("\nCurated alignments are under ALISCORE folder")

      }

    } else {}

  } else {}

  setwd(mainDir)
  cat("\nDone!!")


}
