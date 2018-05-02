update_phylota <-
  function(lineage,
           nsamples = 5,
           genes = NULL,
           MSA = TRUE,
           ALI = TRUE,
           outgroup = NULL,
           correct_db = TRUE,
           delete_all = TRUE,
           c_directory=FALSE) {
    ##New version
    options(warn = -1)

    Clade <- lineage
    clade <- lineage
    fn <- "Unaligned"
    fn2 <- "Aligned"
    fn2.1 <- "Aliscore"
    fn3 <- "phylota.sampling.csv"
    fn4 <- "included.species.csv"
    if (file.exists(fn))
      unlink(fn, recursive = T)
    if (file.exists(fn2))
      unlink(fn2, recursive = T)
    if (file.exists(fn2.1))
      unlink(fn2.1, recursive = T)
    if (file.exists(fn3))
      unlink(fn3, recursive = T)
    if (file.exists(fn4))
      unlink(fn4, recursive = T)
    if (length(grep("Molecular_*", list.files("."), value = T)) == 0) {
    } else {
      unlink(grep("Molecular_*", list.files("."), value = T))
    }



    cat("\n Get PhyLota clusters \n")

    rPhyloTa::get_PHYLOTA(lineage, MSA = F, ALI = F, c_directory= c_directory)
    mainDirect <- ifelse(c_directory == FALSE, tempdir() , getwd())
    mainDir<- ifelse(c_directory == FALSE, tempdir() , getwd())

    subwd <- paste0(mainDirect, "/Unaligned/")

    ##Check for new species in genbank

    cat("Checking for novel species in genbank \n")

    spp_genbank <-
      taxize::downstream(lineage, db = "ncbi", downto = 'species')[[1]]

    spp_sampled <-
      do.call(rbind, lapply(list.files(
        path = subwd,
        pattern = c(".csv"),
        full.names = T
      ), read.csv))
    spp_sampled <- spp_sampled[!duplicated(spp_sampled$name), ]
    new_spp <-
      setdiff(gsub(" ", "_", spp_genbank$childtaxa_name), spp_sampled$name)

    if (length(new_spp) == 0) {
      if (file.exists(fn))
        unlink(fn, recursive = T)

      stop("Nothing to be added")
    } else{
      cat("Tranforming each cluster into a Blast DB \n")

      ##First make DB for each cluster
      files <- list.files(path = subwd, pattern = c("Cluster"))
      files <- Filter(function(x)
        grepl("\\.fasta$", x), files)

      subwd<- if(Sys.info()[['sysname']]=="Windows") {
        gsub("\\", "/", subwd, fixed=T) } else {subwd}

      for (i in 1:length(files)) {
        args <-
          paste0("-in ", subwd, files[i], " -input_type fasta -dbtype nucl")
        system2(
          command = 'makeblastdb',
          args = args,
          wait = FALSE,
          stdout = TRUE
        )
        Sys.sleep(1)
        print(i)
      }



      cat(
        length(new_spp),
        " species can be added to the ",
        dim(spp_sampled)[1],
        "that are already sampled in PhyloTa \n"
      )

      ##Check which genes do I have in PhyloTa clusters. This is based on 5 species per each cluster.

      if (is.null(genes) == T) {
        cat("Looking for genes sampled in PhyLota \n")

        sample_Gi <- do.call(rbind, lapply(list.files(path = subwd,
                                                      pattern = c("csv"), full.names = T), read.csv, nrow = nsamples))

        sampled_genes <- ncbi_byid(gsub("gi", "", as.vector(sample_Gi$gi)))

        choosebank("genbank")
        gene_fin <- list()
        for (i in 1:length(sampled_genes$acc_no)) {
          Dengue1 <-
            query("Dengue1", paste0("AC=", gsub(
              "\\..*", "", sampled_genes$acc_no[i]
            )))
          annots <- getAnnot(Dengue1$req[[1]])
          a <- gsub("/gene=\"", "", grep("/gene=", annots, value = T))
          b <- gsub("\"", "", a)
          gene_fin[[i]] <- gsub("[[:space:]]", "", b)[1]
          cat("New gene symbol found!",
              gsub("[[:space:]]", "", b)[1],
              "\n")

        }

        closebank()
        sampled_gene_names <-levels(factor(unlist(gene_fin)))
        cat("\n Your sampling includes", levels(factor(unlist(gene_fin))), "loci \n")
      } else{
        sampled_gene_names <- genes
      }


      ##Then run each sequence against each cluster
      #some functions
      grab.results <- function (term) {
        # Search for the given term on nuccore. This gives us a list of
        # record IDs.
        ids <- esearch(term, db="nuccore")

        # Grab summaries for the given record IDs, as a sort-of data frame.
        sum <- esummary(ids, db="nuccore")
        data <- content(sum, as="parsed")

        # For some reason, this parser gives us lists of lists instead of a
        # proper data frame (which should be lists of vectors). Return a
        # fixed-up version.
        data.frame(lapply(data, as.character), stringsAsFactors=FALSE)
      }

      #######################################################################

      # Takes a dataframe of Genbank results and a keyword string and returns a
      # string containing an accession number. This is the accession number of
      # the first result in the table that had the keyword in its Title field.

      first.of.type <- function (results, type) {
        # Filter rows whose title contains the given text
        filtered <- results[grepl(type, results$Title, fixed=TRUE),]

        # Return the first accession number
        filtered$OSLT[1]
      }

      #######################################################################

      # Given a list of species and genes, return a dataframe of accession numbers.
      #
      # Takes string vectors as arguments and returns a dataframe. The rows of the
      # dataframe are species. There is a column for each gene. These columns are
      # filled with accession numbers. The species names are put in a column called
      # "Species".

      scrape.genbank <- function (species, genes) {
        # Create dataset skeleton
        data <- data.frame(Species=species)
        for (i in genes) {
          data[,i] <- rep(NA, length(species))
        }

        # Look up data for each species
        for (i in 1:length(species)) {
          n <- species[i]
          print(sprintf("Looking up for %s (%d/%d)...", n, i, length(species)))
          tryCatch({
            r <- grab.results(paste(n, "[Primary Organism]"))
            for (g in genes) {
              data[i,g] <- first.of.type(r, g)
            }
          }, error = function(e) { })
        }

        data
      }

      #######################################################################


      # Genbank scraping functions
      # Copyright (C) 2015 Anusha Beer <anbeer29@gmail.com>
      #
      # Permission to use, copy, modify, and/or distribute this software for
      # any purpose with or without fee is hereby granted, provided that the
      # above copyright notice and this permission notice appear in all copies.
      #
      # THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
      # WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
      # MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
      # ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
      # WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
      # ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
      # OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

      # Search Genbank for search term e.g. species name. Returns a dataframe
      # of all available (nucleotide) sequences on Genbank.
      #
      # Give it a string e.g. species name. Returns a dataframe. Each row is a
      # result and each column is an attribute of the results. Some of the imprtant
      # attributes are "Title", "OSLT". Title contains a description of the
      # sequence. OSLT contains acession number.

      grab.results <- function (term) {
        # Search for the given term on nuccore. This gives us a list of
        # record IDs.
        ids <- esearch(term, db="nuccore")

        # Grab summaries for the given record IDs, as a sort-of data frame.
        sum <- esummary(ids, db="nuccore")
        data <- content(sum, as="parsed")

        # For some reason, this parser gives us lists of lists instead of a
        # proper data frame (which should be lists of vectors). Return a
        # fixed-up version.
        data.frame(lapply(data, as.character), stringsAsFactors=FALSE)
      }

      #######################################################################

      # Takes a dataframe of Genbank results and a keyword string and returns a
      # string containing an accession number. This is the accession number of
      # the first result in the table that had the keyword in its Title field.

      first.of.type <- function (results, type) {
        # Filter rows whose title contains the given text
        filtered <- results[grepl(type, results$Title, fixed=TRUE),]

        # Return the first accession number
        filtered$OSLT[1]
      }

      #######################################################################

      # Given a list of species and genes, return a dataframe of accession numbers.
      #
      # Takes string vectors as arguments and returns a dataframe. The rows of the
      # dataframe are species. There is a column for each gene. These columns are
      # filled with accession numbers. The species names are put in a column called
      # "Species".

      scrape.genbank <- function (species, genes) {
        # Create dataset skeleton
        data <- data.frame(Species=species)
        for (i in genes) {
          data[,i] <- rep(NA, length(species))
        }

        # Look up data for each species
        for (i in 1:length(species)) {
          n <- species[i]
          print(sprintf("Looking up for %s (%d/%d)...", n, i, length(species)))
          tryCatch({
            r <- grab.results(paste(n, "[Primary Organism]"))
            for (g in genes) {
              data[i,g] <- first.of.type(r, g)
            }
          }, error = function(e) { })
        }

        data
      }
      ##Start!
      ##
      new_spp2 <-   if (is.null(outgroup) == T) {
        new_spp
      } else{
        new_spp <- c(new_spp, outgroup)
      }

      data <-
        scrape.genbank(gsub("_", " ", new_spp), sampled_gene_names)



      ###Now, I download each sequence and test where it fits better

      cat("\nFitting the new sequences into the existing clusters \n")

      acce_new <- as.character(na.omit(unlist(data[, -1])))

      if (length(acce_new) == 0) {
        stop("No new sequences found")
      }

      n_clust <- length(list.files(path = subwd, pattern = c("csv")))
      species_included <- list()
      for (i in 1:length(acce_new)) {
        sequence <- read.GenBank(acce_new[i])
        Acc_spp <- names(sequence)
        names(sequence) <- attr(sequence, "species")

        for (j in 1:n_clust) {
          write.dna(sequence, "sequence.fasta", format = "fasta")
          que <-
            paste0("-outfmt 6 -query sequence.fasta -out test.txt -db ",
                   subwd,
                   files[j])
          system2(
            command = 'blastn',
            args = que,
            wait = FALSE,
            stdout = TRUE
          )

          if (file.info("test.txt")$size > 0) {
            cluster <- read.dna(paste0(subwd, files[j]), format = "fasta")
            cluster[length(cluster) + 1] <- sequence
            names(cluster)[length(cluster)] <- names(sequence)
            write.dna(cluster, paste0(subwd, files[j]), format = "fasta")
            species_included[[i]] <-
              data.frame(
                Included_species = names(sequence),
                AN = Acc_spp,
                Cluster = files[j]
              )
            cat("*******Sequence matched!******* \n")
            break ##If sequence is matched, we shold stop.
          } else {
            cat("****Still working \n")
          }

        }

        cat("**Done with sequence",
            i,
            "of",
            length(acce_new),
            "****\n")

      }
      df <- do.call("rbind", species_included)

      if (file.exists("sequence.fasta"))
        unlink("sequence.fasta")
      if (file.exists("test.txt"))
        unlink("test.txt")


    }

    setwd(mainDirect)
    cat("\n Retrieving accession numbers \n")

    filesna <-
      list.files(path = subwd,
                 pattern = c(".csv"),
                 full.names = T)
    filesna2 <- list.files(path = subwd, pattern = c(".csv"))

    everyGi <- lapply(filesna, read.csv)

    dat_acc2 <- list()
    for (m in 1:length(filesna)) {
      dat_acc <- ncbi_byid(gsub("gi", "", as.vector(everyGi[[m]]$gi)))
      dat_acc2[[m]] <-
        cbind.data.frame(taxa = dat_acc$taxon , dat_acc$acc_no)
    }
    merged.data.frame <-
      Reduce(function(...)
        merge(..., by = "taxa", all = TRUE), dat_acc2)
    names(merged.data.frame) <- c("Species", filesna2)
    write.csv(merged.data.frame,
              paste0("Molecular_sampling_", lineage, ".csv"))

    if (file.exists("sequence.fasta"))
      unlink("sequence.fasta")
    if (file.exists("test.txt"))
      unlink("test.txt")
    if (file.exists("Test.csv"))
      unlink("Test.csv")

    ##Checking taxonomy

    if (correct_db == T) {
      speciesbdb <-
        read.csv(paste0("Molecular_sampling_", lineage, ".csv"))

      del <- c()
      for (i in 1:dim(speciesbdb)[1]) {
        del[i] <-
          if (gnr_resolve(names = as.character(speciesbdb[i, 2]))[1, "score"] < 0.98)
            "delete"
        else
          "Ok"
        cat("Checking species", i, "of", dim(speciesbdb)[1], "\n")
      }

      sp_names <- strsplit(as.character(speciesbdb[, 2]), " ")

      del2 <- c()
      for (i in 1:length(sp_names)) {
        del2[i] <- if (length(sp_names[[i]]) > 2)
          "Delete"
        else
          "Ok"
      }


      spp_delete <-as.character(speciesbdb[, 2])[which(del == "Delete" |
                                                         del2 == "Delete")]

      ##Correct DB
      corrected.db <- speciesbdb[!speciesbdb$Species %in% spp_delete, ]
      ##Correct clusters
      for (i in 1:length(files)) {

        cluster <- read.dna(paste0(subwd, files[i]), format = "fasta")
        if(length(which(names(cluster) %in% gsub(" ", "_", spp_delete))) == 0){
          cluster}else{
            rep_la<-which(names(cluster) %in% gsub(" ",
                                                   "_", spp_delete))

            cluster<-cluster[-(rep_la)]
          }

        setwd(file.path(mainDir, "unaligned"))
        if (length(names(cluster)) == 0) {
          write.table(c("Cluster", i, "skipped"),
                      paste0("No_cluster_", files[i]))
          next} else{

            dp<- which(duplicated(names(cluster)))
            if(length(dp)==0){cluster}else{
              cluster<-cluster[-(dp)]
            }
            write.dna(cluster, paste0("corrected_", files[i]), format = "fasta")

          }

      }

    } else{
      cat("Done!")
    }



    setwd(mainDir)
    write.csv(df, "included.species.csv")
    write.csv(corrected.db, "phylota.sampling.csv")
    if (file.exists(grep("Molecular_*", list.files("."), value = T)))
      unlink(grep("Molecular_*", list.files("."), value = T))

    if (delete_all == T) {
      do.call(file.remove, list(setdiff(
        list(list.files(subwd, full.names = TRUE))[[1]], grep("corrected" , list(list.files(
          subwd, full.names = TRUE
        ))[[1]], value = T)
      )))
    } else{
      cat("Done")
    }

    ##

    cat("Creating a summary file of your molecular sampling")

    db1 <- read.csv("included.species.csv")
    db2 <- read.csv("phylota.sampling.csv")
    db3 <- db2[, -c(1:2)]

    for (i in 1:dim(db1)[1]) {
      a <-
        unlist(regmatches(
          as.character(db1[i, "Cluster"]),
          gregexpr("[[:digit:]]+", as.character(db1[i, "Cluster"]))
        ))
      b <-
        unlist(regmatches(names(db3)[-c(1)], gregexpr("[[:digit:]]+", names(db3)[-c(1)])))

      to_add <- NULL
      to_add <-  rep(NA, dim(db3)[2])
      to_add[1] <- as.character(db1$Included_species)[i]
      to_add[which(a == b) + 1] <- as.character(db1$AN)[i]

      db3 <- data.frame(rbind(as.matrix(db3), t(as.matrix(to_add))))
      print(i)
    }

    db3 <- db3[order(db3$Species),]
    db3<- db3[!duplicated(db3$Species),]
    db3<-melt(db3,id='Species',na.rm=TRUE)
    db3<-dcast(db3,Species~as.character(variable),value.var='value')

    write.csv(db3, paste0("Molecular_sampling_", lineage, ".csv"))
    #if (file.exists(fn3)) unlink(fn3,recursive =T)
    #if (file.exists(fn4)) unlink(fn4,recursive =T)

    cat("Please go over the molecular sampling file carefully")

    ##Do MSA--------
    ##Read fasta files
    if(MSA==TRUE){

      cat("\nalignment in process. Please be patient\n")

      setwd(file.path(mainDir, "Unaligned"))

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

        plot.progress <- function(percent) {
          plot(c(0,100), c(0,1), type='n', xlab='', ylab='', yaxt='n')
          rect(0, 0.1, percent*100, 0.9, col='blue')
          title(paste('Progress: ', round(percent*100,2), '%', sep=''))
        }


        aln_aliscored<-list()
        for (i in 1: length(alignments)){
          aln_aliscored[[i]]<-aliscore(alignments[[i]], gaps = "ambiguous", w = 3, path = paste0(mainDir, "/Aliscore_v.2.0"))
          write.dna(aln_aliscored[[i]], paste0(clade, "_", "AliScore", "_Alignment","_Cluster_" ,i, ".fasta"), format="fasta")
          plot.progress(i/length(alignments))
          cat("\nCurated alignments are under ALISCORE folder")

        }
        setwd(mainDirect)
      } else {}
      setwd(mainDirect)
    } else {}

    setwd(mainDir)
    setwd(mainDirect)
    ####Remove subspecies for the last time

    dfs <-
      lapply(list.files(
        pattern = c(".csv"),
        full.names = T
      ), read.csv)

    files_names<-list.files(
      pattern = c(".csv"),
      full.names = T
    )

    dfs_co<-lapply(1:length(dfs), function(x){
     vec<- gsub(" ", "_",as.character(dfs[[x]][,grep("species",colnames(dfs[[x]]),ignore.case=TRUE)]))
      num1<-which(sapply(strsplit(vec, "_"), length)==3)
      df_corr<- if(length(num1)==0){ dfs[[x]]}else{dfs[[x]][-num1,]}
      return(df_corr)
    })


    lapply(1:length(dfs_co), function(x){
      write.csv(dfs_co[[x]], file=files_names[x])
    })


    alns <-
      lapply(list.files(
        pattern = c(".fasta"),
        full.names = T, recursive = T
      ), read.dna, "fasta")

    names_alns <-list.files(
      pattern = c(".fasta"),
      full.names = T, recursive = T
    )


    lapply(1:length(alns), function(x){
      aln_ex<-alns[[x]]
      num1<-which(sapply(strsplit(labels(aln_ex), "_"), length)==3)
      if(length(num1)==0){ }else{
        newaln<- if(is.null(dim(aln_ex)) ==F){
          aln_ex[-(num1),]
        }else{
          aln_ex[-(num1)]
        }
        write.dna(newaln, names_alns[x], format = "fasta")
      }
    })

    ####Re organize datasets

    dfs <-
      lapply(list.files(
        pattern = c(".csv"),
        full.names = T
      ), read.csv)

    files_names<-list.files(
      pattern = c(".csv"),
      full.names = T
    )

    species_correct<- unique(gsub(" ", "_",as.character(dfs[[2]][,grep("species",colnames(dfs[[x]]),ignore.case=TRUE)])))

    dfs_co<-lapply(1:length(dfs), function(x){
      spps<- gsub(" ", "_",as.character(dfs[[x]][,grep("species",colnames(dfs[[x]]),ignore.case=TRUE)]))
      new_df<-dfs[[x]][!duplicated(spps),]
      spps<- gsub(" ", "_",as.character(new_df[,grep("species",colnames(new_df),ignore.case=TRUE)]))
      new_df<-new_df[which(spps %in% species_correct),]
      return(new_df)
    })


    lapply(1:length(dfs_co), function(x){
      write.csv(dfs_co[[x]], file=files_names[x])
    })


    alns <-
      lapply(list.files(
        pattern = c(".fasta"),
        full.names = T, recursive = T
      ), read.dna, "fasta")

    names_alns <-list.files(
      pattern = c(".fasta"),
      full.names = T, recursive = T
    )


    lapply(1:length(alns), function(x){
      aln_ex<-alns[[x]]
      num1<-which(sapply(strsplit(labels(aln_ex), "_"), length)==3)
      num2<-which(! labels(aln_ex) %in% species_correct)

      num3<-c(num1,num2)

      if(length(num3)==0){ }else{
        newaln<- if(is.null(dim(aln_ex)) ==F){
          aln_ex[-(num3),]
        }else{
          aln_ex[-(num3)]
        }
        write.dna(newaln, names_alns[x], format = "fasta")
      }
    })


    cat("\nDone!!")

    options(warn = 0)
  }
