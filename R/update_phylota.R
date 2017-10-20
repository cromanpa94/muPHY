update_phylota<-function(lineage, nsamples=5, database="ncbi"){
  fn <- "Unaligned"
  if (file.exists(fn)) unlink(fn,recursive =T)

  cat("Get PhyLota clusters")

  get_PHYLOTA(lineage,MSA=F, ALI=F)
  subwd<-paste0(getwd(), "/Unaligned/")

  ##Check for new species in genbank

  cat("Check novel species in genbank")

  spp_genbank<-downstream(lineage, db = database, downto = 'species')[[1]]
  spp_sampled <- do.call(rbind,lapply(list.files(path = subwd, pattern = c( "csv"), full.names=T),read.csv))
  spp_sampled <-spp_sampled[!duplicated(spp_sampled$name),]
  new_spp<-setdiff(gsub(" ", "_", spp_genbank$childtaxa_name), spp_sampled$name)

  if(length(new_spp) ==0){ cat("Nothing to add...sorry")
    if (file.exists(fn)) unlink(fn,recursive =T)
        }else{

    cat("Tranform each cluster into a Blast DB")

    ##First make DB for each cluster
    files<-list.files(path = subwd, pattern = c( "Cluster"))
    files<-Filter(function(x) grepl("\\.fasta$", x), files)


    for(i in 1:length(files)){
      args<-paste0("-in ", subwd, files[i], " -input_type fasta -dbtype nucl")
      system2(command = 'makeblastdb', args = args, wait=FALSE,stdout = TRUE)
      Sys.sleep(1)
      print(i)
    }



    cat(length(new_spp)," species can be added to the ", dim(spp_sampled)[1], "sampled in PhyloTa")

    ##Check which genes do I have in PhyloTa clusters. This is based on 5 species per each cluster.

    sample_Gi <- do.call(rbind,lapply(list.files(path = subwd, pattern = c( "csv"), full.names=T),read.csv, nrow=nsamples))
    sampled_genes<-ncbi_byid(gsub("gi","",as.vector(sample_Gi$gi)))

    cat("Lets look for the genes that were sampled in PhyLota")


    choosebank("genbank")
    gene_fin<-list()
    for(i in 1:length(sampled_genes$acc_no)){
      Dengue1<-query("Dengue1", paste0("AC=",gsub("\\..*","",sampled_genes$acc_no[i])))
      annots <- getAnnot(Dengue1$req[[1]])
      a<-gsub("/gene=\"","", grep("/gene=",annots, value = T))
      b<-gsub("\"","",a)
      gene_fin[[i]] <-gsub("[[:space:]]", "", b)[1]
      cat("New gene symbol found!", gsub("[[:space:]]", "", b)[1])
      print(i)
    }

    closebank()

    cat("Your sampling includes", levels(factor(unlist(gene_fin))), "loci")

    sampled_gene_names<- levels(factor(unlist(gene_fin)))

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
    first.of.type <- function (results, type) {
      # Filter rows whose title contains the given text
      filtered <- results[grepl(type, results$Title, fixed=TRUE),]

      # Return the first accession number
      filtered$OSLT[1]
    }
    scrape.genbank <- function (species, genes) {
      # Create dataset skeleton
      data <- data.frame(Species=species)
      for (i in genes) {
        data[,i] <- rep(NA, length(species))
      }

      # Look up data for each species
      for (i in 1:length(species)) {
        n <- species[i]
        print(sprintf("Looking up %s (%d/%d)...", n, i, length(species)))
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
    data <- scrape.genbank(gsub("_", " ",new_spp), sampled_gene_names)

    ###Now, I download each sequence and test where it fits better

    cat("Fitting the new sequences into the existing clusters")

    acce_new<-as.character(na.omit(unlist(data[,-1])))


    n_clust<-length(list.files(path = subwd, pattern = c( "csv")))
    species_included<-list()
    for (i in 1:length(acce_new)){

      sequence <- read.GenBank(acce_new[i])
      Acc_spp<-names(sequence)
      names(sequence)<-attr(sequence, "species")

      for(j in 1:n_clust){
        write.dna(sequence, "sequence.fasta", format="fasta")
        que<-paste0("-outfmt 6 -query sequence.fasta -out test.txt -db ", subwd,files[j])
        system2(command = 'blastn', args = que, wait=FALSE,stdout = TRUE)

        if(file.info("test.txt")$size>0){
          cluster<- read.dna( paste0(subwd,files[j]), format = "fasta")
          cluster[length(cluster)+1]<- sequence
          names(cluster)[length(cluster)]<-names(sequence)
          write.dna(cluster, paste0(subwd,files[j]), format="fasta")
          species_included[[i]] <- data.frame(Included_species=names(sequence), AN=Acc_spp,Cluster=files[j])
          cat("*******Sequence matched!******* \n")
          break ##If sequence is matched, we shold stop.
        }else{
          cat("****Still working \n")
        }

      }

      cat("**Done with sequence", i, "of", length(acce_new), "****\n")

    }
    df <- do.call("rbind", species_included)

    return(df)
  }
}
