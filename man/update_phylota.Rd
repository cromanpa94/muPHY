\name{update_phylota}
\alias{update_phylota}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update Phylota Clusters
%%  ~~function to do ... ~~
}
\description{This function updates a set of phylota clusters using sequences from genbank.

Make sure to have blast installed in your computer (also in your path). For windows users follow the steps outlined in the following link:

https://dmcase.blogspot.com.br/2014/02/how-to-install-blast-on-windows.html

Instructions for MAC users:

https://gist.github.com/knmkr/5393474

In short, this function first download all the available orthologous clusters for a given lineage. Then, it looks for non-overlapping species between the sampled (i.e., PhyloTa) and newly published sequences (i.e., genbank). This function also allows the inclussion of outgroup species.

Each new sequence is then blastes agains each orthologous cluster and merged when a significant match is found. Next, taxonomy is corrected and sampling summaries are exported.

Additional arguments are provided for users interested in performing multiple sequence alignment or using aliscore.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update_phylota(lineage,nsamples,genes,MSA,ALI,outgroup,correct_db,delete_all,c_directory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lineage}{String: A given lineage (e.g. "Sphyraena")}
  \item{nsamples}{numeric: If sequence names are not provided, type the number of sequences to be used when assesing the sequence identity (i.e., more sequence will take more time to process, but the likelihood of ID is higher. I suggest using at least 5)}
  \item{genes}{String: Which genes should the pipeline look for? Use either previous studies to guide yourself or look at PhyloTa database}
  \item{MSA}{Logical: Align sequences?}
  \item{ALI}{Logical: Run Aliscore (only if MSA=T)?}
  \item{outgroup}{String: A vector of strings containing the names of outgroup species}
  \item{correct_db}{Logical: Whether an additional correction step is conducted after the new sequences are included. Subspecies and sp. are removed from the resulting files}
  \item{delete_all}{Logical: Delete all intermediate files}
  \item{c_directory}{Logical: Should the function work in your current working directory?}

%%     ~~Describe \code{x} here~~
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{df}{A dataframe containing all the species that were included}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Sanderson, M. J., Boss, D., Chen, D., Cranston, K. A., & Wehe, A. (2008). The PhyLoTA Browser: processing GenBank for molecular phylogenetics research. Systematic Biology, 57(3), 335-346.
%% ~put references to the literature/web site here ~
}
\author{Cristian Roman-Palacios
%%  ~~who you are~~
}

\examples{
\dontrun{
update_phylota("Alcedinidae",genes=c("RAG-1","ND2"),"Actenoides hombroni",ALI=T,MSA=T,c_directory=T)
}
}
