\name{get_ortho_clusters}
\alias{get_ortho_clusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Obtain orthologous clusters for a given lineage (and outgroup taxa)
%%  ~~function to do ... ~~
}
\description{This function retrieves all sequences from genbank for a given lineages. Next, it performs blast all-vs-all blast and summarizes all orthologous clusters. It is also able to align

}
\usage{
get_ortho_clusters(ingroup, outgroup, T)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ingroup}{String: Ingroup name (e.g. "Sphyraena")}
  \item{Outgroup}{String: Outgroup name}
  \item{MSA}{Logical: align sequences?}


%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\references{

Sanderson, M. J., Boss, D., Chen, D., Cranston, K. A., & Wehe, A. (2008). The PhyLoTA Browser: processing GenBank for molecular phylogenetics research. Systematic Biology, 57(3), 335-346.

Sievers, F., & Higgins, D. G. (2014). Clustal Omega, accurate alignment of very large numbers of sequences. Multiple sequence alignment methods, 105-116.
%% ~put references to the literature/web site here ~
}
\author{Cristian Roman-Palacios
%%  ~~who you are~~
}
\examples{
\dontrun{
get_ortho_clusters("Centrolene", c("Hyalinobatrachium", "Cochranella"), T)
 }
}
