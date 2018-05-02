\name{get_PHYLOTA}
\alias{get_PHYLOTA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Download PHYLOTA clusters
%%  ~~function to do ... ~~
}
\description{This function retrieves all PHYLOTA clusters for a given clade. Make sure that there are
informative clusters in the database. Otherwise it will not work. This function (1) creates a new folder in
your working directory, and (2) saves each cluster into a fasta format file.

Two additional arguments were included into the function to (1) perform multiple sequence alignment of each orthologous cluster and (2) run aliscore on the same set.

I suggest using SuperMatrix function (evobioR package) to concatenate MSAs.

}
\usage{
get_PHYLOTA(clade,MSA=FALSE, ALI=FALSE, c_directory=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clade}{String: Lineage name (e.g. "Sphyraena")}
  \item{MSA}{Logical: Perform a MSA on the clusters (only ClustalOmega available)}
  \item{ALI}{Logical: Run Aliscore on the MSA}
  \item{c_directory}{Logical: Should the function work in your current working directory?}

%%     ~~Describe \code{x} here~~
}

\details{Please select a single clade. No
%%  ~~ If necessary, more details than the description above ~~
}

\references{

Sanderson, M. J., Boss, D., Chen, D., Cranston, K. A., & Wehe, A. (2008). The PhyLoTA Browser: processing GenBank for molecular phylogenetics research. Systematic Biology, 57(3), 335-346.

Misof, B., & Misof, K. (2009). A Monte Carlo approach successfully identifies randomness in multiple sequence alignments: a more objective means of data exclusion. Systematic biology, 58(1), 21-34.

Sievers, F., & Higgins, D. G. (2014). Clustal Omega, accurate alignment of very large numbers of sequences. Multiple sequence alignment methods, 105-116.
%% ~put references to the literature/web site here ~
}
\author{Cristian Roman-Palacios
%%  ~~who you are~~
}
\examples{
\dontrun{
 get_PHYLOTA("Sphyraena", MSA=TRUE, ALI=TRUE, c_directory=T)
 }
}
