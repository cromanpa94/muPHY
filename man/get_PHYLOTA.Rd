\name{get_PHYLOTA}
\alias{get_PHYLOTA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Download PHYLOTA clusters
%%  ~~function to do ... ~~
}
\description{A function fo retrieving PHYLOTA clusters for a given clade. Make sure that there are
informative clusters in the website. Otherwise it will not work. This function (1) creates a new folder in
your working directory, and (2) exports each cluster in fasta format.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_PHYLOTA(clade,MSA=TRUE, ALI=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clade}{Character: A given lineage (e.g. "Trochidae")}
  \item{MSA}{Logical: Perform a MSA on the clusters (only ClustalOmega)}
  \item{ALI}{Logical: Use Aliscore on the MSA}
%%     ~~Describe \code{x} here~~
}

\details{Use a single clade.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
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
\note{Make sure you have the lastest version of R
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
get_PHYLOTA("Trochilidae", MSA=TRUE, ALI=TRUE)
}
