\name{update_phylota}
\alias{update_phylota}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update phylota clusters
%%  ~~function to do ... ~~
}
\description{This function updates phylota clusters using sequences from genbank. It overwrites the clusters that are download. I've used it only on mac.

******Make sure to have blast installed in your computer (also in your path)*******

For windows users: Follow the steps outlined in the following link.

http://dmcase.blogspot.com.br/2014/02/how-to-install-blast-on-windows.html

Don't forget to set the path!


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update_phylota(lineage, nsamples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lineage}{Character: A given lineage (e.g. "Sphyraena")}
  \item{nsamples}{numeric: Number of species sampled from clusters and used to obtain gene symbols}
%%     ~~Describe \code{x} here~~
}

\details{Use a single clade.
%%  ~~ If necessary, more details than the description above ~~
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
\note{Make sure you have the lastest version of R
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
get_PHYLOTA("Sphyraena", MSA=TRUE, ALI=TRUE)
}
