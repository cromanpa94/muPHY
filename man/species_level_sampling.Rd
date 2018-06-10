\name{species_level_sampling}
\alias{species_level_sampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Obtain species-level sequence fasta files.
%%  ~~function to do ... ~~
}
\description{This function reduce sequence files in fasta format to species-level samples (one sample per species per cluster). It runs on the current working directoy and creates a subdirectory with the corrected files.

}
\usage{
species_level_sampling(pattern)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pattern}{Common pattern for every target sequence files (e.g. Cluster)}


%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\references{
%% ~put references to the literature/web site here ~
}
\author{Cristian Roman-Palacios
%%  ~~who you are~~
}
\examples{
\dontrun{
species_level_sampling("Cluster")
}
}
