\name{supermatrix_SF}
\alias{supermatrix_SF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Concatenate files based on a given threshold (missing data)
%%  ~~function to do ... ~~
}
\description{This concatenates fasta files based on a given threshold. It read all fasta files within the working directory that matches matching a given pattern (e.g. Cluster)

}
\usage{
supermatrix_SF(threshold, pattern)

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
supermatrix_SF(threshold=0.2, pattern="Cluster")
}
}
