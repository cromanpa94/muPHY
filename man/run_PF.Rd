\name{run_PF}
\alias{run_PF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run PartitionFinder 1 from R
%%  ~~function to do ... ~~
}
\description{Need to set the name of the folder within your Working directory that contains a set of fasta files (i.e., only MSA are allowed).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
run_PF(folder = "ALI")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder}{ALI: Folder ALISCORE; ali: Alignment folder; Any folder within your working directory.}
%%     ~~Describe \code{folder} here~~

}

\value{
\item{folder}{"ALI"}
}

\references{
Lanfear, R., Calcott, B., Ho, S. Y., & Guindon, S. (2012). PartitionFinder: combined selection of partitioning schemes and substitution models for phylogenetic analyses. Molecular biology and evolution, 29(6), 1695-1701.
}
\author{
Cristian Roman-Palacios}
\examples{
\dontrun{

 run_PF(folder = "ALI")
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~PartitionFinder }% use one of  RShowDoc("KEYWORDS")
