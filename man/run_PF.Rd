\name{run_PF}
\alias{run_PF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run PartitionFinder 1 from R
%%  ~~function to do ... ~~
}
\description{Only set the name of the folder within your Working directory that contains a set of fasta files (msa) and run the code. You can also use the output from get_PHYLOTA function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
run_PF(folder = "WD")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder}{ALI: Folder ALISCORE; ali: Alignment folder; Any folder within your working directory.}
%%     ~~Describe \code{folder} here~~

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{folder}{"ALI"}
}

\references{
Lanfear, R., Calcott, B., Ho, S. Y., & Guindon, S. (2012). PartitionFinder: combined selection of partitioning schemes and substitution models for phylogenetic analyses. Molecular biology and evolution, 29(6), 1695-1701.
}
\author{
Cristian Roman-Palacios}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
run_PF(folder = "WD")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~PartitionFinder }% use one of  RShowDoc("KEYWORDS")
