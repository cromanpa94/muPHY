\name{Phylo_Get_Div}
\alias{Phylo_Get_Div}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Phylo_Get_Div
%%  ~~function to do ... ~~
}
\description{Estimate the ME diversification rates for a given tree
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Phylo_Get_Div(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{ A tree
%%     ~~Describe \code{x} here~~
}

  \item{richness}{ A dataframe containing two columns. The first are the genera names (this could be also changed for families, orders, etc), and a second column with the species richness.

  Make sure that the taxonomy in the first column matches with the tree labels
%%     ~~Describe \code{x} here~~
}


}
\details{
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
