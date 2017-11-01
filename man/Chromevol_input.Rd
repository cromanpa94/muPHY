\name{create_chromevol_unput}
\alias{create_chromevol_unput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{create_chromevol_unput
%%  ~~function to do ... ~~
}
\description{Creates ChromEvol input files for a set of nodes.

This function creates a single input file for each subtree. Useful when analyzing big trees. However, you can only input the root number to obtain a single chromevol input for the entire tree.


Make sure that PIP_control is in your working directory. Note that it doesn't have any extention! Download my file from the following link:

https://www.dropbox.com/s/ny1jvkw3diufoyp/PIP_control?dl=0

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create_chromevol_unput(nodes, tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodes}{A vector containing the nodes of interest
%%     ~~Describe \code{x} here~~
}
  \item{tree}{Your tree name
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
