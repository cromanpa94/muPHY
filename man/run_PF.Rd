\name{run_PF}
\alias{run_PF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Use PartitionFinder 1 from R
%%  ~~function to do ... ~~
}
\description{PF R
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
run_PF(folder = "WD")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder}{
%%     ~~Describe \code{folder} here~~
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
function (folder = "WD")
{
    mainDir <- getwd()
    subDir <- if (folder == "WD") {
        paste0(mainDir, "/Aliscore")
    }
    else {
        paste0(mainDir, "/Aligned")
    }
    setwd(subDir)
    evobiR::SuperMatrix()
    concatenated <- read.FASTA("concatenated.fasta")
    setwd(mainDir)
    download.file("https://github.com/brettc/partitionfinder/archive/v1.1.1.zip",
        "PartitionFinder.tar.gz")
    untar("PartitionFinder.tar.gz")
    partitions <- read.csv(paste0(subDir, "/concatenated.partitions.csv"))
    config_file <- readLines("partitionfinder-1.1.1/examples/nucleotide/partition_finder.cfg")
    block <- list()
    for (i in 1:dim(partitions)[1]) {
        block[[i]] <- paste0("GENE_", i, " ", "=", " ", partitions[i,
            2], " ", "-", " ", partitions[i, 3], ";")
    }
    conf_fi_mod <- config_file[-c(16:24)]
    newconf <- append(conf_fi_mod, unlist(block), after = 15)
    newconf[9] <- paste("models =", "beast;")
    newconf[2] <- paste("alignment =", "concatenated.phy;")
    subDir_4 <- "PartitionFinder_Analyses"
    dir.create(file.path(mainDir, subDir_4))
    setwd(file.path(mainDir, subDir_4))
    writeLines(newconf, "partition_finder.cfg")
    write.dna(concatenated, file = "concatenated.phy", format = "sequential",
        nbcol = -1, indent = 0, colsep = "", blocksep = 0)
    setwd(mainDir)
    tofile <- paste0(getwd(), "/PartitionFinder_Analyses/partition_finder.cfg")
    system(paste0("python ./partitionfinder-1.1.1/PartitionFinder.py",
        " ", tofile), wait = F)
    cat("You can use R while PartitionFinder is running!\n\n Results are written to PartitionFinder_Analyses folder")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
