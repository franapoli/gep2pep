\name{rankPEPsByCols}
\alias{rankPEPsByCols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rankPEPsByCols(peps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peps}{
%%     ~~Describe \code{peps} here~~
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
function (peps) 
{
    rankPEP <- function(PVs, ESs) {
        sorter <- abs(PVs)
        pos <- ESs > 0
        pos[is.na(pos)] <- F
        sorter[pos] = -sorter[pos]
        return(rank(sorter, ties.method = "random", na.last = "keep"))
    }
    PVs <- peps[["PVmat"]]
    ESs <- peps[["ESmat"]]
    x <- sapply(1:ncol(PVs), function(i) rankPEP(PVs[, i], ESs[, 
        i]))
    colnames(x) <- colnames(PVs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
