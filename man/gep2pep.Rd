\name{gep2pep}
\alias{gep2pep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gep2pep(geps, pathw, parallel = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geps}{
%%     ~~Describe \code{geps} here~~
}
  \item{pathw}{
%%     ~~Describe \code{pathw} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
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
function (geps, pathw, parallel = F) 
{
    genemat <- geps
    genes <- rownames(genemat)
    x <- list()
    sets <- sapply(unname(pathw), get, x = "set")
    pb <- txtProgressBar()
    for (j in 1:ncol(genemat)) {
        setTxtProgressBar(pb, (j - 1)/ncol(genemat))
        genematj <- genemat[, j]
        "\%dobest\%" <- if (parallel) 
            get("\%dopar\%")
        else get("\%do\%")
        gres <- foreach(set = sets, .export = c("gsea", "ks.test.2")) \%dobest\% 
            {
                where <- match(set, genes)
                where <- where[!is.na(where)]
                gsea(where, genematj, F)
            }
        x[[j]] <- gres
    }
    setTxtProgressBar(pb, 1)
    close(pb)
    ESmat <- matrix(NA, length(pathw), ncol(genemat))
    PVmat <- matrix(NA, length(pathw), ncol(genemat))
    for (i in 1:ncol(genemat)) {
        PVmat[, i] <- sapply(x[[i]], "get", x = "p")
        ESmat[, i] <- sapply(x[[i]], "get", x = "ES")
    }
    rownames(ESmat) <- rownames(PVmat) <- sapply(pathw, "get", 
        x = "setid")
    colnames(ESmat) <- colnames(PVmat) <- colnames(genemat)
    return(list(ESmat = ESmat, PVmat = PVmat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
