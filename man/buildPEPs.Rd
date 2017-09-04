\name{buildPEPs}
\alias{buildPEPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
buildPEPs(rp, geps, parallel = F, existing = "stop")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{geps}{
%%     ~~Describe \code{geps} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{existing}{
%%     ~~Describe \code{existing} here~~
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
function (rp, geps, parallel = F, existing = "stop") 
{
    okargs <- c("stop", "overwrite", "skip", "append")
    if (length(existing) > 1 || (!existing \%in\% okargs)) 
        say(paste0("existing must be one of: ", paste(okargs, 
            collapse = ", ")), T)
    dbs <- rp$get("DB list")
    for (i in 1:length(dbs)) {
        say(paste0("Working on DB: ", dbs[i], " (", i, "/", length(dbs), 
            ")"))
        if (rp$has(dbs[i])) 
            if (existing != "stop") {
                say(paste0("Existing DB found, action: ", existing))
            }
            else say(paste0("DB already present, check parameter: existing"), 
                T)
        if (rp$has(dbs[i]) && existing == "skip") 
            next
        thisdb <- rp$get(paste0(dbs[i], "_pws"))
        peps <- gep2pep(geps, thisdb, parallel)
        storePEPs(rp, dbs[i], peps, existing)
        cat("\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
