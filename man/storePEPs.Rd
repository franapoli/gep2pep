\name{storePEPs}
\alias{storePEPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
storePEPs(rp, db_id, peps, existing)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{db_id}{
%%     ~~Describe \code{db_id} here~~
}
  \item{peps}{
%%     ~~Describe \code{peps} here~~
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
function (rp, db_id, peps, existing) 
{
    okargs <- c("overwrite", "append", "stop")
    if (length(existing) > 1 || (!existing \%in\% okargs)) 
        say(paste0("existing must be one of: ", paste(okargs, 
            collapse = ", ")), T)
    replace <- F
    if (rp$has(db_id)) {
        if (existing == "stop") 
            say("PEP already exists", T)
        if (existing == "append") {
            say("Merging peps...")
            curpep <- rp$get(db_id)
            peps[["ESmat"]] <- cbind(curpep[["ESmat"]], peps[["ESmat"]])
            peps[["PVmat"]] <- cbind(curpep[["PVmat"]], peps[["PVmat"]])
            replace <- T
        }
        if (existing == "overwrite") 
            replace <- T
    }
    say("Storing pathway expression profiles...")
    rp$put(peps, db_id, paste0("Pathway data for DB ", db_id, 
        ". It contains 2 matrices: 1 for enrichement scores ", 
        "(signed Kolmogorov Smirnov statistic) and one for ", 
        "the corresponding p values."), c("gep2pep", "pep"), 
        replace = replace)
    curids <- rp$get("DB list")
    rp$put(colnames(peps[[1]]), "perturbagens", "Names of perturbagens inducing expression profiles", 
        c("gep2pep", "meta"), replace = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
