\name{buildEmptyDB}
\alias{buildEmptyDB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
buildEmptyDB(path, pws, name = NULL, description = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
%%     ~~Describe \code{path} here~~
}
  \item{pws}{
%%     ~~Describe \code{pws} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{description}{
%%     ~~Describe \code{description} here~~
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
function (path, pws, name = NULL, description = NULL) 
{
    if (is.null(name)) 
        name <- "gep2pep database"
    if (is.null(description)) 
        description <- paste("This database contains pathway information", 
            "and pathway expression profiles created by", "the gep2pep package.")
    if (file.exists(path)) {
        say("Can not create DB in existing folder", T)
    }
    else rp <- repo_open(path, T)
    rp$project(name, description)
    db_ids <- buildDBids(pws)
    subdbs <- unique(db_ids)
    for (i in 1:length(subdbs)) {
        dbi <- subdbs[i]
        say(paste("Storing pathway data for DB:", dbi))
        rp$put(pws[db_ids == dbi], paste0(subdbs[i], "_pws"), 
            paste("Pathway information for DB", subdbs[i]), c("gep2pep", 
                "pws"))
    }
    rp$put(unique(db_ids), "DB list", "IDs of pathway sub-databases included in this repository", 
        c("gep2pep", "meta"))
    return(rp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
