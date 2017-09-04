\name{importMsigDB.xml}
\alias{importMsigDB.xml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
importMsigDB.xml(fname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
%%     ~~Describe \code{fname} here~~
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
function (fname) 
{
    library(XML)
    xml <- xmlTreeParse(fname, useInternalNodes = T)
    sets <- xml["/MSIGDB/GENESET"]
    ids <- sapply(sets, function(x) xmlAttrs(x)[["SYSTEMATIC_NAME"]])
    msigDB <- data.frame(setid = ids, setname = sapply(sets, 
        function(x) xmlAttrs(x)[["STANDARD_NAME"]]), db = sapply(sets, 
        function(x) xmlAttrs(x)[["CATEGORY_CODE"]]), subdb = sapply(sets, 
        function(x) xmlAttrs(x)[["SUB_CATEGORY_CODE"]]), organism = sapply(sets, 
        function(x) xmlAttrs(x)[["ORGANISM"]]), desc = sapply(sets, 
        function(x) xmlAttrs(x)[["DESCRIPTION_BRIEF"]]), desc_full = sapply(sets, 
        function(x) xmlAttrs(x)[["DESCRIPTION_FULL"]]), set = sapply(sets, 
        function(x) xmlAttrs(x)[["MEMBERS_SYMBOLIZED"]]), stringsAsFactors = F)
    msigDB <- apply(msigDB, 1, as.list)
    msigDB <- lapply(msigDB, function(x) {
        x$set <- strsplit(x$set, ",")[[1]]
        x
    })
    names(msigDB) <- ids
    return(msigDB)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
