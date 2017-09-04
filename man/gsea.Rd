\name{gsea}
\alias{gsea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsea(S, ranks_list, check = F, alternative = "two.sided", leadedge = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{ranks_list}{
%%     ~~Describe \code{ranks_list} here~~
}
  \item{check}{
%%     ~~Describe \code{check} here~~
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
}
  \item{leadedge}{
%%     ~~Describe \code{leadedge} here~~
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
function (S, ranks_list, check = F, alternative = "two.sided", 
    leadedge = F) 
{
    S <- S[!(is.na(S))]
    S1 <- ranks_list[S]
    S2 <- ranks_list[-S]
    if (length(S1) < 1 || length(S2) < 1 || all(is.na(S1)) || 
        all(is.na(S2))) 
        return(list(ES = NA, p = NA))
    ks <- ks.test.2(S1, S2, alternative = alternative, maxCombSize = 10^10)
    lead = NA
    if (leadedge) {
        sSm <- sort(Sm)
        if (ks$ES > 0) 
            lead = sSm[sSm <= ks$edge]
        else lead = sSm[sSm >= ks$edge]
    }
    return(list(ES = ks$ES, p = ks$p.value, edge = ks$edge, lead = lead))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
