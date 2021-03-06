% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gep2pep.R
\name{getResults}
\alias{getResults}
\title{Extracts the results matrix from \code{CondSEA} or \code{PathSEA}
output}
\usage{
getResults(analysis, collection)
}
\arguments{
\item{analysis}{The output of either \code{CondSEA} or
\code{PathSEA}.}

\item{collection}{One of the names returned by
\code{getCollections}.}
}
\value{
A 2-columns matrix including ESs and
    p-values (see details) for each pathway database and
    condition.
}
\description{
Extracts the results matrix from \code{CondSEA} or \code{PathSEA}
output
}
\examples{
db <- loadSamplePWS()
repo_path <- file.path(tempdir(), "gep2pepTemp")

rp <- createRepository(repo_path, db)
geps <- loadSampleGEP()
buildPEPs(rp, geps)

pgset <- c("(+)_chelidonine", "(+/_)_catechin")
psea <- CondSEA(rp, pgset)

getResults(psea, "c3_TFT")

unlink(repo_path, TRUE)

}
\seealso{
CondSEA, PathSEA
}
