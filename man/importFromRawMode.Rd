% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gep2pep.R
\name{importFromRawMode}
\alias{importFromRawMode}
\title{Imports PEPs created in raw mode}
\usage{
importFromRawMode(rp, path = file.path(rp$root(), "raw"),
  collections = "all")
}
\arguments{
\item{rp}{A repository created by \code{\link{createRepository}}.}

\item{path}{Path were raw PEPs are stored (default is a "raw"
directory under the repository root folder).}

\item{collections}{A subset of the collection names returned by
\code{getCollections}. If set to "all" (default), all the
collections in \code{rp} will be used.}
}
\value{
Nothing, used for side effects.
}
\description{
Raw mode is meant to deal with large collections of PEPs (like
hundreds of thousands). In this case, problems may arise while
trying to convert GEPs by loading all of them in memory at
once. Raw mode is meant to be used with HDF5 format, which allows
to load subsets of GEPs from the disk. \code{buildPEPs}, when used
in raw mode, can create the corresponding subsets of PEPs, so that
the job can be distributed on a computer
cluster. \code{importFromRawMode} is meant to join the chunks into
HDF5 matrices, which are than stored into the repository. The
\code{.loadPEPs} function can seamlessly load PEPs stored in normal
(RDS) or HDF5 format.
}
\details{
PEPs are expect to be found at the specified \code{path}
    and follow the naming convention as generated by
    \code{buildPEPs}. According to such convention, each file is
    named usign the format
    category_subcategory#chunknumber.RDS. All non-alphanumeric
    characters from the original category and subcategory names are
    replace with an underscore (in rare cases this could create
    ambiguity that should be manually prevented). All chunks for
    the same subcategory are joined together following the chunk
    numbers into a single HDF5 matrix and stored in the repository
    as an "attachment" (see \code{repo} documentation).

    Note that raw PEPs (by default everything at
    repository_root/raw) can be safly removed once they have been
    imported.
}
\examples{
db <- loadSamplePWS()
repo_path <- file.path(tempdir(), "gep2pepTemp")

rp <- createRepository(repo_path, db)

## The following will create PEPs in 2 separate files
geps <- loadSampleGEP()
buildPEPs(rp, geps[,1:2], progress_bar=FALSE,
    rawmode_id=1)
buildPEPs(rp, geps[,3:5], progress_bar=FALSE,
    rawmode_id=2)

## The separate files are then merged into one (possibly big) file
## in HDF5 format

importFromRawMode(rp)

## Now most operations (excluding the addition of new PEPs to
## existing collections) will be available as usual.

unlink(repo_path, TRUE)
}
\seealso{
buildPEPs
}
