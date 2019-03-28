#' gep2pep: creation and analysis of Pathway Expression Profiles
#'
#' Pathway Expression Profiles (PEPs) are based on the expression of
#' pathways (or generic gene sets) belonging to a collection, as
#' opposed to individual genes. \code{gep2pep} supports the conversion
#' of gene expression profiles (GEPs) to PEPs and performs enrichment
#' analysis of both pathways and conditions.
#'
#' @details
#'
#' \code{gep2pep} creates a local repository of gene sets, which can
#' also be imported from the MSigDB [1] database. The local repository
#' is in the \code{repo} format. When a GEP, defined as a ranked list
#' of genes, is passed to \code{\link{buildPEPs}}, the stored database
#' of pathways is used to convert the GEP to a PEP and permanently
#' store the latter.
#'
#' One type of analysis that can be performed on PEPs and that is
#' directly supported by \code{gep2pep} is the Drug-Set Enrichment
#' Analysis (DSEA [2]). It finds pathways that are consistently
#' dysregulated by a set of drugs, as opposed to a background of other
#' drugs. Of course PEPs may refer to non-pharmacological conditions
#' (genetic perturbations, disease states, cell types, etc.) for
#' analogous analyses. See \code{\link{CondSEA}} function.
#'
#' A complementary approach is that of finding conditions that
#' consistently dysregulate a set of pathways. This is the
#' pathway-based version of the Gene Set Enrichment Analysis
#' (GSEA). As an application example, this approach can be used to
#' find drugs mimicking the dysregulation of a gene by looking for
#' drugs dysregulating pathways involving the gene (this has been
#' published as the \code{gene2drug} tool [3]). See
#' \code{\link{PathSEA}}.
#'
#' Both DSEA and gene2drug analyses can be performed using
#' preprocessed data from
#' \url{http://dsea.tigem.it/downloads.php}. The data include
#' Connectivity Map [4] GEPs (drug-induced gene expression profiles)
#' converted to PEPs in the form of a \code{gep2pep} repository.
#' 
#' Naming conventions:
#'
#' \itemize{
#'
#'   \item{pathway: }{any set of gene identifiers (not necessarily
#'   representing a molecular pathway).}
#'
#'   \item{pathway collection: }{a set of pathways.}
#'
#'   \item{pathway database: }{a set of pathway collections, like the
#'   MSigDB.}
#'
#'   \item{Gene Expression Profile (GEP): }{a named vector where names
#'   are gene identifiers of the same type as those in the pathway
#'   database and elements are ranks ranging from 1 to the number of
#'   genes.}
#'
#'   \item{Pathway Expression Profile (PEP): }{a ranked list of
#'   pathways, as converted from a GEP according to a pathway
#'   collection.}
#'
#'   \item{condition: }{any transcriptomic-modelled biological state
#'   (drug treatment, gene knock-out, disease state, cell type, etc.)
#'   characterized by an induced GEP and therefore a PEP.}
#'
#'   \item{gep2pep repository: }{a pathway database and possibly a
#'   related database of PEPs as created by the \code{gep2pep}
#'   package. It is implemented in \code{repo} format.}
#' }
#' @references 
#' [1] Subramanian A. et al. Gene set enrichment analysis: A
#'     knowledge-based approach for interpreting genome-wide
#'     expression profiles. PNAS 102, 15545-15550 (2005).
#'
#' [2] Napolitano F. et al, Drug-set enrichment analysis: a novel tool
#'     to investigate drug mode of action. Bioinformatics 32, 235-241
#'     (2016).
#' 
#' [3] Napolitano, F. et al. gene2drug: a computational tool for
#'     pathway-based rational drug repositioning. Bioinformatics
#'     (2017). https://doi.org/10.1093/bioinformatics/btx800
#' 
#' [4] Lamb, J. et al. The Connectivity Map: Using Gene-Expression
#'     Signatures to Connect Small Molecules, Genes, and Disease. Science
#'     313, 1929-1935 (2006).
#' @docType package
#' @name gep2pep-package
#' @author Francesco Napolitano \email{franapoli@@gmail.com}
#' @aliases gep2pep
NULL

## repo is for storage of repositories
#' @import repo

## hdf5 is for raw-mode (large matrices)
#' @import rhdf5

## For gene sets management
#' @import GSEABase

## digest is needed to cache merged profiles
#' @importFrom digest digest

## iter is used with foreach
#' @importFrom iterators iter

#' @importFrom foreach foreach %dopar% %do%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats ks.test pchisq
#' @importFrom XML xmlTreeParse xmlAttrs
#' @importFrom Biobase mkScalar
#' @importFrom methods is new

NULL

#' A class to contain categorized gene set collection
#'
#' This class is a simple generalization of the
#' \code{BroadCollection} function of \code{GSEABase} to store gene
#' sets having assigned categories and subcategories that can be
#' different from those of the MSigDB.
#'
#' @slot category A character defining the main category that the gene
#'     set belongs to.
#' @slot subCategory A character defining the secondary category that
#'     the gene set belongs to.
#' @name CategorizedCollection-class
#' @rdname CategorizedCollection-class
#' @export
setClass("CategorizedCollection",
         contains = "CollectionType",
         representation = representation(
             category = "character",
             subCategory = "character"),         
         prototype = prototype(
             type = mkScalar("Categorized"),
             category = "uncategorized",
             subCategory = "uncategorized"),
         )

#' Constructor method for objects of class CategorizedCollection.
#'
#' See \code{CategorizedCollection-class}.
#'
#' @param category A character defining the main category that the
#'     gene set belongs to.
#' @param subCategory A character defining the secondary category that
#'     the gene set belongs to.
#' @return An object of class \code{CategorizedCollection}.
#' @export
#' @examples
#'
#' library(GSEABase)
#' gs1 <- GeneSet(setName="set1", setIdentifier="101")
#' collectionType(gs1) <- CategorizedCollection()
CategorizedCollection <- function(category="uncategorized",
                                  subCategory="uncategorized") {
    if (length(category)!=1)
        stop("category must be a scalar (length = 1)")
    
    new("CategorizedCollection",
        category=category,
        subCategory=subCategory)
}

#' Loads sample Gene Expression Profiles
#'
#' @return Sample gene expression data
#' @export
#' @examples
#'
#' geps <- loadSampleGEP()
#' 
loadSampleGEP <- function() {
    return(
        readRDS(system.file("extdata", "testgep.RDS",
                            package="gep2pep"))
    )
}

#' Loads sample pathway collections
#'
#' @return Sample pathway collections in \code{GeneSetCollection}
#'     format
#' @export
#' @examples
#'
#' geps <- loadSampleGEP()
#'
loadSamplePWS <- function() {
  db <- readRDS(system.file("extdata", "testgmd.RDS",
                            package="gep2pep"))
  db <- as.CategorizedCollection(db)
  return(db)
}


setMethod("show",
          signature=signature(object="CategorizedCollection"),
          function(object) {
              cat("collectionType: ", collectionType(object), "\n",
                  "  category: ",
                  attr(object, "category"), "\n",
                  "  subCategory: ",
                  attr(object, "subCategory") ,"\n", sep="")
          })


#' Converts GeneSetCollection objects to CategorizedCollection objects.
#'
#' @param GScollection An object of class \code{GeneSetCollection}.
#' @param category The name of the category that all the gene sets
#'     will be assigned to (see details).
#' @param subCategory The name of the subcategory that all the gene
#'     sets will be assigned to (see details).
#' @return A CategorizedCollection object
#' @details This function sets the \code{CollectionType} for each set
#'     in the collection to {CategorizedCollection}. If
#'     \code{GScollection} contains \code{BroadCollection} gene sets,
#'     their fields \code{category} and \code{subcategory} will be
#'     used. Otherwise the \code{category} and \code{subcategory}
#'     fields will be used.
#' @examples
#' \dontrun{
#' 
#' ## To run this example, first obtain the MSigDB database in XML
#' ## format (see
#' ## http://software.broadinstitute.org/gsea/downloads.jsp). It is
#' ## assumed that the database is locally available as the file
#' ## "msigdb_v6.0.xml".
#'
#' The \code{importMSigDB.xml} function is just a shortcut to the
#' following:
#'
#' db <- getBroadSets("msigdb_v6.1.xml")
#' db <- as.CategorizedCollection(db)
#'
#' ## The database is now in an acceptable format to create a local
#' ## repository using createRepository
#' }
#'
#' ## A small sample of the MSigDB as imported by importMSigDB.xml is
#' ## included in gep2pep. The following creates (and deletes) a
#' ## gep2pep repository.
#'
#' db_sample <- loadSamplePWS()
#' 
#' ## The function \code{as.CategorizedCollection} can also be used to
#' ## create arbitrary gene set collections specifying the categories
#' ## and subcategories once for all the sets:
#'
#' library(GSEABase)
#' mysets <- as.CategorizedCollection(
#'               GeneSetCollection(
#'                   list(GeneSet(c("g1", "g2"), setName="set1"),
#'                        GeneSet(c("g3", "g4"), setName="set2"))
#'                   ),
#'               category="mycategory",
#'               subCategory="mysubcategory"
#'               )
#' newCollection <- GeneSetCollection(c(db_sample, mysets))
#'
#' ## The created repository will include both the sample gene sets
#' ## and the two sets just created:
#' 
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#' rp <- createRepository(repo_path, newCollection)
#' 
#' ## removing temporary repository
#' unlink(repo_path, TRUE)
#' @export
as.CategorizedCollection <- function(GScollection,
                                     category="uncategorized",
                                     subCategory="uncategorized")
{
    if(!is(GScollection, "GeneSetCollection"))
        say("collection must be an object of class GeneSetCollection",
            "error")

    y <- GeneSetCollection(
        sapply(GScollection,
               function(g) {
                   
                   if(is(collectionType(g), "BroadCollection")) {
                       coll <- attributes(collectionType(g))
                       originalID <- setIdentifier(g)
                       collectionType(g) <-
                           CategorizedCollection(coll$category,
                                                 coll$subCategory)
                       setIdentifier(g) <- originalID
                   } else if(is(collectionType(g), "CategorizedCollection")) {
                       ## nothing to do
                   } else {
                       collectionType(g) <- CategorizedCollection(
                           category, subCategory)
                   }
                   g
               })
    )

}

#' Converts gene set IDs to gene set names
#'
#' @param sets An object of class GeneSetCollection
#' @param ids character vector of gene set IDs to be converted to set
#' names
#' @return A vector of gene set names
#' @seealso CondSEA, PathSEA
#' @examples
#' collection <- loadSamplePWS()
#' setId2setName(collection, c("M3128", "M11607"))
#' @export
setId2setName <- function(sets, ids) {
    allids <- sapply(sets, setIdentifier)
    snm <- sapply(sets, setName)
    if(missing(ids))
        return(snm)
    return(snm[match(ids, allids)])
}


#' Imports PEPs created in raw mode
#'
#' Raw mode is meant to deal with large collections of PEPs (like
#' hundreds of thousands). In this case, problems may arise while
#' trying to convert GEPs by loading all of them in memory at
#' once. Raw mode is meant to be used with HDF5 format, which allows
#' to load subsets of GEPs from the disk. \code{buildPEPs}, when used
#' in raw mode, can create the corresponding subsets of PEPs, so that
#' the job can be distributed on a computer
#' cluster. \code{importFromRawMode} is meant to join the chunks into
#' HDF5 matrices, which are than stored into the repository. The
#' \code{.loadPEPs} function can seamlessly load PEPs stored in normal
#' (RDS) or HDF5 format.
#'  
#' @inheritParams dummyFunction
#' @param path Path were raw PEPs are stored (default is a "raw"
#'     directory under the repository root folder).
#' @return Nothing, used for side effects.
#' @seealso buildPEPs
#' @details PEPs are expect to be found at the specified \code{path}
#'     and follow the naming convention as generated by
#'     \code{buildPEPs}. According to such convention, each file is
#'     named usign the format
#'     category_subcategory#chunknumber.RDS. All non-alphanumeric
#'     characters from the original category and subcategory names are
#'     replace with an underscore (in rare cases this could create
#'     ambiguity that should be manually prevented). All chunks for
#'     the same subcategory are joined together following the chunk
#'     numbers into a single HDF5 matrix and stored in the repository
#'     as an "attachment" (see \code{repo} documentation).
#'
#'     Note that raw PEPs (by default everything at
#'     repository_root/raw) can be safly removed once they have been
#'     imported.
#' @export
#' @examples
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#'
#' ## The following will create PEPs in 2 separate files
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps[,1:2], progress_bar=FALSE,
#'     rawmode_id=1)
#' buildPEPs(rp, geps[,3:5], progress_bar=FALSE,
#'     rawmode_id=2)
#'
#' ## The separate files are then merged into one (possibly big) file
#' ## in HDF5 format
#'
#' importFromRawMode(rp)
#'
#' ## Now most operations (excluding the addition of new PEPs to
#' ## existing collections) will be available as usual.
#'
#' unlink(repo_path, TRUE)
importFromRawMode <- function(rp, path=file.path(rp$root(), "raw"),
                              collections="all")
{
    
    allfiles <- list.files(path)
    say(paste0("Found ", length(allfiles), " raw files."))
    dbs <- gsub("#.+[.]RDS", "", allfiles)
    say(paste0("Found ", length(unique(dbs)), " collection names."))
    
    ## real collection names (not stored in filenames)
    collnames <- getCollections(rp)
    cleanedcollnames <- gsub("[^[:alnum:]]", "_", collnames)
    
    for(i in 1:length(unique(dbs))) {
        dbi <- unique(dbs)[i]
        
        collname <- collnames[match(dbi, cleanedcollnames)]

        if(!(length(collections==1) && collections=="all")) {
            if(! collname %in% collections) {
                say(paste0("Collection ", collname,
                           " was not selected and will be skipped."))
                next
            }
        }      

        ## if(rp$has(collname)) {
        ##     say(paste0("A collection named ", collname,
        ##                " is already in the repository ",
        ##                "and will be skipped in raw mode"), "warning")
        ##     next
        ## }
        
        say(paste0("Working on collection: ", collname))        
        
        subfiles <- allfiles[grep(dbi, allfiles)]
        ids <- as.numeric(gsub(".+#|[.]RDS", "", subfiles))
        say(paste0("Found ", length(unique(ids)), " unique IDs."))

        ## extracting chunk size
        fname <- paste0(dbi, "#", min(ids), ".RDS")
        Nchunk <- ncol(readRDS(file.path(path, fname))$ES)
        fname <- paste0(dbi, "#", max(ids), ".RDS")
        NlastChunk <- ncol(readRDS(file.path(path, fname))$ES)
        Ncol <- Nchunk*(length(unique(ids))-1) + NlastChunk
        say(paste0("Expected profiles: ", Ncol, " in chunks of size: ", Nchunk))

        say(paste0("Creating a repository entry."))
        fl <- tempfile()
        h5createFile(fl)
        rp$put(fl, collname, asattach=TRUE, tags=c("pep", "#hdf5"))

        ## attachments are hidden by default, but unhiding them
        ## triggers a warning because of a bug in repo
        curwarn <- options()$warn
        options(warn=-1)
        rp$untag(collname, "hide")
        options(warn=curwarn)
        
        fl <- rp$get(collname)

        fname <- paste0(dbi, "#", min(ids), ".RDS")
        x <- readRDS(file.path(path, fname))$ES
        Nrow <- nrow(x)
        say(paste0("Creating 2 HDF5 dataset of size: ", Nrow, "x", Ncol))

        h5createDataset(fl, "ES", c(Nrow, Ncol),
                        maxdims=c(Nrow*10, Ncol*10),
                        chunk=c(Nrow, Nchunk))
        h5createDataset(fl, "PV", c(Nrow, Ncol),
                        maxdims=c(Nrow*10, Ncol*10),                        
                        chunk=c(Nrow, Nchunk))
        h5createDataset(fl, "rownames", Nrow,
                        maxdims=c(Nrow*10),
                        storage.mode="character",
                        size=256, chunk=Nrow)
        h5createDataset(fl, "colnames", Ncol,
                        maxdims=c(Ncol*10),
                        storage.mode="character",
                        size=256, chunk=Nchunk)
        h5write(rownames(x), fl, "rownames")      
        
        say("Adding chunks...")
        uids <- sort(unique(ids))
        for(j in 1:length(uids)) {
            fsize <- .format.object_size(file.size(fl), "auto")
            fname <- paste0(dbi, "#", uids[j], ".RDS")
            ifile <- file.path(path, fname)
            x <- readRDS(ifile)
            ifsize <- .format.object_size(file.size(ifile), "auto")
            startCol <- (j-1)*Nchunk+1

            h5write(x$ES, fl, "ES", start=c(1,startCol),
                    createnewfile=FALSE)
            h5write(x$PV, fl, "PV", start=c(1,startCol),
                    createnewfile=FALSE)
            h5write(colnames(x$ES), fl, "colnames", start=startCol)

            cat("Chunk ", j, " of ",
                length(uids), ", ", ifsize,
                " (current file size: ", fsize, ")\r",
                sep="")
        }
        say("Done.")
    }

    H5close()
}

#' Imports pathways data from an MSigDB XML file.
#'
#' Creates a \code{GeneSetCollection} object using the XML
#' distribution of the MSigDB (see references). The returned object
#' can be passed to \code{createRepository}.
#'
#' @param fname Path to an XML file downloaded from MSigDB.
#' @param organism Select only gene sets for a given organism. Default
#' is "Homo Sapiens".
#' @return A CategorizedCollection object
#' @references
#'     \url{http://software.broadinstitute.org/gsea/downloads.jsp}
#' @details This function now just calls \code{getBroadSets(fname)}
#'     from the \code{GSEABase} package. However, it is left for
#'     backward compatibility and as an entry point to package
#'     functionalities.
#' @examples
#' \dontrun{
#' 
#' ## To run this example, first obtain the MSigDB database in XML
#' ## format (see
#' ## http://software.broadinstitute.org/gsea/downloads.jsp). It is
#' ## assumed that the database is locally available as the file
#' ## "msigdb_v6.0.xml".
#' 
#' db <- importMSigDB.xml("msigdb_v6.0.xml")
#'
#' ## The database is now in an acceptable format to create a local
#' ## repository using createRepository
#' }
#'
#' ## A small excerpt from the MSigDB is included in gep2pep. The
#' ## following creates (and then deletes) a gep2pep repository.
#'
#' db_sample <- loadSamplePWS()
#' 
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#' rp <- createRepository(repo_path, db_sample)
#'
#' ## removing temporary repository
#' unlink(repo_path, TRUE)
#' @export
importMSigDB.xml <- function(fname, organism="Homo Sapiens") {

    say("Loading gene sets...")

    result = tryCatch({
        sets <- getBroadSets(fname)
        say(paste0("Loaded ", length(sets), " sets."))
        orgs <- sapply(sets, function(s) attributes(s)$organism)        
        if(organism != "all") {
            w <- tolower(orgs) == tolower(organism)
            sets <- sets[w]
            say(paste0("Selected ", length(sets), " sets for: ", organism))
        }
        say("Converting gene sets...")
        as.CategorizedCollection(sets)    
    }, error = function(e) {        
        say(paste("GSEABase::getBroadSets failed with the error:"))
        print(e)
        say("Parsing XML with fallback code...")
        xml <- xmlTreeParse(fname, useInternalNodes=TRUE)   
        sets <- xml["/MSIGDB/GENESET"]

        ids <- sapply(sets, function(x) xmlAttrs(x)[["SYSTEMATIC_NAME"]])
        msigDB <- data.frame(
            id = ids,
            name = sapply(sets, function(x)
                xmlAttrs(x)[["STANDARD_NAME"]]),
            category = sapply(sets, function(x)
                xmlAttrs(x)[["CATEGORY_CODE"]]),
            subcategory = sapply(sets, function(x)
                xmlAttrs(x)[["SUB_CATEGORY_CODE"]]),
            organism = sapply(sets, function(x)
                xmlAttrs(x)[["ORGANISM"]]),
            desc = sapply(sets, function(x)
                xmlAttrs(x)[["DESCRIPTION_BRIEF"]]),
            desc_full = sapply(sets, function(x)
                xmlAttrs(x)[["DESCRIPTION_FULL"]]),
            set = sapply(sets, function(x)
                xmlAttrs(x)[["MEMBERS_SYMBOLIZED"]]),
            stringsAsFactors=FALSE
        )

        say("Converting gene sets...")
        gs <- list()
        k <- 1
        for(i in seq_len(nrow(msigDB))) {
            if(tolower(msigDB[[i]]$organism) == organism) {
                gs[[i]] <- GeneSet(strsplit(msigDB$set[i], ",")[[1]],
                                   shortDescription = msigDB$desc[i],
                                   longDescription = msigDB$desc_full[i],
                                   setName = msigDB$name[i],
                                   setIdentifier = msigDB$id[i],
                                   organism = msigDB$organism[i],
                                   collectionType = CategorizedCollection(
                                       category=msigDB$category[i],
                                       subCategory=msigDB$subcategory[i]
                                   ))
                k <- k+1
            }
        }
        GeneSetCollection(gs)
    })

    say("done.")
    
    return(result)
}



#' Dummy function for parameter inheritance
#' @param rp A repository created by \code{\link{createRepository}}.
#' @param rp_peps A repository created with
#'     \code{\link{createRepository}}, and containing PEPs created
#'     with \code{\link{buildPEPs}}.
#' @param collections A subset of the collection names returned by
#'     \code{getCollections}. If set to "all" (default), all the
#'     collections in \code{rp} will be used.
#' @param collection One of the names returned by
#'     \code{getCollections}.
#' @return Nothing
#' @keywords internal
dummyFunction <- function(rp, rp_peps, collections, collection) {}


#' Check an existyng repository for consistency
#'
#' Check both repository data consistency (see \code{repo_check} from
#' the \code{repo} package) and specific gep2pep data consistency.
#' @inheritParams dummyFunction
#' @return Nothing.
#' @examples
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' checkRepository(rp)
#'
#' unlink(repo_path, TRUE)
#' @export
checkRepository <- function(rp) {
    say("Checking repository consistency...")
    rp$check()

    message()
    say("Cheching for gep2pep data consistency...")
    message()
    
    dbs <- getCollections(rp)
    perts <- rp$get("conditions")
    if(is.null(perts[[1]])) ## this happens due to a bug
        perts <- perts[-1]

    say("Checking conditions list...")
    problems <- FALSE
    
    off <- setdiff(names(perts), dbs)
    if(length(off>0)) {
        say(paste("The following collections are in conditions",
                  "list but not in repository collections:"),
            "warning", off)
        problems <- TRUE
    }

    off <- setdiff(names(perts), dbs)
    if(length(off>0)) {
        say(paste("The following collections are in repository",
                  "collections but not in conditions list:"),
            "warning", off)
        problems <- TRUE
    }

    w <- sapply(lapply(rp$entries(), get, x="tags"), `%in%`, x="pep")
    pepitems <- names(w[w])

    off <- setdiff(pepitems, dbs)
    if(length(off>0)) {
        say(paste("The following collections have PEPs but not",
                  "pathways in the repository:"),
            "warning", off)
        problems <- TRUE
    }

    if(!problems)
        say("Ok.")

    if(length(perts)>0) {

        say("Checking PEPs...")
        for(i in seq_along(dbs)) {
            problems <- FALSE
            
            say(paste0("Checking collection: ", dbs[i], "..."))
            if(rp$has(dbs[i])) {
                peps <- rp$get(dbs[i])

                if(!identical(colnames(peps$ES), perts[[dbs[i]]])) {
                    say(paste("Column names in the PEP matrix differ from",
                              "those in the conditions repository item:",
                              "this is a serious inconsistency!"),
                        "warning")
                    problems <- TRUE
                }

                sets <- .loadCollection(rp, dbs[i])

                if(!identical(colnames(peps$ES), colnames(peps$PV))) {
                    say(paste("Column names of the ES matrix are not",
                              "identical to column names of PV matrix:",
                              "this is a serious inconsistency!"),
                        "warning")
                    problems <- TRUE
                }

                if(!identical(rownames(peps$ES), rownames(peps$PV))) {
                    say(paste("Row names of the ES matrix are not",
                              "identical to row names of PV matrix!",
                              "this is a serious inconsistency!"),
                        "warning")
                    problems <- TRUE
                }

                if(!setequal(names(sets), rownames(peps$PV))) {
                    say(paste("There are pathways in the repository that",
                              "are not in the PEP matrix."),
                        "warning")
                    problems <- TRUE
                }

                if(!problems)
                    say("ok.")
            }
        }

        say(paste("Checking wether different PEP collections",
                  "include the same conditions..."))
        out <- outer(perts, perts,
                     Vectorize(
                         function(a,b) length(intersect(a,b))
                     ))
        if(length(unique(as.vector(out)))!=1) {
            say(paste("Not all collections include the same conditions.",
                      "The intersection matrix follows."))
            print(out)
        } else say("All PEP collections have the same conditions")
        
    } else say("No conditions to check.")
    
}



#' Loads the matrix of Enrichment Scores for a collection
#'
#' @inheritParams dummyFunction
#' @return The matrix of Enrichment Scores (ES) of the
#'     Kolmogorov-Smirnov statistic for the pathway collection, if
#'     previously computed with \code{buildPEPs}. The entry \code{i,j}
#'     reports the ES for the pathway \code{i}, condition{j}. If
#'     \code{buildPEPs} was not run, throws an error.
#' @examples
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' loadESmatrix(rp, "c3_TFT")[1:5,1:2]
#'
#' unlink(repo_path, TRUE)
#'
#' @export
loadESmatrix <- function(rp, collection)
{
    if(!is(collection, "character") || length(collection) > 1)
        say("Please provide a single collection name", "error")
    
    res <- getPEPs(rp, collection)$ES
    return(res)
}


#' Loads the matrix of p-values for a collection
#'
#' @inheritParams dummyFunction
#' @return The matrix of p-values (PV) of the
#'     Kolmogorov-Smirnov statistic for the pathway collection, if
#'     previously computed with \code{buildPEPs}. The entry \code{i,j}
#'     reports the PV for the pathway \code{i}, condition{j}. If
#'     \code{buildPEPs} was not run, throws an error.
#' @examples
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' loadPVmatrix(rp, "c3_TFT")
#'
#' unlink(repo_path, TRUE)
#'
#' @export
loadPVmatrix <- function(rp, collection)
{
    if(!is(collection, "character") || length(collection) > 1)
        say("Please provide a single collection name", "error")
    
    res <- getPEPs(rp, collection)$PV
    return(res)
}

#' Loads a collection of pathways from the repository
#'
#' @inheritParams dummyFunction
#' @return a \code{GeneSetCollection} object loaded from the
#'     repository \code{rp}.
#' @examples
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' 
#' loadCollection(rp, "c3_TFT")
#'
#' unlink(repo_path, TRUE)
#'
#' @export
loadCollection <- function(rp, collection)
{
    if(!is(collection, "character") || length(collection) > 1)
        say("Please provide a single collection name", "error")

    res <- rp$get(paste0(collection, "_sets"))

    return(res)
}



#' Returns the names of the pathway collections in a repository.
#'
#' Given a \code{gep2pep} repository, returns the names of the
#' stored collections by looking at appropriate repository item
#' names.
#' 
#' @inheritParams dummyFunction
#' @return Vector of collection names (see details).
#' @details Each collection in a database has a "category" and a
#'     "subcategory" assigned, which are used to build the collection
#'     identifier as "category_subcategory". This function obtains the
#'     identifiers by looking at data stored in the repository
#'     \code{rp} (entries that are tagged with "sets").
#' @examples
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for collection: c3_TFT
#' ## [15:45:06] Storing pathway data for collection: c3_MIR
#' ## [15:45:06] Storing pathway data for collection: c4_CGN
#'
#' getCollections(rp)
#' ## [1] "c3_TFT" "c3_MIR" "c4_CGN"
#'
#' unlink(repo_path, TRUE)
#'
#' @export
getCollections <- function(rp)
{    
    w <- sapply(lapply(rp$entries(), get, x="tags"), `%in%`, x="sets")
    res <- gsub("_sets", "", names(w[w]))
    return(res)
}


#' Creates a repository of pathway collections.
#'
#' Given a database of collections, stores them in a local repository
#' to be used by \code{gep2pep} functions.
#'
#' @param path Path to a non-existing directory where the repository
#'     will be created.
#' @param sets An object of class \code{CategorizedCollection}.
#' @param name Name of the repository. Defaults to \code{NULL} (a
#'     generic name will be given).
#' @param description Description of the repository. If NULL
#'     (default), a generic description will be given.
#' @return An object of class \code{repo} that can be passed to
#'     \code{gep2pep} functions.
#' @details \code{sets} can be created by
#'     \code{\link{importMSigDB.xml}} or using \code{GSEABase}
#'     \code{GeneSetCollection} class and then converting it to
#'     CategorizedCollection. See examples.
#' @seealso buildPEPs
#' @examples
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for collection: c3_TFT
#' ## [15:45:06] Storing pathway data for collection: c3_MIR
#' ## [15:45:06] Storing pathway data for collection: c4_CGN
#'
#' rp
#' ##         ID    Dims     Size
#' ## c3_TFT_sets   10 18.16 kB
#' ## c3_MIR_sets   10 17.25 kB
#' ## c4_CGN_sets   10   6.9 kB
#'
#' unlink(repo_path, TRUE)
#' @export
createRepository <- function(path, sets, name=NULL, description=NULL)
{
    if(!is.null(sets) && !is(sets, "GeneSetCollection"))
        say("sets must be an object of class GeneSetCollection", "error")

    types <- sapply(sets, collectionType)
    allCat <- all(sapply(types, is, class2="CategorizedCollection"))
    if(!allCat) {
        say("all sets must be of type CategorizedCollection",
            "error")
    }
    
    if(file.exists(path)) {
        say("Can not create repository in existing folder", "error")
    } else rp <- repo_open(path, TRUE)

    if(is.null(name))
        name <- "gep2pep repository"
    if(is.null(description))
        description <- paste("This repository contains pathway information",
                             "and possibly Pathway Expression Profiles", 
                             "created with the gep2pep package.")    
    
    
    rp$project(name, description)

    ## un-hiding the project item
    ## silencing warning because of a bug in repo
    curwarn <- options()$warn
    options(warn=-1)
    rp$untag(name, "hide")
    options(warn=curwarn)

    if(!is.null(sets)) {
        db_ids <- makeCollectionIDs(sets)
        subdbs <- unique(db_ids)
        
        for(i in seq_along(subdbs))
        {
            dbi <- subdbs[i]
            say(paste("Storing pathway data for collection:", dbi))
            rp$put(sets[which(db_ids == dbi)],
                   paste0(subdbs[i], "_sets"),
                   paste("Pathway information for collection", subdbs[i]),
                   c("gep2pep", "sets"),
                   prj = name)
        }
    }

    rp$put(list(NULL), "conditions",
           "Condition lists for PEPs",
           c("gep2pep", "perts"),
           prj = name)
    
    return(invisible(rp))
}


#' Opens an existing repository of pathway collections.
#'
#' The repository must have been created by
#' \code{\link{createRepository}}. Provides an R object to interact
#' with the repository.
#'
#' @param path Path to a directory where the repository has been
#'     created with \code{\link{createRepository}}.
#' @return An object of class \code{repo} that can be passed to
#'     \code{gep2pep} functions.
#' @details This function only calls the \code{repo_open} function
#'     from the \code{repo} package on \code{path}. It is meant to
#'     allow users not to explicitly load the \code{repo} library,
#'     unless they want to access advanced features.
#' @seealso createRepository
#' @examples
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' rp2 <- openRepository(repo_path)
#'
#' ## rp and rp2 point to the same data:
#' identical(rp$entries(), rp2$entries())
#' ## > [1] TRUE
#'
#' unlink(repo_path, TRUE)
#' @export
openRepository <- function(path)
{
    if(!file.exists(path)) {
        say(paste("path must point to an existing directory containing",
                  "a repository created with createRepository"), "error")
    }
    rp <- repo_open(path, TRUE)
}



#' Creates a collection label for each pathway.
#'
#' Given a database, uses "category" and "subcategory" entries to
#' create a vector of collection identifiers. Useful to extract a
#' collection from a database.
#'
#' @param sets A pathway database in the same format as output by
#'     \code{importMSigDB.xml}.
#' @return A vector of identifiers, one per pathway, with the format:
#'     "category_subcategory".
#' @seealso importMSigDB.xml
#' @examples
#' db <- loadSamplePWS()
#' ids <- makeCollectionIDs(db)
#'
#' unique(ids)
#' ## [1] "c3_TFT" "c3_MIR" "c4_CGN"
#'
#' db <- db[ids=="c3_MIR"]
#'
#' length(db)
#' ## [1] 10
#'
#' @export
makeCollectionIDs <- function(sets) {
    if(!is(sets, "GeneSetCollection"))
        say("sets must be an object of class GeneSetCollection")
    sets <- convertFromGSetClass(sets)
    
    dbs <- sapply(sets, get, x="category")
    subdbs <- sapply(sets, get, x="subcategory")
    w <- which(subdbs=="" | is.na(subdbs)) 
    if(length(w)>0)
        subdbs[w] <- dbs[w]
    db_ids <- paste(dbs, subdbs, sep="_")
    return(db_ids)
}


#' Build PEPs from GEPs and stores them in the repository.
#'
#' Given a matrix of ranked lists of genes (GEPs) and a \code{gep2pep}
#' repository, converts GEPs to PEPs and stores the latter in the
#' repository.
#'
#' @inheritParams dummyFunction
#' @param geps A matrix of ranks where each row corresponds to a gene
#'     and each column to a condition. Each column must include all
#'     ranks from 1 to the number of rows. Row and column names must
#'     be defined. Row names will be matched against gene identifiers
#'     in the pathways collections, and unrecognized gene names will
#'     not be used.
#' @param min_size An integer representing the minimum number of genes
#'     that must be included in a set before the KS statistic is
#'     computed. Smaller gene sets will get ES=NA and p=NA. Default is
#'     3. Ignored for SGE mode (see \code{addSingleGeneSets}).
#' @param max_size An integer representing the maximum number of genes
#'     that must be included in a set before the KS statistic is
#'     computed. Larger gene sets will get ES=NA and p=NA. Default is
#'     500.
#' @param parallel If TRUE, gene sets will be processed in
#'     parallel. Requires a parallel backend.
#' @param replace_existing What to do if PEPs, identified by column
#'     names of \code{geps} are already present in the repository. If
#'     set to TRUE, they will be replaced, otherwise they will be
#'     skipped and only PEPs of new conditions will be computed and
#'     added. Either ways, will throw a warning.
#' @param progress_bar If set to TRUE (default) will show a progress
#'     bar updated after coversion of each column of \code{geps}.
#' @param rawmode_id An integer to be appended to files produced in
#'     raw mode (see details). If set to NULL (default), raw mode is
#'     turned off.
#' @param rawmode_outdir A charater vector specifying the destination
#'     path for files produced in raw mode (by the fault it is
#'     ROOT/raw, where ROOT is the root of the repository). Ignored if
#'     \code{rawmode_id} is NULL.
#' @return Nothing. The computed PEPs will be available in the
#'     repository.
#' @seealso buildPEPs
#' @details By deault, output is written to the repository as new
#'     items named using the collection name. However, it is possible
#'     to avoid the repository and write the output to regular files
#'     turning 'raw mode' on through the \code{rawmode_id} and
#'     \code{rawmode_outdir} parameters. This is particuarly useful
#'     when dealing with very large corpora of GEPs, and conversions
#'     are split into independent jobs submitted to a scheduler. At
#'     the end, the data will need to be reconstructed and put into
#'     the repository using \code{importFromRawMode} in order to
#'     perform \code{CondSEA} or \code{PathSEA} analysis.
#' @examples
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for collection: c3_TFT
#' ## [15:45:06] Storing pathway data for collection: c3_MIR
#' ## [15:45:06] Storing pathway data for collection: c4_CGN
#'
#' rp
#' ##          ID   Dims     Size
#' ## c3_TFT_sets   10   18.16 kB
#' ## c3_MIR_sets   10   17.25 kB
#' ## c4_CGN_sets   10     6.9 kB
#'
#' ## Loading sample gene expression profiles
#' geps <- loadSampleGEP()
#'
#' geps[1:3,1:3]
#' ##       (+)_chelidonine (+)_isoprenaline (+/_)_catechin
#' ## AKT3               88              117            417
#' ## MED6              357              410             34
#' ## NR2E3             383              121            453
#'
#' buildPEPs(rp, geps)
#'
#' rp
#' ##           ID  Dims     Size
#' ##   c3_TFT_sets   10 18.16 kB
#' ##   c3_MIR_sets   10 17.25 kB
#' ##   c4_CGN_sets   10   6.9 kB
#' ##        c3_TFT    2  1.07 kB
#' ##        c3_MIR    2  1.07 kB
#' ##        c4_CGN    2  1.04 kB
#'
#' unlink(repo_path, TRUE)
#'
#' @export
buildPEPs <- function(rp, geps, min_size=3, max_size=500,
                      parallel=FALSE, collections="all",
                      replace_existing=FALSE, donotstore=FALSE,
                      progress_bar=TRUE,
                      rawmode_id=NULL,
                      rawmode_outdir=file.path(rp$root(), "raw"))
{
    checkGEPsFormat(geps)
    perts <- rp$get("conditions")
    rawmode <- !is.null(rawmode_id)

    if(replace_existing && donotstore)
        say("either replace_existing or donotstore can be true",
            type="error")

    ## the following is to force recomputing existing PEPs
    if(donotstore)
        replace <- existing <- TRUE 
    
    if(rawmode && rawmode_id %% 1 != 0)
        say("rawmode_id must be an integer number", type="error")
    
    if(length(collections) == 1 && collections == "all") {
        dbs <- getCollections(rp)
    } else dbs <- collections

    if(donotstore)
        res <- list()
    
    for(i in seq_along(dbs))
    {
        say(paste0("Working on collection: ", dbs[i],
                   " (", i, "/", length(dbs), ")" ))

        if(rp$has(dbs[i]))
            curpeps <- perts[[dbs[i]]] else curpeps <- NULL
        newpeps <- setdiff(colnames(geps), curpeps)
        oldpeps <- intersect(colnames(geps), curpeps)

        if(length(oldpeps > 0)) {
            if(rawmode) {
                say(paste0(length(oldpeps),
                           " existing PEPs found, ",
                           "but this will be ignored in rawmode: ",
                           paste(oldpeps, collapse=", ")))
                newpeps <- colnames(geps)
            } else {
                msg <- paste0(length(oldpeps),
                              " existing PEPs will be skipped: ",
                              paste(oldpeps, collapse=", "))
                if(replace_existing) {
                    msg <- gsub("skipped", "replaced", msg)
                    newpeps <- colnames(geps)
                }
                say(msg, type="warning")
            }
        }

        if(length(newpeps) > 0) {
            gepsi <- geps[, newpeps, drop=FALSE]
            thisdb <- .loadCollection(rp, dbs[i])

            SGEmode <- "SGE" %in% rp$tags(paste0(dbs[i], "_sets"))
            if(SGEmode) {
                say("Running in Single Gene mode")
                peps <- gep2pep(gepsi, thisdb, min_size=1, max_size,
                                parallel, progress_bar, SGEmode=TRUE)
            } else {
                peps <- gep2pep(gepsi, thisdb, min_size, max_size,
                                parallel, progress_bar, SGEmode=FALSE)
            }

            if(donotstore) {
                res[[dbs[i]]] <- peps
            } else {
                storePEPs(rp, dbs[i], peps, rawmode_id,
                          rawmode_outdir)
            }
        }
    }
    if(donotstore)
        return(res)
}



#' Export CondSEA or PathSEA results to XLS format
#'
#' The XLS output includes the full CondSEA or PathSEA results,
#' together with additional gene set information for the CondSEA. If
#' the PathSEA or CondSEA analysis was performed with
#' \code{details=TRUE}, details will be reported in the XLS file. This
#' function requires the WriteXLS library.
#'
#' @inheritParams dummyFunction
#' @param results The output of \code{CondSEA} or \code{PathSEA}.
#' @param outname Name of the XLS file to be created.
#' @return Nothing.
#' @seealso CondSEA, PathSEA
#' @examples
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' pgset <- c("(+)_chelidonine", "(+/_)_catechin")
#' psea <- CondSEA(rp, pgset)
#'
#' \dontrun{
#' exportSEA(rp, psea)
#' }
#'
#' unlink(repo_path, TRUE)
#'
#' @export
exportSEA <- function(rp, results, outname=NULL)
{
    
    type <- names(results)[[1]]
    if(! tolower(type) %in% c("condsea", "pathsea"))
        say("type must be on of: CondSEA, PathSEA", "error")

    if(is.null(outname))
        outname <- paste0(type, ".xls")
    
    if (requireNamespace("WriteXLS", quietly = TRUE)) {
        sheets <- attachInfo(rp, results)
        names(sheets) <- gsub(":","_",names(sheets))
        WriteXLS::WriteXLS(sheets, outname, AutoFilter=TRUE,
                           BoldHeaderRow=TRUE, FreezeRow=1)
    } else {
        stop("The suggested package WriteXLS is not installed.")
    }              
}


#' Extracts the results matrix from \code{CondSEA} or \code{PathSEA}
#' output
#'
#' @inheritParams dummyFunction
#' @param analysis The output of either \code{CondSEA} or
#'     \code{PathSEA}.
#' @return A 2-columns matrix including ESs and
#'     p-values (see details) for each pathway database and
#'     condition.
#' @seealso CondSEA, PathSEA
#' @examples
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' pgset <- c("(+)_chelidonine", "(+/_)_catechin")
#' psea <- CondSEA(rp, pgset)
#'
#' getResults(psea, "c3_TFT")
#'
#' unlink(repo_path, TRUE)
#'
#' @export
getResults <- function(analysis, collection)
{
    if(!is(collection, "character") || length(collection) > 1)
        say("Please provide a single collection name", "error")

    what <- c("CondSEA", "PathSEA") %in% names(analysis)
    if(
        !is(analysis, "list") ||
        !any(what) ||
        !("details" %in% names(analysis))
    ) say("analysis does not look like the output of CondSEA or PathSEA",
          "error")

    wanalysis <- c("CondSEA", "PathSEA")[what]
    if(! collection %in% names(analysis[[wanalysis]]))
        say(paste("There are no results for the collection:", collection),
            "error")

    return(analysis[[wanalysis]][[collection]])
}


#' Extracts the details matrix from \code{CondSEA} or \code{PathSEA}
#' output
#' @inheritParams dummyFunction
#' @param analysis The output of either \code{CondSEA} or
#'     \code{PathSEA}.
#' @return A matrix including the ranks of each pathway (over rows)
#'     and each condition (over columns) used as input to
#'     \code{CondSEA} or \code{PathSEA}.
#' @seealso CondSEA, PathSEA
#' @examples
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' pgset <- c("(+)_chelidonine", "(+/_)_catechin")
#' psea <- CondSEA(rp, pgset)
#'
#' getDetails(psea, "c3_TFT")
#'
#' unlink(repo_path, TRUE)
#'
#' @export
getDetails <- function(analysis, collection)
{
    if(!is(collection, "character") || length(collection) > 1)
        say("Please provide a single collection name", "error")

    what <- c("CondSEA", "PathSEA") %in% names(analysis)
    if(
        !is(analysis, "list") ||
        !any(what) ||
        !("details" %in% names(analysis))
    ) say("analysis does not look like the output of CondSEA or PathSEA",
          "error")

    wanalysis <- c("CondSEA", "PathSEA")[what]
    if(! collection %in% names(analysis[[wanalysis]]))
        say(paste("There are no results for the collection:", collection),
            "error")    
    
    return(analysis[["details"]][[collection]])
}


.loadPerts <- function(rp, coll) {
    ish5 <- "#hdf5" %in% rp$tags(coll)

    if(ish5){
        perts <- h5read(rp$get(coll), "colnames")
    } else perts <- colnames(rp$get(coll)$ES)
    return(perts)
}

.loadPEPs <- function(rp, coll, subset) {
    ish5 <- "#hdf5" %in% rp$tags(coll)

    if(ish5) {
        fname <- rp$get(coll)
        hcolnames <- h5read(fname, "colnames")
        
        if(missing(subset))
            subset <- hcolnames

        if(is.character(subset))
            subset <- match(subset, hcolnames)
        
        peps <- list(
            ES = h5read(fname, "ES", index=list(NULL, subset)),
            PV = h5read(fname, "PV", index=list(NULL, subset))
        )
        rownames(peps$ES) <- rownames(peps$PV) <-
            h5read(fname, "rownames")
        colnames(peps$ES) <- colnames(peps$PV) <-
            hcolnames[subset]
    } else {
        peps <- list(
            ES = rp$get(coll)$ES[, subset, drop=FALSE],
            PV = rp$get(coll)$PV[, subset, drop=FALSE]
        )
    }
    
    return(peps)
}


#' Performs Condition Set Enrichment Analysis
#'
#' Condition Set Enrichment Analysis (CondSEA) can be seen as a
#' Gene-SEA performed over rows (as opposed to columns) of a matrix of
#' GEPs. It tells how much a pathway is consistently dysregulated
#' under a set of conditions (such as a set of drug treatments,
#' disease states, cell types, etc.) when compared to a statistical
#' background of other conditions.
#'
#' @inheritParams dummyFunction
#' @param pgset A vector of names of conditions. Corresponding PEPs
#'     must exist in all the pathway collections currently in
#'     \code{rp}.
#' @param bgset The background against which to compare
#'     \code{pgset}. If set to \code{all} (default), all the remaining
#'     PEPs will be used. If provided, the corresponding PEPs must
#'     exist in all the pathway collections currently in \code{rp}.
#' @param details If TRUE (default) rank details will be reported for
#'     each condition in \code{pgset}.
#' @param rankingFun The function used to rank PEPs column-wise. By
#' default \code{rankPEPsByRows.ES} is used, which ranks using gene
#' set enrichment scores (see details).
#' @param usecache If set to TRUE, the computed ranked matrix will be
#' stored in the the repository (see details). FALSE by default.
#' @param sortoutput If TRUE (default) the output gene sets will be
#' sorted in order of increasing p-value.
#' @return A list of 2, by names "CondSEA" and "details". The
#' "CondSEA" entry is a 2-columns matrix including ESs and p-values
#' (see details) for each pathway database and condition. The
#' "details" entry reports the rank of each condition in \code{pgset}
#' for each pathway.
#' 
#' @details For each pathway, all conditions are ranked by how much
#' they dysregulate it (from the most UP-regulating to the most
#' DOWN-regulating). Then, a Kolmogorov-Smirnov (KS) test is performed
#' to compare the ranks assigned to conditions in \code{pgset} against
#' the ranks assigned to conditions in \code{bgset}. A positive
#' (negative) Enrichment Score (ES) of the KS test indicates whether
#' each pathway is UP- (DOWN-) regulated by \code{pgset} as compared
#' to \code{bgset}. A p-value is associated to the ES.
#'
#' When PEPs are obtained from drug-induced gene expression profiles,
#' \code{PathSEA} is the Drug-Set Enrichment Analysis [1].
#'
#' The \code{rankingFun} must take in input PEPs like those loaded
#' from the repository and return a matrix of row-wise ranks. Each row
#' must contains ranks from 1 to the number of PEPs minus the number
#' of NAs in the row.
#' 
#' When \code{usecache=TRUE}, the ranked matrix is permanently stored
#' in HDF5 format, and subsequent calls to \code{CondSEA} will load
#' from the disk the necessary ranks (not the whole matrix). The
#' correct cached data is identified by the alphabetically sorted set
#' \code{union(pgset, bgset)}, by the collection name, and by the
#' ranking function. Additional alls to CondSEA with variations of
#' these inputs will create additional cache. Cached data is hidden in
#' the repository by default and can be printed with
#' \code{rp_peps$print(all=TRUE)}, and cleared with
#' \code{clearCache(rp_peps)}.
#'
#'     
#' @seealso getResults, getDetails, clearCache
#' @references
#'
#' [1] Napolitano F. et al, Drug-set enrichment analysis: a
#'     novel tool to investigate drug mode of action. Bioinformatics
#'     32, 235-241 (2016).
#'
#' @examples
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' pgset <- c("(+)_chelidonine", "(+/_)_catechin")
#' psea <- CondSEA(rp, pgset)
#'
#' res <- getResults(psea, "c3_TFT")
#'
#' ## getting the names of the top pathways
#'
#' setId2setName(loadCollection(rp, "c3_TFT"), rownames(res))
#' 
#' unlink(repo_path, TRUE)
#'
#' @export
CondSEA <- function(rp_peps, pgset, bgset="all", collections="all",
                    details=TRUE, rankingFun = rankPEPsByRows.ES,
                    usecache=FALSE, sortoutput=TRUE)
{
    dbs <- collections

    if(!is.character(pgset) || !is.character(bgset))
        say("pgset and bgset must be of class character",
           "error")
   
    if(length(dbs) == 1 && dbs=="all") {
        dbs <- getCollections(rp_peps)
    } else {
        off <- setdiff(dbs, getCollections(rp_peps))
        if(length(off)>0)
            say(paste0("The following collections could not be found: ",
                       paste(off, collapse=", ")), "error")
    }

    if(details)
        thedetails <- list() else thedetails <- NULL
    
    res <- list()
    
    for(i in seq_along(dbs)) {
        say(paste0("Working on collection: ", dbs[i]))

        allperts <- .loadPerts(rp_peps, dbs[i])

        checkMissingPerts <- function(testset, backset) {
            if(!all(testset %in% backset))
                say(paste("The following conditions could not be found:",
                          paste(
                              setdiff(testset, backset),
                              collapse = ", ")),
                    "error")
        }

        checkMissingPerts(pgset, allperts)
        if(length(bgset) == 1 && bgset=="all") {
            bgset <- allperts
        }
        checkMissingPerts(bgset, allperts)                

        if(length(intersect(pgset, bgset))>0) {
            bgset <- setdiff(bgset, pgset)
            say("Common conditions removed from bgset")
        }
        
        rankingset <- union(bgset, pgset)
        
        ranked <- .getPEPRanks(
            rp_peps, dbs[i], rankingset,
            subcols=pgset,
            rankingFun=rankingFun,
            usecache=usecache
        )
        
        say(paste0("Computing enrichments"))
        
        ks <- .doPSEA(ranked, length(rankingset))

        PVs <- sapply(ks, get, x="p.value")
        
        if(sortoutput) {
            sorter <- order(PVs)
        } else sorter <- seq_along(PVs)
        res[[dbs[i]]] <- data.frame(ES = sapply(ks, get, x="ES"),
                                    PV = PVs)[sorter, ]
        if(details)
            thedetails[[dbs[i]]] <- ranked[sorter, ]
        say("Done.")
    }

    return(list(CondSEA=res, details=thedetails))
}

.doPSEA <- function(ranked, backgroundSize, parallel=FALSE) {
    nas <- attributes(ranked)$nas

    '%dobest%' <- if (parallel) get('%dopar%') else get('%do%')
    set <- NULL ## to cope with R CMD check NOTE
 
    ##    ks <- sapply(1:nrow(ranked), function(i) {
    ks <- foreach(row = iter(ranked, by="row"),
                  nas = iter(nas)) %dobest% {
                      innas <- is.na(row)
                      inset <- row[!innas]
                      maxr <- backgroundSize - nas
                      outset <- (1:maxr)[-inset]
                      if(length(inset)>=1 && length(outset)>=1) {
                          res <- ks.test.2(inset, outset, maxCombSize=10^10)
                      } else res <- list(ES=NA, p.value=NA)
                      res
                  }
                      
    names(ks) <- rownames(ranked)
    return(ks)
}

#' Clear cached ranked matrices
#'
#' @inheritParams dummyFunction
#' @details This will clear everything in the repository tagged with
#' "stashed", which by default includes only matrices ranked by some
#' gep2pep functions such as \code{CondSEA}.
#' @seealso CondSEA
#' @return Nothing, used for side effects
#' @examples
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' pgset <- c("(+)_chelidonine", "(+/_)_catechin")
#' psea <- CondSEA(rp, pgset, usecache=TRUE)
#'
#' ## the repository contains cached data
#' print(rp, all=TRUE)
#'
#' clearCache(rp)
#'
#' unlink(repo_path, TRUE)
#' @export
clearCache <- function(rp_peps)
{
    if("stash" %in% rp_peps$tags())
        suppressWarnings(rp_peps$stashclear(TRUE))
}

.storeCachedRanks <- function(rp_peps, name, ranks) {
    fl <- tempfile()
    h5createFile(fl)
    rp_peps$put(fl, name, asattach=TRUE, tags=c("stash", "#hdf5"))        
    fl <- rp_peps$get(name)
    Nrow <- nrow(ranks)
    Ncol <- ncol(ranks)
    ranks[is.na(ranks)] <- -1
    h5createDataset(fl, "ranks", c(Nrow, Ncol),
                    maxdims=c(Nrow*10, Ncol*10),
                    chunk=c(Nrow, Ncol))
    l <- length(attr(ranks, "nas"))
    h5createDataset(fl, "nas", l,
                    maxdims=l*10,
                    chunk=length(attr(ranks, "nas")))
    h5write(ranks, fl, "ranks")
    h5write(attr(ranks, "nas"), fl, "nas")
    h5write(colnames(ranks), fl, "colnames")
    h5write(rownames(ranks), fl, "rownames")
    H5close()
}

.getPEPRanks <- function(rp_peps, coll, rankingset,
                         subrows, subcols,
                         rankingFun, usecache=FALSE)
{
    if(!is.character(rankingset) ||
       !is.character(subcols))
        say("rankingset and subcols must be of class character",
            "error")
    meta <- list(coll, sort(rankingset), deparse(substitute(rankingFun)))
    md5 <- digest(meta)
    if(rp_peps$has(md5) && usecache) {        
        say("Loading cached ranks")
        if(missing(subrows)) subrows <- NULL
        if(missing(subcols)) subcols <- NULL
        if(!is.null(subrows)) {
            rnames <- h5read(rp_peps$get(md5), "rownames")
            subrows <- match(subrows, rnames)
        }
        if(!is.null(subcols)) {
            cnames <- h5read(rp_peps$get(md5), "colnames")
            subcols <- match(subcols, cnames)
        }
        ranked <- h5read(rp_peps$get(md5), "ranks",
                         index=list(subrows, subcols))
        ranked[ranked==-1] <- NA
        nas <- h5read(rp_peps$get(md5), "nas")
        attr(ranked, "nas") <- nas
    } else {
        peps <- .loadPEPs(rp_peps, coll, rankingset)
        say(paste0("Ranking collection"))
        ranked <- rankingFun(peps)

        if(usecache) {
            say(paste0("Caching ranks"))
            .storeCachedRanks(rp_peps, md5, ranked)
        }
        ## clean this
        nas <- attr(ranked, "nas")
        ranked <- ranked[subrows, subcols, drop=FALSE]
        attr(ranked, "nas") <- nas
    }
    return(ranked)
}

#' Performs Pathway Set Enrichment Analysis (PSEA)
#'
#' PathSEA is analogous to the Gene Set Enrichment Analysis (GSEA),
#' but for pathways instead of single genes. It can therefore be used
#' to look for conditions under which a given set of pathways is
#' consistently UP- or DOWN-regulated.
#'
#' @inheritParams dummyFunction
#' @param pathways A database of pathways in the same format as input
#'     to \code{\link{createRepository}}. PSEA will be performed for
#'     each database separately.
#' @param bgsets Another list like \code{pathways}, representing the
#'     statistical background for each database. If set to "all" (the
#'     default), all pathways that are in the repository and not in
#'     \code{pathways} will be used.
#' @param subset Character vector including PEP names to be
#'     considered (all by default, which may take time).
#' @param details If TRUE (default) details will be reported for each
#'     condition in \code{pgset}.
#' @param rankingFun The function used to rank PEPs column-wise. By
#'     default \code{rankPEPsByCols.ES} is used, which uses gene set
#'     enrichment scores (see details).
#' @return A list of 2, by names "PathSEA" and "details". The
#'     "PathSEA" entry is a 2-columns matrix including ESs and
#'     p-values for each collection and condition. The "details" entry
#'     reports the rank of each pathway in \code{pathways} for each
#'     condition.
#' @details For each condition, all pathways are ranked by how much
#'     they are dysregulated by it (from the most UP-regulated to the
#'     most DOWN-regulatied, according to the corresponding
#'     p-values). Then, a Kolmogorov-Smirnov (KS) test is performed to
#'     compare the ranks assigned to pathways in \code{pathways}
#'     against the ranks assigned to pathways in \code{bgsets}. A
#'     positive (negative) Enrichment Score (ES) of the KS test
#'     indicates whether each pathway is UP- (DOWN-) regulated by
#'     \code{pgset} as compared to \code{bgset}. A p-value is
#'     associated to the ES.
#'
#'     When PEPs are obtained from drug-induced gene expression
#'     profiles, \code{PathSEA} can be used together with
#'     \code{gene2pathways} to perform gene2drug [1] analysis, which
#'     predicts which drugs may target a gene of interest (or mimick
#'     such effect).
#' 
#'     The \code{rankingFun} must take in input PEPs like those loaded
#'     from the repository and return a matrix of column-wise
#'     ranks. Each column must contain ranks from 1 to the number of
#'     gene sets minus the number of NAs in the column.
#'
#' @seealso getResults, getDetails
#' @references
#'     [1] Napolitano, F. et al. gene2drug: a computational tool for
#'         pathway-based rational drug repositioning. Bioinformatics
#'         (2017). https://doi.org/10.1093/bioinformatics/btx800
#' @examples
#' library(GSEABase)
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' pathways <- c("M11607", "M10817", "M16694",         ## from c3_TFT
#'               "M19723", "M5038", "M13419", "M1094") ## from c4_CGN
#' w <- sapply(db, setIdentifier) %in% pathways
#'
#' psea <- PathSEA(rp, db[w])
#' ## [15:35:29] Working on collection: c3_TFT
#' ## [15:35:29] Common pathway sets removed from bgset.
#' ## [15:35:29] Column-ranking collection...
#' ## [15:35:29] Computing enrichments...
#' ## [15:35:29] done.
#' ## [15:35:29] Working on collection: C4_CGN
#' ## [15:35:29] Common pathway sets removed from bgset.
#' ## [15:35:29] Column-ranking collection...
#' ## [15:35:29] Computing enrichments...
#' ## [15:35:29] done.
#'
#' getResults(psea, "c3_TFT")
#' ##                         ES        PV
#' ## (_)_mk_801       0.7142857 0.1666667
#' ## (_)_atenolol     0.7142857 0.1666667
#' ## (+)_isoprenaline 0.5714286 0.4000000
#' ## (+/_)_catechin   0.5714286 0.4000000
#' ## (+)_chelidonine  0.3333333 0.9333333
#'
#' unlink(repo_path, TRUE)
#'
#' @export
PathSEA <- function(rp_peps, pathways, bgsets="all",
                    collections="all", subset="all", details=TRUE,
                    rankingFun=rankPEPsByCols.SPV)
{
    checkSets(rp_peps, pathways)
    
    pathways <- convertFromGSetClass(pathways)
    
    pathways <- pwList2pwStruct(pathways)

    if(length(bgsets)==1 && bgsets != "all") {
        checkSets(bgsets)
        bgsets <- pwList2pwStruct(bgsets)
    }

    
    if(!(length(collections) == 1 && collections=="all")) {

        if(!all(collections %in% names(pathways)))
            say(paste("There is at least one selected collections for which",
                      "no pathway has been provided"), "warning")
        
        offcols <- setdiff(collections, getCollections(rp_peps))
        if(length(offcols) > 0)
            say("The following collections are not in the repository:",
                "error", offcols)
    } else collections <- getCollections(rp_peps)

    for(i in seq_along(pathways)) {
        dbi <- names(pathways)[i]
        if(dbi %in% collections && !rp_peps$has(dbi))
            say("Could not find PEPs: ", "error", names(pathways)[i])
    }
    
    collections <- intersect(names(pathways), collections)
    
    if(length(setdiff(names(pathways), collections)>1)) {
        say(paste("Removing pathways not in specified collections"))
        pathways <- pathways[collections]
    }
    
    if(details)
        thedetails <- list() else thedetails <- NULL     
    
    res <- list()
    for(i in seq_along(pathways))
    {
        say(paste0("Working on collection: ", collections[i]))
        gmd <- names(pathways[[collections[i]]])

        allsets <- names(.loadCollection(rp_peps, collections[i]))
        
        if(length(bgsets) == 1 && bgsets=="all") {
            bgset <- allsets
        } else bgset <- names(pwList2pwStruct(bgsets)[[i]])

        if(length(intersect(gmd, bgset)) > 0) {
            bgset <- setdiff(bgset, gmd)
            say("Common pathway sets removed from bgset")
        }
        rankingset <- c(gmd, bgset)
        
        if(length(subset)==1 && subset == "all") {
            peps <- .loadPEPs(rp_peps, collections[i])
        } else peps <- .loadPEPs(rp_peps, collections[i], subset)

        notok <- rankingset[rankingset %in% rownames(peps)]
        if(length(notok)>0)
            say(paste0("Pathway set ids not found in ", collections[i], ": ",
                       paste(notok, collapse=", ")), "error")

        say(paste0("Column-ranking collection"))
        ranked <- rankingFun(peps, rankingset)

        say(paste0("Computing enrichments"))
        ks <- apply(ranked, 2, function(col) {
            inset <- col[gmd]
            inset <- inset[!is.na(inset)]
            outset <- col[bgset]
            outset <- outset[!is.na(outset)]
            if(length(inset)>1 && length(outset)>1) {
                res <- ks.test.2(inset, outset, maxCombSize=10^10)
            } else res <- list(ES=NA, p.value=NA)
            res
        })
        
        PVs <- sapply(ks, get, x="p.value")
        sorter <- order(PVs)
        
        res[[collections[i]]] <- data.frame(
            ES = sapply(ks, get, x="ES")[sorter],
            PV = PVs[sorter]
        )

        if(details)
            thedetails[[collections[i]]] <- ranked[gmd, sorter]
        
        say("done")
    }

    return(list(PathSEA=res, details=thedetails))
}


#' Finds pathways including a given gene.
#'
#' Given a gene, find the set of pathways that involve it in each
#' collection of the repository. This can be used to define a set of
#' pathways for the \code{\link{PathSEA}}.
#'
#' @inheritParams dummyFunction
#' @param genes A vector of gene identifiers of the same type as that
#'     used to create the repository.
#' @param and If set to TRUE (default), will return sets containing
#'     all of \code{genes}. Otherwise will return the sets containing
#'     any of \code{genes}.
#' @return A database of pathways suitable as input to
#'     \code{\link{PathSEA}}.
#' @seealso createRepository, PathSEA
#' @examples
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#'
#' ## Finding all pathways containing "FAM126A":
#' subpw <- gene2pathways(rp, "FAM126A")
#'
#' print(names(subpw))
#'
#' unlink(repo_path, TRUE)
#'
#' @export
gene2pathways <- function(rp, genes, and=TRUE)
{
    dbs <- getCollections(rp)
    mods <- list()
    for(i in seq_along(dbs)) {
        db <- loadCollection(rp, dbs[i])
        if(and) {
            w <- sapply(db, function(s) all(genes %in% geneIds(s)))
        } else w <- sapply(db, function(s) any(genes %in% geneIds(s)))
        mods <- GeneSetCollection(c(mods, db[w]))
    }

    return(mods)        
}


## convert2SummarizedExperiment <- function(rp, res, bg)
## {
##     restype <- names(res)[[1]]
##     SEs <- list()
##     for(i in 1:length(res)) {
##         dfi <- t(data.frame(res$CondSEA[[i]], res$details[[i]]))
##         conds <- data.frame(originalNames=c(
##                                 colnames(res$CondSEA[[i]]),
##                                 colnames(res$details[[i]]))
##                             )                            
##         coll <- loadCollection(rp, names(df)[i])
##         w <- match(colnames(dfi), sapply(coll, setIdentifier))
##         nms <- data.frame(setNames=sapply(coll, setName)[w])
##         data <- list(list(dfi))
##         names(data[[1]]) <- names(df)[i]
##         se <- SummarizedExperiment(data,
##                                    rowData=conds,
##                                    colData=nms)
##         metadata(se) <- list(
##             set = colnames(res$details[[i]]),
##             background = bg,
##             sysinfo=Sys.info()
##         )
##         SEs[[i]] <- se
##     }       
## }

gep2pep <- function(geps, sets, min_size, max_size, parallel=FALSE,
                    pbar=TRUE, SGEmode=FALSE) {

    pathw <- sets
    genemat <- geps
    genes <- rownames(genemat)

    x <- list()
    sets <- sapply(unname(pathw), get, x="set")

    ES <- matrix(NA, length(pathw), ncol(genemat))
    PV <- matrix(NA, length(pathw), ncol(genemat))

    if(pbar)
        pb <- txtProgressBar()
    
    for(j in seq_len(ncol(genemat)))
    {
        if(pbar)
            setTxtProgressBar(pb, (j-1)/ncol(genemat))
        genematj <- genemat[,j]

        if(!SGEmode) {
            '%dobest%' <- if (parallel) get('%dopar%') else get('%do%')
            set <- NULL ## to cope with R CMD check NOTE
            gres <- foreach(set = sets,
                            .export=c("gsea","ks.test.2")) %dobest%
                {
                    where <- match(set, genes)
                    where <- where[!is.na(where)]
                    gsea(where, genematj, min_size, max_size)
                }

            x[[j]] <- gres
            
            if(all(is.na(sapply(gres, "get", x="ES"))))
                say(paste0("All NAs in PEP for profile: ",
                           colnames(genemat)[j]), "warning")

        } else { ## SGE mode
            nona <- genematj[!is.na(genematj)]
            pv <- .fastKSp(nona)
            es <- .fastKSES(nona)

            if(all(is.na(pv)) || all(is.na(es)))
                say(paste0("All NAs in PEP for profile: ",
                           colnames(genemat)[j]), "warning")
            PV[,j] <- pv
            ES[,j] <- es
        }
    }
        
    if(pbar) {
        setTxtProgressBar(pb, 1)
        close(pb)
    }

    if(!SGEmode)
        for(i in seq_len(ncol(genemat))){
            PV[,i] <- sapply(x[[i]], "get", x="p")
            ES[,i] <- sapply(x[[i]], "get", x="ES")
        }

    rownames(ES) <- rownames(PV) <-
        sapply(pathw, "get", x="id")
    colnames(ES) <- colnames(PV) <- colnames(genemat)
    
    return(list(ES=ES, PV=PV))
}


gsea <- function(S, ranks_list, min_size, max_size, alternative = "two.sided")
{
    S <- S[!(is.na(S))]
    S1 <- ranks_list[S]
    S2 <- ranks_list[-S]

    if(length(S1) < min_size || length(S1) > max_size ||
       length(S2) < min_size ||
       all(is.na(S1)) || all(is.na(S2)))
        return(list(ES=NA, p=NA))

    ks <- ks.test.2(S1, S2, alternative=alternative, maxCombSize=10^10)

    return(list(ES=ks$ES, p=ks$"p.value", edge=ks$edge));
}


storePEPs <- function(rp, db_id, peps, rawmode_suffix,
                      rawmode_outdir)
{
    rawmode <- !is.null(rawmode_suffix)
    
    if(rp$has(db_id) && !rawmode) {
        curmat <- rp$get(db_id)

        ## checking what's new and what's old
        curpeps <- colnames(curmat$ES)
        newpeps <- setdiff(colnames(peps$ES), curpeps)
        oldpeps <- intersect(colnames(peps$ES), curpeps)

        ## replacing existing PEPs
        curmat[["ES"]][, oldpeps] <- peps$ES[, oldpeps]
        curmat[["PV"]][, oldpeps] <- peps$PV[, oldpeps]

        ## adding new PEPs
        peps$ES <- cbind(curmat$ES, peps$ES[, newpeps, drop=FALSE])
        peps$PV <- cbind(curmat$PV, peps$PV[, newpeps, drop=FALSE])
    }

    if(!rawmode) {        
        say("Storing PEPs to the repository...")
        
        rp$put(peps, db_id,
               paste0("Pathway data for collection ", db_id,
                      ". It contains 2 matrices: 1 for enrichement scores ",
                      "(signed Kolmogorov Smirnov statistic) and one for ",
                      "the corresponding p-values."),
               c("gep2pep", "pep"), replace=TRUE,
               depends = paste0(db_id, "_sets"),
               prj = get_repo_prjname(rp))

        say("Storing condition information...")
        perts <- rp$get("conditions")
        perts[[db_id]] <- colnames(peps$ES)
        rp$set("conditions", perts)
    } else {
        fdb <- gsub("[^[:alnum:]]", "_", db_id)
        fname <- paste0(fdb,
                        "#",
                        format(rawmode_suffix, scientific=FALSE),
                        ".RDS")
        fname <- file.path(rawmode_outdir, fname)

        if(!dir.exists(rawmode_outdir))
            dir.create(rawmode_outdir)
        
        say(paste0("Storing PEPs to: ", fname))
        saveRDS(peps, fname)
    }

    say("Done.")
}

say <- function(txt, type="diagnostic", names=NULL) {
    ## can be error or warning
    msg <- paste0(
        "[",
        format(Sys.time(), format="%H:%M:%S"),
        "] ",
        txt,
        paste(names, collapse = ", ")
    )

    if(type=="error") {
        stop(msg, call.=FALSE)
    } else if (type=="warning") {
        warning(msg, call.=FALSE, immediate.=TRUE)
    } else message(msg)
}

ks.sign <- function (x, y)
{
    n.x <- as.double(length(x))
    n.y <- length(y)

    w <- c(x, y)        
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))

    return(sign(z[which.max(abs(z))]))
}


ks.test.2 <- function(x, y, ...) {
    ks <- ks.test(x, y, ...)
    ks[["ES"]] <- unname(ks.sign(x, y)*ks$statistic)
    return(ks)
}


rankPEPsByCols.SPV <- function(peps, rankingset="all")
{
    rankPEP <- function(PVs, ESs)
    {
        sorter <- abs(1-PVs)
        pos <- ESs > 0
        pos[is.na(pos)] <- FALSE
        sorter[pos] = -sorter[pos]        
        return(rank(sorter, ties.method = "random", na.last="keep"))
    }

    if(length(rankingset) == 1 && rankingset == "all")
        rankingset <- seq_len(nrow(peps[["ES"]]))
    
    PVs <- peps[["PV"]][rankingset, ]
    ESs <- peps[["ES"]][rankingset, ]
    x <- sapply(seq_len(ncol(PVs)), function(i) rankPEP(PVs[,i], ESs[,i]))
    colnames(x) <- colnames(PVs)
    rownames(x) <- rownames(PVs)
    attr(x, "nas") <- apply(x, 2, function(x) sum(is.na(x)))
    return(x)
}

rankPEPsByCols.NES <- function(peps, rankingset="all")
{
    if(length(rankingset) == 1 && rankingset == "all")
        rankingset <- seq_len(nrow(peps[["ES"]]))

    ESs <- t(scale(t(peps$ES)))
    attr(ESs, "scaled:center") <- NULL
    attr(ESs, "scaled:scale") <- NULL
    x <- apply(-ESs, 2, rank, ties.method = "random", na.last="keep")
    return(x)
}

rankPEPsByCols.ES <- function(peps, rankingset="all")
{
    if(length(rankingset) == 1 && rankingset == "all")
        rankingset <- seq_len(nrow(peps[["ES"]]))
  
    x <- apply(-peps$ES, 2, rank, ties.method = "random", na.last="keep")
    return(x)
}


rankPEPsByRows.ES <- function(peps)
{
    ESs <- peps[["ES"]]
    x <- t(apply(-ESs, 1, rank, ties.method = "random", na.last="keep"))
    attr(x, "nas") <- apply(x, 1, function(x) sum(is.na(x)))
    return(x)
}

attachInfo <- function(rp, results)
{
    type <- names(results)[[1]]
    dbs <- names(results[[type]])
    newres <- list()
    for(i in seq_along(results[[type]])) {
        db <- .loadCollection(rp, dbs[i])
        resmat <- results[[type]][[i]]

        if(type == "CondSEA") {
            modIDs <- rownames(resmat)
            newres[[i]] <- cbind(
                Pathway = sapply(db[modIDs], get, x="name"),
                Description = sapply(db[modIDs], get, x="desc"),
                results[[type]][[i]]
            )
            if(!is.null(results[["details"]]))
                newres[[i]] <- cbind(newres[[i]],
                                     results[["details"]][[i]])

        } else {
            res <- cbind(
                Condition = rownames(results[[type]][[i]]),
                results[[type]][[i]]
            )
            
            if(!is.null(results[["details"]])) {
                modIDs <- rownames(results[["details"]][[i]])
                modnames <- sapply(db[modIDs], get, x="name")
                details <- t(results[["details"]][[i]])

                colnames(details) <- paste0(modnames, " (",
                                            rownames(res$details[[i]]),
                                            ")")
                newres[[i]] <- cbind(res, details)
            } else {
                newres[[i]] <- res
            }
        }

    }
    names(newres) <- names(results[[type]])
    return(newres)
}


checkGEPsFormat <- function(geps)
{
    dims <- dim(geps)
    if(length(dims) != 2)
        say("GEPs must be a two-dimensional matrix", "error")
    if(dims[[1]]==1)
        say(paste("GEPs must have genes over rows,",
                  "while input matrix seems to have 1 row"),
            "error")

    gnames <- rownames(geps)
    if(is.null(gnames))
        say("GEPs must have rownames (as gene IDs)", "error")
    if(any(duplicated(gnames)))
        say("GEP rownames must be unique", "error")
    
    pnames <- colnames(geps)
    if(is.null(pnames))
        say("GEPs must have colnames (as condition IDs)", "error")
    if(any(duplicated(pnames)))
        say("GEP colnames must be unique", "error")
    
    not_unique <- apply(geps, 2, function(x) {
        any(duplicated(x))})
    mins <- apply(geps, 2, min, na.rm=TRUE)
    maxs <- apply(geps, 2, max, na.rm=TRUE)
    
    if(any(mins != 1) || any(maxs != dims[1]) || any(not_unique))
        say(paste("GEP columns must be ranks. Check",
                  "that each column is made of integer numbers from 1",
                  "to the number of rows."), "error")
    
}

get_repo_prjname <- function(rp) {
    prjname <- names(
        rp$entries()[sapply(lapply(rp$entries(), get, x="tags"),
                            `%in%`, x="#project")]
    )[[1]] ## take the first for robustness, should be only 1
}

## splits a flat list of pathways into sublists according to
## collections
pwList2pwStruct <- function(db, collids) {
    collids <- .makeCollectionIDs(db)
    colls <- unique(collids)
    out <- list()
    for(i in seq_along(colls))
        out[[colls[[i]]]] <- db[collids == colls[[i]]]
    return(out)
}


getPEPs <- function(rp, id) {
    if(! id %in% getCollections(rp))
        say(paste("Collection is not one of getCollections():", id), "error")

    if(! rp$has(id))
        say(paste("No PEP found for collection:", id), "error")
    
    return(rp$get(id))
}


checkSets <- function(rp, sets) {        
    coll_ids <- makeCollectionIDs(sets)
    ucoll_ids <- unique(coll_ids)
    
    off <- setdiff(ucoll_ids, getCollections(rp))
    if(length(off)>0)
        say(paste0("The following collections could not be found: ",
                   paste(off, collapse=", ")), "error")

    for(i in seq_along(ucoll_ids)) {
        sub <- sapply(sets[which(coll_ids == ucoll_ids[i])], setIdentifier)
        coll <- .loadCollection(rp, ucoll_ids[i])
        off <- setdiff(sub, names(coll))
        if(length(off) > 0)
            say(paste0("The following pathways could not be found ",
                       "in collection ", ucoll_ids[i], ": "), "error", off)
    }
}

convertFromGSetClass <- function(gsets) {
    res <- list()
    for(i in seq_along(gsets)) {
        set <- gsets[[i]]
        res[[setName(set)]] <- list(
            id = setIdentifier(set),
            name = setName(set),
            category = attributes(collectionType(set))$category,
            subcategory = attributes(collectionType(set))$subCategory,
            organism = organism(set),
            desc = description(set),
            set = geneIds(set)
        )
    }
    names(res) <- sapply(res, get, x="id")
    return(res)
}

.loadCollection <- function(rp, db) {
    thisdb <- rp$get(paste0(db, "_sets"))
    return(convertFromGSetClass(thisdb))
}

## loadCollection <- function(rp, db) {
##     thisdb <- rp$get(paste0(db, "_sets"))
##     return(thisdb)
## }


## this calls makeCollectionIDs with the old format
.makeCollectionIDs <- function(sets) {
    dbs <- sapply(sets, get, x="category")
    subdbs <- sapply(sets, get, x="subcategory")
    w <- which(subdbs=="" | is.na(subdbs))
    if(length(w)>0)
        subdbs[w] <- dbs[w]
    db_ids <- paste(dbs, subdbs, sep="_")
    return(db_ids)
}


#' Merge multiple PEPs to build a repository of consensus PEPs
#' @inheritParams dummyFunction
#' @param rpIn_path path to existing gep2pep repository
#' @param rpOut_path path where the new merged repository will be
#' created
#' @param mergestr a named list of character vectors, each one
#' including a set of PEP names. For each list entry, a consensus PEP
#' will be created and assigned the entry name.
#' @param progressBar if TRUE, show a progress bar
#' @details The merging is performed as follows. Given N PEPs, the
#' corresponding consensus PEP will get as enrichement score the
#' average enrichment scores of the N PEPs, and as p-value the
#' composition of the N PEP p-values by Fisher's method.
#' @return Nothing, used for side effects.
#' @export
#' @examples
#'
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#' rp <- createRepository(repo_path, db)
#' geps <- loadSampleGEP()
#' buildPEPs(rp, geps)
#'
#' mergestr <- list(
#'     che_iso = c("(+)_chelidonine", "(+)_isoprenaline"),
#'     cat_mk = c("(+/_)_catechin", "(_)_mk_801")
#' )
#'
#' merged_path <- file.path(tempdir(), "gep2pepTempMerged")
#'
#' createMergedRepository(repo_path, merged_path, mergestr)
#'
#' unlink(repo_path, TRUE)
#' unlink(merged_path, TRUE)
createMergedRepository <- function(rpIn_path, rpOut_path, mergestr,
                                   progressBar=TRUE, collections="all") {
  ## the following are hidden parameters
  mergeFunc <- "mean"
  parallel <- FALSE

    if(!missing(mergestr)) {
      if(!is.list(mergestr)) {
        say(
            paste("mergestr parameter must be a list of character",
                  "(column names) or numbers (column indices)"),
          "error"
          )
      }        
        if(is.null(names(mergestr)) || any(duplicated(names(mergestr)))) {
            say(paste("mergestr entries must be uniquely named"), "error")
        }
    }

    ## if(any(duplicated(unlist(mergestr)))) {
    ##     say("mergestr includes duplicated entries", "error")
    ## }

    
    rpin <- openRepository(rpIn_path)
    if(!file.exists(rpOut_path)) {
        rpout <- createRepository(rpOut_path, NULL)
    } else {
        say("Output repository already exists")
        rpout <- openRepository(rpOut_path)
        if(missing(mergestr)) {
            say("Loading merge structure from the repository")
            mergestr <- rpout$get("merging structure")
        }
    }

    colls <- getCollections(rpin)
    if(!(length(collections)==1 && collections=="all")) {
        colls <- intersect(colls, collections)
    }

    if(length(colls)>0) {
        for(i in 1:length(colls)) {
            say(paste0("Copying sets for collection: ", colls[i]))
            rpin$copy(rpout, paste0(colls[i], "_sets"))        
            say("Merging PEPs")
            mergedPEPs <- .mergePEPs(rpin, colls[i], mergestr,
                                     progressBar, mergeFunc, parallel)
            suppressMessages(storePEPs(rpout, colls[i], mergedPEPs, NULL, NULL))
            say("Merged PEPs stored")
        }
        say("Storing merging structure")
        
        if(rpout$has("merging structure")) {
            say("Existing merging structure will not be overwritten",
                "warning")
        } else {
            rpout$put(mergestr, "merging structure",
                      "list of single PEPs used to make each merged PEP")
        }
        
        say("Done.")
    } else say("None of the specified collections is available", "warning")
}





.mergePEPs <- function(rpin, coll, mergestr, progressBar,
                       mergeFunc = "mean", parallel=FALSE) {
    
    if(!is.character(unlist(mergestr)))
        say("mergestr entries must be of class character", "error")
    ish5 <- "#hdf5" %in% rpin$tags(coll)
    
    nrows <- length(loadCollection(rpin, coll))
    
    outpeps <- list(
        ES=matrix(NA, nrows, length(mergestr)),
        PV=matrix(NA, nrows, length(mergestr))
    )

    if(progressBar)
        pb <- txtProgressBar()

    if(mergeFunc == "mean") {
        if(parallel)
            say("Parallelization will not be used for mean method")
        if(!ish5)
            allpeps <- .loadPEPs(rpin, coll)

        for(i in 1:length(mergestr)) {
            if(progressBar)
                setTxtProgressBar(pb, i/length(mergestr))

            if(ish5) {
                peps <- .loadPEPs(rpin, coll, mergestr[[i]])
            } else {
                peps <- list(ES=allpeps$ES[, mergestr[[i]], drop=FALSE],
                             PV=allpeps$PV[, mergestr[[i]], drop=FALSE])
            }

            mpeps <- list(
                ES = apply(peps$ES, 1, mean),
                PV = apply(peps$PV, 1, function(ps) {
                    S=-2*sum(log(ps))
                    cump <- pchisq(S, 2*length(ps))
                    return(cump)
                })
            )

            outpeps$ES[,i] <- mpeps$ES
            outpeps$PV[,i] <- mpeps$PV
        }

        if(ish5) {
            rownames(outpeps$ES) <-
                rownames(outpeps$PV) <-
                rownames(peps$ES)
        } else {
            rownames(outpeps$ES) <-
                rownames(outpeps$PV) <-
                rownames(allpeps$ES)
        }

    } else if(mergeFunc == "CondSEA") {
        for(i in 1:length(mergestr)) {
            if(progressBar)
                setTxtProgressBar(pb, i/length(mergestr))
        
            ranked <- suppressMessages(
                .getPEPRanks(rpin, coll,
                             unlist(mergestr),
                             subcols = mergestr[[i]],
                             rankingFun = rankPEPsByRows.ES
                             )
                )

            ks <- .doPSEA(ranked, length(unlist(mergestr)), parallel)
                                  
            outpeps$PV[,i] <- sapply(ks, get, x="p.value")
            outpeps$ES[,i] <- sapply(ks, get, x="ES")            
        }

        rownames(outpeps$ES) <-
            rownames(outpeps$PV) <-
            rownames(ranked)

    } else say("Wrong merging method", "error")

    if(progressBar)
        close(pb)    
    
    colnames(outpeps$ES) <-
        colnames(outpeps$PV) <-
        names(mergestr)

    return(outpeps)
}


.buildGene2SetIndex <- function(rp, collection_name)
{
    collection <- loadCollection(rp, collection_name)
    allsets <- sapply(collection, geneIds)
    setids <- sapply(collection, setIdentifier)
    allgenes <- unique(unlist(allsets))
    
    res <- list()
    pb <- txtProgressBar()
    for(i in 1:length(allgenes)) {
        setTxtProgressBar(pb, i/length(allgenes))
        w <- sapply(allsets, `%in%`, x=allgenes[i])
        res[[i]] <- c(
            gene = allgenes[i],
            sets = paste(setids[w], collapse="; ")
        )
    }
    return(res)
}

.exportCollectionsInfo <- function(rp, outfile) {
    colls <- getCollections(rp)

    dbs <- ids <- names <- vector("character")
    for(i in 1:length(colls)) {
        coll <- loadCollection(rp, colls[i])
        dbs <- c(dbs, rep(colls[i], length(coll)))
        ids <- c(ids, sapply(coll, setIdentifier))
        ## names <- c(names, paste0("[", colls[i], "] ", sapply(coll, setName)))
        names <- c(names, sapply(coll, setName))
    }

    return(cbind(dbs, ids, names))
}

## This function is copy-pasted from the utils packages, as it was not
## exported and the check command would complain about using the `:::`
## operator. Original name: utils:::.format.object_size
.format.object_size <- function (x, units = "b",standard = "auto",
                                 digits = 1L, ...)
{
    known_bases <- c(legacy = 1024, IEC = 1024, SI = 1000)
    known_units <- list(SI = c("B", "kB", "MB", "GB", "TB", "PB", 
        "EB", "ZB", "YB"), IEC = c("B", "KiB", "MiB", "GiB", 
        "TiB", "PiB", "EiB", "ZiB", "YiB"), legacy = c("b", "Kb", 
        "Mb", "Gb", "Tb", "Pb"), LEGACY = c("B", "KB", "MB", 
        "GB", "TB", "PB"))
    units <- match.arg(units, c("auto", unique(unlist(known_units), 
        use.names = FALSE)))
    standard <- match.arg(standard, c("auto", names(known_bases)))
    if (standard == "auto") {
        standard <- "legacy"
        if (units != "auto") {
            if (grepl("iB$", units)) 
                standard <- "IEC"
            else if (grepl("b$", units)) 
                standard <- "legacy"
            else if (units == "kB") 
                stop("For SI units, specify 'standard = \"SI\"'")
        }
    }
    base <- known_bases[[standard]]
    units_map <- known_units[[standard]]
    if (units == "auto") {
        power <- if (x <= 0) 
            0L
        else min(as.integer(log(x, base = base)), length(units_map) - 
            1L)
    }
    else {
        power <- match(toupper(units), toupper(units_map)) - 
            1L
        if (is.na(power)) 
            stop(gettextf("Unit \"%s\" is not part of standard \"%s\"", 
                sQuote(units), sQuote(standard)), domain = NA)
    }
    unit <- units_map[power + 1L]
    if (power == 0 && standard == "legacy") 
        unit <- "bytes"
    paste(round(x/base^power, digits = digits), unit)
}

#' Adds a collection of single-gene psuedo-sets.
#'
#' This function can be used to add single-gene (as opposed to
#' pathway) -based collections. Sets including a single gene don't need
#' to go through normal Kolmogorov-Smirnov statistic computation and
#' are treated differently for performance.
#'  
#' @inheritParams dummyFunction
#' @param organism Character vector used to annotate the created
#'     sets. "Homo Sapiens" by default.
#' @return Nothing, used for side effects.
#' @seealso buildPEPs
#' @details
#'
#' Enrichment Scores and p-values for sets including a single gene are
#' computed with dedicated (fast) routines. Although a statistic based
#' on a single gene is not efficient per se, it is useful to have data
#' in the same format as pathway-based profiles. \code{buildPEPs}
#' internally calls single gene dedicated routines whenever a gene set
#' collection is tagged (see repo function \code{tag}) with "SGE"
#' ("Single Gene Expression"), which is done automatically by
#' \code{addSingleGeneSets}. In that case, the \code{min_size}
#' parameter is ignored.
#' 
#' @export
#' @examples
#' 
#' db <- loadSamplePWS()
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#'
#' ## The following will create PEPs in 2 separate files
#' geps <- loadSampleGEP()
#' addSingleGeneSets(rp, rownames(geps))
#'
#' unlink(repo_path, TRUE)
#' 
addSingleGeneSets <- function(rp, genes, organism = "Homo Sapiens")
{
gs <- list()
descr <- "Pseudo-gene-set containing only gene "
for(i in 1:length(genes)) {
    g <- genes[i]
    gs[[i]] <- GeneSet(g,
                       shortDescription = paste0(descr, g),
                       longDescription = paste0(descr, g),
                       setName = g,
                       setIdentifier = paste0("SGE-", g),
                       organism = organism,
                       collectionType = CategorizedCollection(
                           category="SGE",
                           subCategory="SGE"
                       ))
    }
sge <- GeneSetCollection(gs)
rp$put(sge, "SGE_sets", "Pathway information for collection SGE",
       tags=c("gep2pep", "sets", "SGE"))
    
}

.fastKSp <- function(v) {
    ## converts a vector of ranks to corresponding KS p-values
    ## used for single-gene sets
    n <- length(v)-sum(is.na(v))
    up <- seq(1, floor(n/2))/(n/2)
    if(length(v)%%2==1)
        ret <- c(up, 1, rev(up)) else ret <- c(up, rev(up))
    return(ret[v])
}

.fastKSES <- function(v) {
    ## converts a vector of ranks to corresponding KS Enrichment
    ## Scores. Used for single-gene sets
    n <- length(v)
    up <- -1/(n-1)*(1:ceiling(n/2)-1) + 1
    if(length(v)%%2==1)
        ret <- c(up[-length(up)], -rev(up)) else ret <- c(up, -rev(up))
    return(ret[v])
}


exportPEPs <- function(rp, PEPname, outfile=paste0(PEPname, ".xls"),
                       collections="all")
{
    colls <- getCollections(rp)

    if(!(length(collections)==1 && collections=="all"))
        colls <- intersect(collections, colls)

    newres <- list()
    
    for(i in 1:length(colls)) {
        db <- .loadCollection(rp, colls[i])
        PEPs <- .loadPEPs(rp, colls[i], PEPname)
        ord <- order(PEPs$PV)
        PEPs$ES <- PEPs$ES[ord,,drop=F]
        PEPs$PV <- PEPs$PV[ord,,drop=F]
        
        modIDs <- rownames(PEPs$ES)

        newres[[i]] <- data.frame(
            Pathway = sapply(db[modIDs], get, x="name"),
            Description = sapply(db[modIDs], get, x="desc"),
            ES = PEPs$ES[,PEPname],
            p.value = PEPs$PV[,PEPname],
            FDR = p.adjust(PEPs$PV, "fdr", n=sum(!is.na(PEPs$PV))),
            check.names = FALSE
        )
    }
    names(newres) <- gsub(":", "_", colls)

    WriteXLS(newres, outfile, BoldHeaderRow = TRUE, FreezeRow = 1,
             AutoFilter = TRUE)
}



do_and_show_gsea <- function(gep, set) {
    setgenes <- geneIds(set)    
    
    x <- gep[setgenes]
    x <- x[!is.na(x)]
    y <- setdiff((1:sum(!is.na(gep))), x)
    n.x <- as.double(length(x))
    n.y <- length(y)
    w <- c(x, y)        
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))

    m <- which.max(abs(z))
    s <- sign(z[m])    

    ks <- ks.test.2(x, y, maxCombSize=10^10)
    ES <- unname(s*ks$statistic)
    p <- ks$p.value

    if(s > 0) {
        leadset <- x[x <= m]
    } else leadset <- x[x >= m]

    fname <- paste(drug, dbi, pwi, sep="_")

    say("Creating PNG device")
    pngf <- paste0(outdir, "/", fname, ".png")
    library(Cairo)
    CairoPNG(pngf, 480*2, 480)
    ##png(pngf, 480*2, 480)

    say("Creating plot")
    plot(z, type="l", ylab="ES", xlab="Ranks",
         main=paste0(drugname, "\n", G_allsubdbs[dbi], " - ", pwname),
         lwd=1.5, col="gray")
    points(m, z[m], pch=21, cex=2)
    points(x, z[x], cex=.5)
    if(s>0) {
        points(x[x <= m], z[x[x <= m]], bg="green", col="green", pch=21)
    } else {
        points(x[x >= m], z[x[x >= m]], bg="red", col="red", pch=21)
    }    
    dev.off()


    ## say("Creating table")
    ## tabf <- paste0(outdir, "/", fname, ".txt")
    ## tab <- cbind(profile[set,,drop=F], lead=F)
    ## tab[names(leadset), 2] <- T
    ## tab <- tab[order(tab[,1]),]
    ## tab <- rbind(stats = c(ES,p), tab)
    ## say("Saving table")
    ## write.table(tab, tabf, quote=F, col.names=F)
    ## say("Table saved")
}
