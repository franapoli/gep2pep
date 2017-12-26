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
#' [2] Napolitano F. et al, Drug-set enrichment analysis: a novel tool
#'     to investigate drug mode of action. Bioinformatics 32, 235-241
#'     (2016).
#' [3] Napolitano F. et al, gene2drug: a Computational Tool for
#'     Pathway-based Rational Drug Repositioning, bioRxiv (2017)
#'     192005; doi: https://doi.org/10.1101/192005
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

## For gene sets management
#' @import GSEABase

#' @importFrom foreach foreach %dopar% %do%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats ks.test
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
#' @param category A character defining the main category that the gene
#'     set belongs to.
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
#' dim(geps)
#' ## [1] 485   5
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
#' geps
#' ## GeneSetCollection
#' ##   names: AAANWWTGC_UNKNOWN, AAAYRNCTG_UNKNOWN, ..., MORF_BUB3 (30 total)
#' ##   unique identifiers: MEF2C, ATP1B1, ..., SLBP (5778 total)
#' ##   types in collection:
#' ##     geneIdType: SymbolIdentifier (1 total)
#' ##     collectionType: BroadCollection (1 total)
#'
loadSamplePWS <- function() {
    return(
        readRDS(system.file("extdata", "testgmd.RDS",
                            package="gep2pep"))
    )
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
#' #' ## A small sample of the MSigDB as imported by importMSigDB.xml is
#' ## included in gep2pep. The following creates (and deletes) a
#' ## gep2pep repository.
#'
#' db_sample <- loadSamplePWS()
#' db_sample <- as.CategorizedCollection(db_sample)
#' 
#' ## The function can also be used to create arbitrary gene set
#' ## collections specifying the categories and subcategories once for
#' ## all the sets:
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

#' Compose raw PEPs into a matrix
#'
#' @return Sample gene expression data
#' @export
#' @examples
importFromRawMode <- function(rp, path=file.path(rp$root(), "raw")) {
  allfiles <- list.files(path)
  say(paste0("Found ", length(allfiles), " raw files."))
  ids <- as.numeric(gsub(".+_|[.]RDS", "", allfiles))
  say(paste0("Found ", length(unique(ids)), " unique IDs."))
  dbs <- gsub("_.[.]RDS", "", outfiles)
  say(paste0("Found ", length(unique(dbs)), " collection names."))

  library(rhdf5)

  ## extracting chunk size
  fname <- paste0(dbs[1], "_", min(ids), ".RDS")
  Nchunk <- ncol(readRDS(file.path(path, fname))$ES)
  fname <- paste0(dbs[1], "_", max(ids), ".RDS")
  NlastChunk <- ncol(readRDS(file.path(path, fname))$ES)
  Ncol <- Nchunk*(length(unique(ids))-1) + NlastChunk
  say(paste0("Expected profiles: ", Ncol, " in ", Nchunk, " chunks."))
  
  for(i in 1:length(unique(dbs))) {
      dbi <- unique(dbs)[i]

      if(rp$has(dbi)) {
          say(paste0("A collection named ", dbi,
                     " is already in the repository ",
                     "and will be skipped in raw mode"), "warning")
          next
      }
      
      say(paste0("Working on collection: ", dbi))

      fl <- tempfile()
      say(paste0("Using temporary file: ", fl))
      print(fl)
      h5createFile(fl)

      fname <- paste0(dbi, "_", min(ids), ".RDS")
      x <- readRDS(file.path(path, fname))$ES
      Nrow <- nrow(x)
      h5createDataset(fl, "ES", c(Nrow, Ncol))
      h5createDataset(fl, "PV", c(Nrow, Ncol))
      h5createDataset(fl, "rownames", Nrow, storage.mode="character", size=256)
      h5createDataset(fl, "colnames", Ncol, storage.mode="character", size=256)
      h5write(rownames(x), fl, "rownames")

      pb <- txtProgressBar()
      for(j in 1:length(unique(ids))) {
          fname <- paste0(dbi, "_", unique(ids)[j], ".RDS")
          x <- readRDS(file.path(path, fname))
          startCol <- (j-1)*Nchunk+1
          h5write(x$ES, fl, "ES", start=c(1,startCol))
          h5write(x$PV, fl, "PV", start=c(1,startCol))
          h5write(colnames(x$ES), fl, "colnames", start=startCol)
          setTxtProgressBar(pb, j/length(unique(ids)))
      }

      say("\nStoring into the repository...")
      rp$put(fl, dbi, asattach=T, tags=c("pep", "#hdf5"))
      #rp$untag(dbi, "hide") ## <- repo bug, gives error
      say("Clearing temporary file...")
      file.remove(fl)
      say("Done.")
  }
}

#' Imports pathways data from an MSigDB XML file.
#'
#' Creates a \code{GeneSetCollection} object using the XML
#' distribution of the MSigDB (see references). The returned object
#' can be passed to \code{createRepository}.
#'
#' @param fname Path to an XML file downloaded from MSigDB.
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
#' db_sample <- as.CategorizedCollection(db_sample)
#' 
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#' rp <- createRepository(repo_path, db_sample)
#'
#' ## removing temporary repository
#' unlink(repo_path, TRUE)
#' @export
importMSigDB.xml <- function(fname) {

    say("Loading gene sets...")

    result = tryCatch({
        sets <- getBroadSets(fname)
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
        for(i in seq_len(nrow(msigDB))) {
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
dummyFunction <- function(rp, rp_peps, collections) {}


#' Check an existyng repository for consistency
#'
#' Check both repository data consistency (see \code{repo_check} from
#' the \code{repo} package) and specific gep2pep data consistency.
#' @inheritParams dummyFunction
#' @return Nothing.
#' @examples
#' db <- loadSamplePWS()
#' db <- as.CategorizedCollection(db)
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
    say("Checing for gep2pep data consistency...")
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

        say("Summary of common conditions across collections:")
        out <- outer(perts, perts,
                     Vectorize(
                         function(a,b) length(intersect(a,b))
                     ))
        print(out)
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
#' db <- as.CategorizedCollection(db)
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
#' db <- as.CategorizedCollection(db)
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
#' db <- as.CategorizedCollection(db)
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
#' db <- as.CategorizedCollection(db)
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
#' db <- as.CategorizedCollection(db)
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
    if(!is(sets, "GeneSetCollection"))
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
#' db <- as.CategorizedCollection(db)
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
#' db <- as.CategorizedCollection(db)
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
    subdbs[subdbs==""] <- dbs[subdbs==""]
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
#' @param parallel If TRUE, gene sets will be processed in
#'     parallel. Requires a parallel backend.
#' @param replace_existing What to do if PEPs, identified by column
#'     names of \code{geps} are already present in the repository. If
#'     set to TRUE, they will be replaced, otherwise they will be
#'     skipped and only PEPs of new conditions will be computed and
#'     added. Either ways, will throw a warning.
#' @param progress_bar If set to TRUE (default) will show a progress
#'     bar updated after coversion of each column of \code{geps}.
#' @param rawmode_suffix A character vector to be appended to files
#'     produced in raw mode (see details). All non-alphanumeric
#'     characters will be replaced by an underscore (_). If set to
#'     NULL (default), raw mode is turned off.
#' @param rawmode_outdir A charater vector specifying the destination
#'     path for files produced in raw mode (by the fault it is
#'     ROOT/raw, where ROOT is the root of the repository). Ignored if
#'     \code{rawmode_suffix} is NULL.
#' @return Nothing. The computed PEPs will be available in the
#'     repository.
#' @seealso buildPEPs
#' @details By deault, output is written to the repository as new
#'     items named using the collection name. However, it is possible
#'     to avoid the repository and write the output to regular files
#'     turning 'raw mode' on through the \code{rawmode_suffix} and
#'     \code{rawmode_outdir} parameters. This is particuarly useful
#'     when dealing with very large corpora of GEPs, and conversions
#'     are split into independent jobs submitted to a scheduler. At
#'     the end, the data will need to be reconstructed and put into
#'     the repository manually in order to perform \code{CondSEA} or
#'     \code{PathSEA} analysis.
#' @examples
#' db <- loadSamplePWS()
#' db <- as.CategorizedCollection(db)
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
buildPEPs <- function(rp, geps, parallel=FALSE, collections="all",
                      replace_existing=FALSE, progress_bar=TRUE,
                      rawmode_suffix=NULL,
                      rawmode_outdir=file.path(rp$root(), "raw"))   
{
    checkGEPsFormat(geps)
    perts <- rp$get("conditions")
    rawmode <- !is.null(rawmode_suffix)
    
    if(length(collections) == 1 && collections == "all") {
        dbs <- getCollections(rp)
        } else dbs <- collections

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
                say(paste0("Existing PEPs found, ",
                           "but this will be ignored in rawmode: ",
                           paste(oldpeps, collapse=", ")))
                newpeps <- colnames(geps)
            } else {
                msg <- paste0("Existing PEPs will be replaced: ",
                              paste(oldpeps, collapse=", "))
                if(!replace_existing)
                    msg <- gsub("replaced", "skipped", msg)
                say(msg, type="warning")
            }
        }

        if(length(newpeps) > 0) {
            gepsi <- geps[, newpeps, drop=FALSE]
            thisdb <- .loadCollection(rp, dbs[i])
            peps <- gep2pep(gepsi, thisdb, parallel, progress_bar)
            storePEPs(rp, dbs[i], peps, rawmode_suffix,
                      rawmode_outdir)
        }
    }
}



#' Export CondSEA or PathSEA results to XLS format
#'
#' The XLS output includes the full CondSEA or PathSEA results,
#' together with additional gene set information for the CondSEA.
#'
#' @inheritParams dummyFunction
#' @param results The output of \code{CondSEA} or \code{PathSEA}.
#' @param outname Name of the XLS file to be created.
#' @return Nothing.
#' @seealso CondSEA, PathSEA
#' @examples
#'
#' db <- loadSamplePWS()
#' db <- as.CategorizedCollection(db)
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
#' db <- as.CategorizedCollection(db)
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
#' db <- as.CategorizedCollection(db)
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
    return(colnames(rp$get(coll)$ES))
}

.loadPEPs <- function(rp, coll, subset) {
    peps <- list(
        ES = rp$get(coll)$ES[, subset, drop=F],
        PV = rp$get(coll)$PV[, subset, drop=F]
    )
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
#' @return A list of 2, by names "CondSEA" and "details". The
#'     "CondSEA" entry is a 2-columns matrix including ESs and
#'     p-values (see details) for each pathway database and
#'     condition. The "details" entry reports the rank of each
#'     condition in \code{pgset} for each pathway.
#' @details For each pathway, all conditions are ranked by how much
#'     they dysregulate it (from the most UP-regulating to the most
#'     DOWN-regulating). Then, a Kolmogorov-Smirnov (KS) test is
#'     performed to compare the ranks assigned to conditions in
#'     \code{pgset} against the ranks assigned to conditions in
#'     \code{bgset}. A positive (negative) Enrichment Score (ES) of
#'     the KS test indicates whether each pathway is UP- (DOWN-)
#'     regulated by \code{pgset} as compared to \code{bgset}. A
#'     p-value is associated to the ES.
#'
#'     When PEPs are obtained from drug-induced gene expression
#'     profiles, \code{PathSEA} can be used to perform Drug-Set
#'     Enrichment Analysis [1].
#' @seealso getResults, getDetails
#' @references
#'
#' [1] Napolitano F. et al, Drug-set enrichment analysis: a
#'     novel tool to investigate drug mode of action. Bioinformatics
#'     32, 235-241 (2016).
#'
#' @examples
#' db <- loadSamplePWS()
#' db <- as.CategorizedCollection(db)
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
CondSEA <- function(rp_peps, pgset, bgset="all", collections="all",
                    details=TRUE)
{
    dbs <- collections
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

        allperts <- .loadPerts(rp, dbs[i])

        if(length(bgset) == 1 && bgset=="all")
            bgset <- allperts
        
        if(length(intersect(pgset, bgset))>0) {
            bgset <- setdiff(bgset, pgset)
            say("Common conditions removed from bgset")
        }
        
        rankingset <- c(bgset, pgset)

        if(!all(rankingset %in% allperts))
            say(paste("The following conditions could not be found:",
                      paste(
                          setdiff(rankingset, colnames(peps$ES)),
                          collapse = ", ")), "error")

        peps <- .loadPEPs(rp_peps, dbs[i], rankingset)
        say(paste0("Row-ranking collection"))
        ranked <- rankPEPsByRows(peps)
        say(paste0("Computing enrichments"))
        
        ks <- apply(ranked, 1, function(row) {
            inset <- row[pgset]
            inset <- inset[!is.na(inset)]
            outset <- row[bgset]
            outset <- outset[!is.na(outset)]
            if(length(inset)>1 && length(outset)>1) {
                res <- ks.test.2(row[pgset], row[bgset], maxCombSize=10^10)
            } else res <- list(ES=NA, p.value=NA)
            return(res)
        })
        PVs <- sapply(ks, get, x="p.value")
        sorter <- order(PVs)
        res[[dbs[i]]] <- data.frame(ES = sapply(ks, get, x="ES"),
                                    PV = PVs)[sorter, ]
        if(details)
            thedetails[[dbs[i]]] <- ranked[sorter, pgset]
        say("done")
    }

    return(list(CondSEA=res, details=thedetails))
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
#' @param details If TRUE (default) details will be reported for each
#'     condition in \code{pgset}.
#' @return A list of 2, by names "PathSEA" and "details". The
#'     "PathSEA" entry is a 2-columns matrix including ESs and
#'     p-values for each collection and condition. The "details"
#'     entry reports the rank of each pathway in \code{pathways} for
#'     each condition.
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
#' @seealso getResults, getDetails
#' @references [1] Napolitano F. et al, gene2drug: a Computational
#'     Tool for Pathway-based Rational Drug Repositioning, bioRxiv
#'     (2017) 192005; doi: https://doi.org/10.1101/192005
#' @examples
#' library(GSEABase)
#'
#' db <- loadSamplePWS()
#' db <- as.CategorizedCollection(db)
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
PathSEA <- function(rp_peps, pathways, bgsets="all", collections="all",
                    details=TRUE)
{
    checkSets(rp_peps, pathways)
    
    pathways <- convertFromGSetClass(pathways)
    
    pathways <- pwList2pwStruct(pathways)

    for(i in seq_along(pathways))
        if(!rp_peps$has(names(pathways)[i]))
            say("Cold not find PEPs: ", "error", names(pathways)[i])

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
        peps <- .loadPEPs(rp_peps$get, collections[i])
        notok <- rankingset[rankingset %in% rownames(peps)]
        if(length(notok)>0)
            say(paste0("Pathway set ids not found in ", collections[i], ": ",
                       paste(notok, collapse=", ")), "error")

        say(paste0("Column-ranking collection"))
        ranked <- rankPEPsByCols(peps, rankingset)
        say(paste0("Computing enrichments"))

        ks <- apply(ranked, 2, function(col) ks.test.2(col[gmd], col[bgset]))

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
#' @param gene A gene identifier of the same type as that used to
#'     create the repository.
#' @return A database of pathways suitable as input to
#'     \code{\link{PathSEA}}.
#' @seealso createRepository, PathSEA
#' @examples
#' db <- loadSamplePWS()
#' db <- as.CategorizedCollection(db)
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
gene2pathways <- function(rp, gene)
{
    dbs <- getCollections(rp)
    mods <- list()
    for(i in seq_along(dbs)) {
        db <- loadCollection(rp, dbs[i])
        w <- sapply(db, function(s) gene %in% geneIds(s))
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

gep2pep <- function(geps, sets, parallel=FALSE, pbar=TRUE) {

    pathw <- sets
    genemat <- geps
    genes <- rownames(genemat)

    x <- list()
    sets <- sapply(unname(pathw), get, x="set")

    if(pbar)
        pb <- txtProgressBar()
    for(j in seq_len(ncol(genemat)))
    {
        if(pbar)
            setTxtProgressBar(pb, (j-1)/ncol(genemat))
        genematj <- genemat[,j]

        '%dobest%' <- if (parallel) get('%dopar%') else get('%do%')
        set <- NULL ## to cope with R CMD check NOTE
        gres <- foreach(set = sets,
                        .export=c("gsea","ks.test.2")) %dobest%
        {
            where <- match(set, genes)
            where <- where[!is.na(where)]
            gsea(where, genematj, FALSE)
        }
        x[[j]] <- gres
    }
    if(pbar) {
        setTxtProgressBar(pb, 1)
        close(pb)
    }

    ES <- matrix(NA, length(pathw), ncol(genemat))
    PV <- matrix(NA, length(pathw), ncol(genemat))
    
    for(i in seq_len(ncol(genemat))){
        PV[,i] <- sapply(x[[i]], "get", x="p")
        ES[,i] <- sapply(x[[i]], "get", x="ES")
    }

    rownames(ES) <- rownames(PV) <-
        sapply(pathw, "get", x="id")
    colnames(ES) <- colnames(PV) <- colnames(genemat)
    
    return(list(ES=ES, PV=PV))
}


gsea <- function(S, ranks_list, check=FALSE, alternative = "two.sided")
{
    S <- S[!(is.na(S))]
    S1 <- ranks_list[S]
    S2 <- ranks_list[-S]

    if(length(S1)<1 || length(S2)<1 || all(is.na(S1)) || all(is.na(S2)))
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
    
    say("Storing pathway expression profiles")

    if(!rawmode) {
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
        rawmode_suffix <- gsub("[^[:alnum:]]", "_", rawmode_suffix)
        fdb <- gsub("[^[:alnum:]]", "_", db_id)
        fname <- paste0(fdb, rawmode_suffix, ".RDS")
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


rankPEPsByCols <- function(peps, rankingset="all")
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
    return(x)
}


rankPEPsByRows <- function(peps)
{
    ESs <- peps[["ES"]]
    x <- t(apply(-ESs, 1, rank, ties.method = "random", na.last="keep"))
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
                results[[type]][[i]],
                results[["details"]][[i]]
            )
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
                  "that each column is made of numbers from 1",
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
    subdbs[subdbs==""] <- dbs[subdbs==""]
    db_ids <- paste(dbs, subdbs, sep="_")
    return(db_ids)
}

.extractWorkingPEPs <- function(rp, coll, fgset, bgset) {
    ishdf5 <- "#rhdf5" %in% rp$tags(coll)

    
}
