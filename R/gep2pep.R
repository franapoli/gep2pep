
#' gep2pep: creation and analysis of Pathway Expression Profiles
#'
#' Pathway Expression Profiles (PEPs) are based on the collective
#' expression of genes belonging to a collection of pathways (or
#' generic gene sets) as opposed to individual genes. \code{gep2pep}
#' supports the conversion of gene expression profiles (GEPs) to PEPs
#' and performs enrichment analysis of pathways or perturbagens.
#'
#' @details
#'
#' \code{gep2pep} creates a local repository of gene sets, which can
#' also be imported from the MSigDB database. The local repository is
#' in the \code{repo} format. When a GEP, defined as a ranked list of
#' genes, is passed to \code{\link{buildPEPs}}, the stored database of
#' pathways is used to convert the GEP to a PEP and permanently store
#' the latter.
#'
#' One type of analysis that can be performed on PEPs and that is
#' directly supported by \code{gep2pep} is the Drug-Set Enrichment
#' Analysis (DSEA, see reference section). It finds pathways that
#' are consistently dysregulated by a set of drugs, as opposed to a
#' background of other drugs. Of course PEPs may refer to
#' non-pharmacological perturbagens (genetic perturbations, disease
#' states, etc.) for analogous analyses. See \code{\link{PertSEA}}
#' function.
#'
#' A complementary approach is that of finding perturbagens that
#' consistently dysregulate a set of pathways. This is the
#' pathway-based version of the Gene Set Enrichment Analysis
#' (GSEA). See \code{\link{PathSEA}}.
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
#'   \item{perturbagen: }{any experimental condition (drug, gene
#'   knock-out, disease state, cell type, etc.) characterized by an
#'   induced GEP and therefore a PEP.}
#'
#'   \item{gep2pep repository: }{a pathway database and possibly a
#'   related database of PEPs as created by the \code{gep2pep}
#'   package. It is implemented in \code{repo} format.}
#'
#' }
#' @references Napolitano F. et al, Drug-set enrichment analysis: a
#'     novel tool to investigate drug mode of action. Bioinformatics
#'     32, 235-241 (2016).
#' 
#' @docType package
#' @name gep2pep-package
#' @author Francesco Napolitano \email{franapoli@@gmail.com}
#' @aliases gep2pep
NULL


## repo is for storage of repositories
#' @import repo
## XML is to import MSigDB data
#' @import XML
## foreach is for easy parallelization support
#' @import foreach
## utils is for txtProgressBar
#' @import utils
## stats is for ks.test
#' @import stats
NULL


#' Dummy function for parameter inheritance
#' @param rp A repository created by \code{\link{createRepository}}.
#' @param rp_peps A repository created with
#'     \code{\link{createRepository}}, and containing PEPs created
#'     with \code{\link{buildPEPs}}.
#' @param collections A subset of the collection names returned by
#'     \code{getCollections}. If set to "all" (default), all the
#'     collections in \code{rp} will be used.
#' @return Nothing
dummyFunction <- function(rp, rp_peps, collections) {}


#' Imports pathways data from an MSigDB XML file.
#'
#' Creates a \code{gep2mep} repository of pathway data using the XML
#' distribution of the MSigDB (see references).
#'
#' @param fname Path to an XML file downloaded from MSigDB.
#' @return A list of pathway entries (see
#'     \code{\link{createRepository}}).
#' @references
#'     \url{http://software.broadinstitute.org/gsea/downloads.jsp}
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
#'
#' length(db)
#' ## [1] 18643
#'
#' head(names(db))
#' ## [1] "M3128"  "M11607" "M12599" "M5067"  "M10817" "M16694"
#'
#' str(db[[1]], nchar.max=20)
#' ## List of 8
#' ##  $ id          : chr "M3128"
#' ##  $ name        : chr "AAANWWTGC_UNKNOWN"
#' ##  $ category    : chr "C3"
#' ##  $ subcategory : chr "TFT"
#' ##  $ organism    : chr "Homo sapiens"
#' ##  $ desc        : chr "Genes having at lea"| __truncated__
#' ##  $ desc_full   : chr "Comprehensive ident"| __truncated__
#' ##  $ set         : chr [1:193] "MEF2C" "ATP1B1" "RORA" "CITED2" ...
#'
#' }
#'
#' ## A small sample of the MSigDB as imported by importMSigDB.xml is
#' ## included in gep2pep:
#'
#' db_sample <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#'
#' str(db_sample[[1]], nchar.max=20)
#' ## List of 8
#' ##  $ id          : chr "M3128"
#' ##  $ name        : chr "AAANWWTGC_UNKNOWN"
#' ##  $ db          : chr "C3"
#' ##  $ subcategory : chr "TFT"
#' ##  $ organism    : chr "Homo sapiens"
#' ##  $ desc        : chr "Genes having at lea"| __truncated__
#' ##  $ desc_full   : chr "Comprehensive ident"| __truncated__
#' ##  $ set         : chr [1:193] "MEF2C" "ATP1B1" "RORA" "CITED2" ...
#' @export
importMSigDB.xml <- function(fname) {

    xml <- xmlTreeParse(fname, useInternalNodes=TRUE)
    sets <- xml["/MSIGDB/GENESET"]

    ids <- sapply(sets, function(x) xmlAttrs(x)[["SYSTEMATIC_NAME"]])
    msigDB <- data.frame(
        id = ids,
        name = sapply(sets, function(x) xmlAttrs(x)[["STANDARD_NAME"]]),
        category = sapply(sets, function(x) xmlAttrs(x)[["CATEGORY_CODE"]]),
        subcategory = sapply(sets, function(x) {
            xmlAttrs(x)[["SUB_CATEGORY_CODE"]]}),
        organism = sapply(sets, function(x) xmlAttrs(x)[["ORGANISM"]]),
        desc = sapply(sets, function(x) xmlAttrs(x)[["DESCRIPTION_BRIEF"]]),
        desc_full = sapply(sets, function(x) xmlAttrs(x)[["DESCRIPTION_FULL"]]),
        set = sapply(sets, function(x) xmlAttrs(x)[["MEMBERS_SYMBOLIZED"]]),
        stringsAsFactors=FALSE
    )

    msigDB <- apply(msigDB, 1, as.list)
    msigDB <- lapply(msigDB, function(x) {
        x$set <- strsplit(x$set, ",")[[1]]; x})
    names(msigDB) <- ids
    
    return(msigDB)    
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
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for collection: C3_TFT
#' ## [15:45:06] Storing pathway data for collection: C3_MIR
#' ## [15:45:06] Storing pathway data for collection: C4_CGN
#'
#' getCollections(rp)
#' ## [1] "C3_TFT" "C3_MIR" "C4_CGN"
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
#' @param sets A database of pathways (or generic gene sets), see
#'     details.
#' @param name Name of the repository. Defaults to \code{NULL} (a
#'     generic name will be given).
#' @param description Description of the repository. If NULL
#'     (default), a generic description will be given.
#' @return An object of class \code{repo} that can be passed to
#'     \code{gep2pep} functions.
#' @details \code{sets} must be in the same format as output by
#'     \code{\link{importMSigDB.xml}}. It is a list where each item
#'     includes the following fields:
#'
#' \itemize{
#'
#'   \item{id: }{a unique identifier of the gene set}
#'
#'   \item{name: }{a descriptive name of the gene set}
#'
#'   \item{category: }{an ID for the category this gene set belongs to
#'   (for example "GO")}
#'
#'   \item{subcategory: }{an ID for the sub-category this gene set
#'   belongs to (for example "BP")}
#'
#'   \item{organism: }{name of the organism (for example
#'   "homo sapiens")}
#'
#'   \item{desc: }{text description of the gene set (typically one
#'   short sentence)}
#'
#'   \item{desc_full: }{a long, detailed description of the gene set}
#'
#'   \item{set: }{genes in the set, as a vector of (typically) gene
#'   symbols}
#' }
#' @seealso buildPEPs
#' @examples
#'
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for collection: C3_TFT
#' ## [15:45:06] Storing pathway data for collection: C3_MIR
#' ## [15:45:06] Storing pathway data for collection: C4_CGN
#'
#' rp
#' ##         ID    Dims     Size
#' ## C3_TFT_sets   10 18.16 kB
#' ## C3_MIR_sets   10 17.25 kB
#' ## C4_CGN_sets   10   6.9 kB
#'
#' unlink(repo_path, TRUE)
#' @export
createRepository <- function(path, sets, name=NULL, description=NULL)
{
    if(is.null(name))
        name <- "gep2pep repository"
    if(is.null(description))
        description <- paste("This repository contains pathway information",
                             "and possibly Pathway Expression Profiles", 
                             "created with the gep2pep package.")    
    
    if(file.exists(path)) {
        say("Can not create repository in existing folder", TRUE) 
    } else rp <- repo_open(path, TRUE)
    
    rp$project(name, description)

    ## un-hiding the project item
    ## silencing warning because of a bug in repo
    curwarn <- options()$warn
    options(warn=-1)
    rp$untag(name, "hide")
    options(warn=curwarn)
    
    db_ids <- makeCollectionIDs(sets)
    subdbs <- unique(db_ids)

    for(i in 1:length(subdbs))
    {
        dbi <- subdbs[i]
        say(paste("Storing pathway data for collection:", dbi))
        rp$put(sets[db_ids == dbi],
               paste0(subdbs[i], "_sets"),
               paste("Pathway information for collection", subdbs[i]),
               c("gep2pep", "sets"),
               prj = name)
    }    
    
    return(rp)
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
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
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
                  "a repository created with createRepository", TRUE))
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
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' ids <- makeCollectionIDs(db)
#'
#' unique(ids)
#' ## [1] "C3_TFT" "C3_MIR" "C4_CGN"
#'
#' db <- db[ids=="C3_MIR"]
#'
#' length(db)
#' ## [1] 10
#'
#' @export
makeCollectionIDs <- function(sets) {
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
#'     and each column to a perturbagen. Each column must include all
#'     ranks from 1 to the number of rows. Row and column names must
#'     be defined. Row names will be matched against gene identifiers
#'     in the pathways collections, and unrecognized gene names will
#'     not be used.
#' @param parallel If TRUE, gene sets will be processed in
#'     parallel. Requires a parallel backend.
#' @param replace_existing What to do if PEPs, identified by column
#'     names of \code{geps} are already present in the repository. If
#'     set to TRUE, they will be replaced, otherwise they will be
#'     skipped and only PEPs of new perturbagens will be computed and
#'     added.
#' @return Nothing. The computed PEPs will be available in the
#'     repository.
#' @seealso buildPEPs
#' @examples
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for collection: C3_TFT
#' ## [15:45:06] Storing pathway data for collection: C3_MIR
#' ## [15:45:06] Storing pathway data for collection: C4_CGN
#'
#' rp
#' ##          ID   Dims     Size
#' ## C3_TFT_sets   10   18.16 kB
#' ## C3_MIR_sets   10   17.25 kB
#' ## C4_CGN_sets   10     6.9 kB
#'
#' ## Loading sample gene expression profiles
#' geps <- readRDS(system.file("testgep.RDS", package="gep2pep"))
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
#' ##   C3_TFT_sets   10 18.16 kB
#' ##   C3_MIR_sets   10 17.25 kB
#' ##   C4_CGN_sets   10   6.9 kB
#' ##        C3_TFT    2  1.07 kB
#' ##        C3_MIR    2  1.07 kB
#' ##        C4_CGN    2  1.04 kB
#'
#' unlink(repo_path, TRUE)
#'
#' @export
buildPEPs <- function(rp, geps, parallel=FALSE, collections="all",
                      replace_existing=FALSE)
{
    checkGEPsFormat(geps)
    
    if(length(collections) == 1 && collections == "all") {
        dbs <- getCollections(rp)
        } else dbs <- collections

    for(i in 1:length(dbs))
    {
        say(paste0("Working on collection: ", dbs[i],
                   " (", i, "/", length(dbs), ")" ))

        if(rp$has(dbs[i]))
            curpeps <- colnames(rp$get(dbs[i])[[1]]) else curpeps <- NULL
        newpeps <- setdiff(colnames(geps), curpeps)
        oldpeps <- intersect(colnames(geps), curpeps)

        if(length(oldpeps > 0)) {
            msg <- paste0("Existing PEPs will be replaced: ",
                          paste(oldpeps, collapse=", "))
            if(!replace_existing)
                msg <- gsub("replaced", "skipped", msg)
            say(msg)
        }

        gepsi <- geps[, newpeps, drop=FALSE]
        
        if(length(newpeps) > 0) {
            thisdb <- rp$get(paste0(dbs[i], "_sets"))
            peps <- gep2pep(gepsi, thisdb, parallel)
            storePEPs(rp, dbs[i], peps)
        }

        cat("\n")        
    }
    
}


#' Export PertSEA or PathSEA results to XLS format
#'
#' The XLS output includes the full PertSEA or PathSEA results,
#' together with additional gene set information for the PertSEA.
#'
#' @inheritParams dummyFunction
#' @param results The output of \code{PertSEA} or \code{PathSEA}.
#' @param outname Name of the XLS file to be created.
#' @return Nothing.
#' @seealso PertSEA, PathSEA
#' @examples
#'
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- readRDS(system.file("testgep.RDS", package="gep2pep"))
#' buildPEPs(rp, geps)
#'
#' pgset <- c("(+)_chelidonine", "(+/_)_catechin")
#' psea <- PertSEA(rp, pgset)
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
    if(! tolower(type) %in% c("pertsea", "pathsea"))
        say("type must be on of: PertSEA, PathSEA", TRUE)

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

#' Performs Perturbagen Set Enrichment Analysis
#'
#' Perturbagen Set Enrichment Analysis (PertSEA) can be seen as a
#' Gene-SEA performed over rows (as opposed to columns) of a matrix of
#' GEPs. It tells how much a pathway is consistently dysregulated by a
#' set of perturbagens when compared to a statistical background of
#' perturbagens.
#'
#' @inheritParams dummyFunction
#' @param pgset A vector of names of perturbagens. Corresponding PEPs
#'     must exist in all the pathway collections currently in
#'     \code{rp}.
#' @param bgset The background against which to compare
#'     \code{pgset}. If set to \code{all} (default), all the remaining
#'     PEPs will be used. If provided, the corresponding PEPs must
#'     exist in all the pathway collections currently in \code{rp}.
#' @param details If TRUE (default) rank details will be reported for
#'     each perturbagen in \code{pgset}.
#' @return A list of 2, by names "PertSEA" and "details". The
#'     "PertSEA" entry is a 2-columns matrix including ESs and
#'     p-values (see details) for each pathway database and
#'     perturbagen. The "details" entry reports the rank of each
#'     perturbagen in \code{pgset} for each pathway.
#' @details For each pathway, all perturbagens are ranked by how much
#'     they dysregulate it (from the most UP-regulating to the most
#'     DOWN-regulating). Then, a Kolmogorov-Smirnov (KS) test is
#'     performed to compare the ranks assigned to perturbagens in
#'     \code{pgset} against the ranks assigned to perturbagens in
#'     \code{bgset}. A positive (negative) Enrichment Score (ES) of
#'     the KS test indicates whether each pathway is UP- (DOWN-)
#'     regulated by \code{pgset} as compared to \code{bgset}. A
#'     p-value is associated to the ES. See reference for further
#'     details.
#' @references Napolitano F. et al, Drug-set enrichment analysis: a
#'     novel tool to investigate drug mode of action. Bioinformatics
#'     32, 235-241 (2016).
#' @examples
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- readRDS(system.file("testgep.RDS", package="gep2pep"))
#' buildPEPs(rp, geps)
#'
#' pgset <- c("(+)_chelidonine", "(+/_)_catechin")
#' psea <- PertSEA(rp, pgset)
#' ## [16:26:58] Common perturbagens removed from bgset.
#' ## [16:26:58] Working on collection: C3_TFT
#' ## [16:26:58] Row-ranking collection...
#' ## [16:26:58] Computing enrichments...
#' ## [16:26:58] done.
#' ## [16:26:58] Working on collection: C3_MIR
#' ## [16:26:58] Row-ranking collection...
#' ## [16:26:58] Computing enrichments...
#' ## [16:26:58] done.
#' ## [16:26:58] Working on collection: C4_CGN
#' ## [16:26:58] Row-ranking collection...
#' ## [16:26:58] Computing enrichments...
#' ## [16:26:58] done.
#'
#' psea$PertSEA$C3_TFT
#' ##                ES  PV
#' ## M5067   1.0000000 0.2
#' ## M2501   1.0000000 0.2
#' ## M3128  -0.6666667 0.6
#' ## M11607  0.6666667 0.6
#' ## M16694 -0.6666667 0.6
#' ## M14463  0.6666667 0.6
#' ## M14686 -0.6666667 0.6
#' ## M12599  0.5000000 0.9
#' ## M10817 -0.5000000 0.9
#' ## M4831   0.3333333 1.0
#'
#' unlink(repo_path, TRUE)
#'
#' @export
PertSEA <- function(rp_peps, pgset, bgset="all", collections="all",
                    details=TRUE)
{
    dbs <- collections
    if(length(dbs) == 1 && dbs=="all") {
        dbs <- getCollections(rp_peps)
    } else {
        off <- setdiff(dbs, getCollections)
        if(length(off)>0)
            say(paste0("The following collections could not be found: ",
                       paste(off, collapse=", ")), TRUE)
    }

    if(details)
        thedetails <- list() else thedetails <- NULL
    
    res <- list()
    for(i in 1:length(dbs)) {        
        say(paste0("Working on collection: ", dbs[i]))

        peps <- rp_peps$get(dbs[i])

        if(length(bgset) == 1 && bgset=="all")
            bgset <- colnames(peps[[1]])
        
        if(length(intersect(pgset, bgset))>0) {
            bgset <- setdiff(bgset, pgset)
            say("Common perturbagens removed from bgset.")
        }
        
        rankingset <- c(bgset, pgset)

        if(!all(rankingset %in% colnames(peps$ES)))
            say(paste("The following perturbagens could not be found:",
                      paste(
                          setdiff(rankingset, colnames(peps)),
                          collapse = ", ")), TRUE)                      
        
        say(paste0("Row-ranking collection..."))
        ranked <- rankPEPsByRows(peps, rankingset)
        say(paste0("Computing enrichments..."))
        
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
        say("done.")
    }

    return(list(PertSEA=res, details=thedetails))
}


#' Performs Pathway Set Enrichment Analysis (PSEA)
#'
#' PathSEA is analogous to the Gene Set Enrichment Analysis (GSEA),
#' but for pathways instead of single genes. It can therefore be used
#' to look for perturbagens that consistently UP- or DOWN-regulate a
#' set of pathways.
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
#'     perturbagen in \code{pgset}.
#' @return A list of 2, by names "PathSEA" and "details". The
#'     "PathSEA" entry is a 2-columns matrix including ESs and
#'     p-values for each collection and perturbagen. The "details"
#'     entry reports the rank of each pathway in \code{pathways} for
#'     each perturbagen.
#' @details For each perturbagen, all pathways are ranked by how much
#'     they are dysregulated by it (from the most UP-regulated to the
#'     most DOWN-regulatied, according to the corresponding
#'     p-values). Then, a Kolmogorov-Smirnov (KS) test is performed to
#'     compare the ranks assigned to pathways in \code{pathways}
#'     against the ranks assigned to pathways in \code{bgsets}. A
#'     positive (negative) Enrichment Score (ES) of the KS test
#'     indicates whether each pathway is UP- (DOWN-) regulated by
#'     \code{pgset} as compared to \code{bgset}. A p-value is
#'     associated to the ES.
#' @examples
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- readRDS(system.file("testgep.RDS", package="gep2pep"))
#' buildPEPs(rp, geps)
#'
#' pathways <- c("M11607", "M10817", "M16694",         ## from C3_TFT
#'               "M19723", "M5038", "M13419", "M1094") ## from C4_CGN
#' subdb <- db[pathways]
#'
#' psea <- PathSEA(rp, subdb)
#' ## [15:35:29] Working on collection: C3_TFT
#' ## [15:35:29] Common pathway sets removed from bgset.
#' ## [15:35:29] Col-ranking collection...
#' ## [15:35:29] Computing enrichments...
#' ## [15:35:29] done.
#' ## [15:35:29] Working on collection: C4_CGN
#' ## [15:35:29] Common pathway sets removed from bgset.
#' ## [15:35:29] Col-ranking collection...
#' ## [15:35:29] Computing enrichments...
#' ## [15:35:29] done.
#'
#' psea$PathSEA$C3_TFT
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
    pathways <- pwList2pwStruct(pathways)
    if(bgsets != "all")
        bgsets <- pwList2pwStruct(bgsets)
       
    dbs <- collections ## just renaming
    if(length(dbs) == 1 && dbs=="all") {
        dbs <- getCollections(rp_peps)
    } else {
        off <- setdiff(dbs, getCollections)
        if(length(off)>0)
            say(paste0("The following collections could not be found: ",
                       paste(off, collapse=", ")), TRUE)
    }
    
    if(! all(names(pathways) %in% dbs))
        say(paste("Names of pathways should match",
                  "pathway collection names in the repository.", TRUE))

    if(details)
        thedetails <- list() else thedetails <- NULL     
    
    dbs <- names(pathways)
    res <- list()
    for(i in 1:length(pathways))
    {
        say(paste0("Working on collection: ", dbs[i]))
        gmd <- names(pathways[[dbs[i]]])

        if(length(bgsets) == 1 && bgsets=="all") {
            bgset <- names(rp_peps$get(paste0(dbs[i], "_sets")))
        } else bgset <- names(bgsets[[i]])

        if(length(intersect(gmd, bgset)) > 0) {
            bgset <- setdiff(bgset, gmd)
            say("Common pathway sets removed from bgset.")
        }
        rankingset <- c(gmd, bgset)

        peps <- rp_peps$get(dbs[i])
        notok <- rankingset[rankingset %in% rownames(peps)]
        if(length(notok)>0)
            say(paste0("Pathway set ids not found in ", dbs[i], ": ",
                       paste(notok, collapse=", ")), TRUE)

        say(paste0("Col-ranking collection..."))
        ranked <- rankPEPsByCols(peps, rankingset)
        say(paste0("Computing enrichments..."))

        ks <- apply(ranked, 2, function(col) ks.test.2(col[gmd], col[bgset]))

        PVs <- sapply(ks, get, x="p.value")
        sorter <- order(PVs)
        
        res[[dbs[i]]] <- data.frame(ES = sapply(ks, get, x="ES")[sorter],
                                    PV = PVs[sorter])

        if(details)
            thedetails[[dbs[i]]] <- ranked[gmd, sorter]
        
        say("done.")        
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
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
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
    for(i in 1:length(dbs)) {
        db <- rp$get(paste0(dbs[i], "_sets"))
        w <- sapply(db, function(x) gene %in% x$set)
        mods <- c(mods, db[w])
    }

    return(mods)        
}



gep2pep <- function(geps, sets, parallel=FALSE) {

    pathw <- sets
    genemat <- geps
    genes <- rownames(genemat)

    x <- list()
    sets <- sapply(unname(pathw), get, x="set")
    
    pb <- txtProgressBar()
    for(j in 1:ncol(genemat))
    {
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
    setTxtProgressBar(pb, 1)
    close(pb)

    ES <- matrix(NA, length(pathw), ncol(genemat))
    PV <- matrix(NA, length(pathw), ncol(genemat))
    
    for(i in 1:ncol(genemat)){
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


storePEPs <- function(rp, db_id, peps) {
    
    if(rp$has(db_id)) {
        curmat <- rp$get(db_id)

        ## checking what's new and what's old
        curpeps <- colnames(curmat$ES)
        newpeps <- setdiff(colnames(peps$ES), curpeps)
        oldpeps <- intersect(colnames(peps$ES), curpeps)

        ## replacing existing PEPs
        curmat[["ES"]][, oldpeps] <- peps$ES[, oldpeps]
        curmat[["PV"]][, oldpeps] <- peps$PV[, oldpeps]

        ## adding new PEPs
        peps$ES <- cbind(curmat$ES, peps$ES[, newpeps])
        peps$PV <- cbind(curmat$PV, peps$PV[, newpeps])        
    }

    
    say("Storing pathway expression profiles...")        
    rp$put(peps, db_id,
           paste0("Pathway data for collection ", db_id,
                  ". It contains 2 matrices: 1 for enrichement scores ",
                  "(signed Kolmogorov Smirnov statistic) and one for ",
                  "the corresponding p-values."),
           c("gep2pep", "pep"), replace=TRUE,
           depends = paste0(db_id, "_sets"),
           prj = get_repo_prjname(rp))

    say("Done.")
}

say <- function(txt, stopping=FALSE) {
    msg <- paste0("[",
                  format(Sys.time(), format="%H:%M:%S"),
                  "] ",
                  txt)
    
    if(stopping)
        stop(msg, call.=FALSE) else message(msg)
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
        rankingset <- 1:nrow(peps[["ES"]])
    
    PVs <- peps[["PV"]][rankingset, ]
    ESs <- peps[["ES"]][rankingset, ]
    x <- sapply(1:ncol(PVs), function(i) rankPEP(PVs[,i], ESs[,i]))
    colnames(x) <- colnames(PVs)
    rownames(x) <- rownames(PVs)
    return(x)
}


rankPEPsByRows <- function(peps, rankingset="all")
{
    if(length(rankingset) == 1 && rankingset == "all")
        rankingset <- 1:ncol(peps[["ES"]])

    ESs <- peps[["ES"]][, rankingset, drop=FALSE]
    x <- t(apply(-ESs, 1, rank, ties.method = "random", na.last="keep"))
    return(x)
}

attachInfo <- function(rp, results)
{
    type <- names(results)[[1]]
    dbs <- names(results[[type]])
    newres <- list()
    for(i in 1:length(results[[type]])) {
        db <- rp$get(paste0(dbs[i], "_sets"))
        resmat <- results[[type]][[i]]

        if(type == "PertSEA") {
            modIDs <- rownames(resmat)
            newres[[i]] <- cbind(
                Pathway = sapply(db[modIDs], get, x="name"),
                Description = sapply(db[modIDs], get, x="desc"),
                results[[type]][[i]],
                results[["details"]][[i]]
            )
        } else {
            res <- cbind(
                Perturbagen = rownames(results[[type]][[i]]),
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
        say("GEPs must be a two-dimensional matrix.", TRUE)
    if(dims[[1]]==1)
        say(paste("GEPs must have genes over rows,",
                  "while input matrix seems to have 1 row."),
            TRUE)

    gnames <- rownames(geps)
    if(is.null(gnames))
        say("GEPs must have rownames (as gene IDs).", TRUE)
    if(any(duplicated(gnames)))
        say("GEP rownames must be unique.", TRUE)
    
    pnames <- colnames(geps)
    if(is.null(pnames))
        say("GEPs must have colnames (as perturbagen IDs).", TRUE)
    if(any(duplicated(pnames)))
        say("GEP colnames must be unique.", TRUE)

    if(any(is.na(geps)))
        say("GEPs can not contain NAs.", TRUE)

    not_unique <- apply(geps, 2, function(x) {
        any(duplicated(x))})
    mins <- apply(geps, 2, min)
    maxs <- apply(geps, 2, max)    
    
    if(any(mins != 1) || any(maxs != dims[1]) || any(not_unique))
        say(paste("GEP columns must be ranks. Check",
                  "that each column is made of numbers from 1",
                  "to the number of rows."), TRUE)
        
}

get_repo_prjname <- function(rp) {
    prjname <- names(
        rp$entries()[sapply(lapply(rp$entries(), get, x="tags"),
                            `%in%`, x="#project")]
    )[[1]] ## take the first for robustness, should be only 1
}

## splits a flat list of pathways into sublists according to
## collections
pwList2pwStruct <- function(db) {
    collids <- makeCollectionIDs(db)
    colls <- unique(collids)
    out <- list()
    for(i in 1:length(colls))
        out[[colls[[i]]]] <- db[collids == colls[[i]]]
    return(out)
}
