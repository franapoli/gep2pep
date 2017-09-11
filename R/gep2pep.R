
#' gep2pep: creation and analysis of Pathway Expression Profiles
#'
#' @details
#' Pathway Expression Profiles (PEPs) are based on the expression of
#' pathways (defined as sets of genes) as opposed to individual
#' genes. \code{gep2pep} supports the convertion of gene expression
#' profiles (GEPs) to PEPs and performs enrichment analysis of
#' pathways or perturbagens. \code{gep2pep} creates a local repository
#' of gene sets, which can also be imported from the MSigDB
#' database. The local repository is in the \code{repo} format. When a
#' GEP is passed to \code{\link{buildPEPs}}, it refers to the stored
#' database of pathways to convert the GEP to a PEP and permanently
#' store the latter.
#'
#' One type of analysis that can be performed on PEPs and that is
#' directly supported by \code{gep2pep} is the Drug-Set Enrichment
#' Analysis (DSEA, see reference below). It finds pathways that
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
#'   \item{Pathway Expression Profile (PEP): }{a ranked list of
#'   pathways, as converted from a gene expression profile (GEP)
#'   according to a pathway collection.}
#'
#'   \item{perturbagen }{any condition inducing a GEP and therefore a
#'   PEP.}
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
#'     \code{\link{createRepository}}, and containing PEPs created with
#'     \code{\link{buildPEPs}}.
#' @param dbs A set of pathway collection names as returned by
#'     \code{getCollections}. If set to "all" (default), all the collections
#'     in \code{rp} will be used.
#' @param collections The set of collections to use. If set to "all",
#'     all collections in the repository will be used.
#' @return Nothing
dummyFunction <- function(rp, rp_peps, dbs, collections) {}


#' Imports pathways data from an MsigDB XML file.
#'
#' @param fname Path to an XML file downloaded from MsigDB.
#' @return A list of pathway entries (see
#'     \code{link{createRepository}}).
#'
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
#' ##  $ id       : chr "M3128"
#' ##  $ name     : chr "AAANWWTGC_UNKNOWN"
#' ##  $ db       : chr "C3"
#' ##  $ subdb    : chr "TFT"
#' ##  $ organism : chr "Homo sapiens"
#' ##  $ desc     : chr "Genes having at lea"| __truncated__
#' ##  $ desc_full: chr "Comprehensive ident"| __truncated__
#' ##  $ set      : chr [1:193] "MEF2C" "ATP1B1" "RORA" "CITED2" ...
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
#' ##  $ id       : chr "M3128"
#' ##  $ name     : chr "AAANWWTGC_UNKNOWN"
#' ##  $ db       : chr "C3"
#' ##  $ subdb    : chr "TFT"
#' ##  $ organism : chr "Homo sapiens"
#' ##  $ desc     : chr "Genes having at lea"| __truncated__
#' ##  $ desc_full: chr "Comprehensive ident"| __truncated__
#' ##  $ set      : chr [1:193] "MEF2C" "ATP1B1" "RORA" "CITED2" ...
#' @export
importMSigDB.xml <- function(fname) {

    xml <- xmlTreeParse(fname, useInternalNodes=TRUE)
    sets <- xml["/MSIGDB/GENESET"]

    ids <- sapply(sets, function(x) xmlAttrs(x)[["SYSTEMATIC_NAME"]])
    msigDB <- data.frame(
        id = ids,
        name = sapply(sets, function(x) xmlAttrs(x)[["STANDARD_NAME"]]),
        db = sapply(sets, function(x) xmlAttrs(x)[["CATEGORY_CODE"]]),
        subdb = sapply(sets, function(x) xmlAttrs(x)[["SUB_CATEGORY_CODE"]]),
        organism = sapply(sets, function(x) xmlAttrs(x)[["ORGANISM"]]),
        desc = sapply(sets, function(x) xmlAttrs(x)[["DESCRIPTION_BRIEF"]]),
        desc_full = sapply(sets, function(x) xmlAttrs(x)[["DESCRIPTION_FULL"]]),
        set = sapply(sets, function(x) xmlAttrs(x)[["MEMBERS_SYMBOLIZED"]]),
        stringsAsFactors=FALSE
    )


    msigDB <- apply(msigDB, 1, as.list)
    msigDB <- lapply(msigDB,
                     function(x) {x$set <- strsplit(x$set, ",")[[1]]; x})
    names(msigDB) <- ids
    
    return(msigDB)    
}



#' Returns the names of the pathway collections in a repository.
#' @inheritParams dummyFunction
#' @return Vector of collection names (see details).
#' @details Each collection in a database has a "db" name and a
#'     "subdb" name assigned, which are used to build the collection
#'     identifier as "db_subdb". This function obtains the identifiers
#'     by looking at data stored in the repository \code{rp} (entries
#'     that are tagged with "sets").
#' @examples
#'
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for DB: C3_TFT
#' ## [15:45:06] Storing pathway data for DB: C3_MIR
#' ## [15:45:06] Storing pathway data for DB: C4_CGN
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
#' @param path Path to a non-existing directory where the repository
#'     will be created.
#' @param sets A database of pathways, see details.
#' @param name Name of the repository. Defaults to \code{NULL} (a
#'     generic name will be given).
#' @param description Description of the repository. Defaults to
#'     \code{NULL} (a generic repository will be given).
#' @return An object of class \code{repo}.
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
#'   \item{db: }{an ID for the database this gene set belongs to (for
#'   example "GO")}
#'
#'   \item{subdb: }{an ID for the sub-database this gene set belongs
#'   to (for example "BP")}
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
#' ## [15:45:06] Storing pathway data for DB: C3_TFT
#' ## [15:45:06] Storing pathway data for DB: C3_MIR
#' ## [15:45:06] Storing pathway data for DB: C4_CGN
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
        name <- "gep2pep database"
    if(is.null(description))
        description <- paste("This database contains pathway information",
                             "and pathway expression profiles created by",
                             "the gep2pep package.")
    
    if(file.exists(path)) {
        say("Can not create repository in existing folder", TRUE) 
    } else rp <- repo_open(path, TRUE)
    
    rp$project(name, description)
    
    db_ids <- makeCollectionIDs(sets)
    subdbs <- unique(db_ids)

    for(i in 1:length(subdbs))
    {
        dbi <- subdbs[i]
        say(paste("Storing pathway data for DB:", dbi))
        rp$put(sets[db_ids == dbi],
               paste0(subdbs[i], "_sets"),
               paste("Pathway information for DB", subdbs[i]),
               c("gep2pep", "sets"))
    }    
    
    return(rp)
}



#' Creates a vector of collection labels for each pathway.
#' @param sets A pathway database in the same format as created by
#'     \code{importMSigDB.xml}.
#' @return A vector of identifiers, one per pathway, with the format:
#'     "db_subdb".
#' @details This function is useful to subset a database of pathway by
#'     collections.
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
    dbs <- sapply(sets, get, x="db")
    subdbs <- sapply(sets, get, x="subdb")
    subdbs[subdbs==""] <- dbs[subdbs==""]
    db_ids <- paste(dbs, subdbs, sep="_")
    return(db_ids)
}


#' Build PEPs from GEPs and stores them in the repository.
#' @inheritParams dummyFunction
#' @param geps A matrix of ranks where each row corresponds to a gene
#'     and each column to a perturbagen. Each column must include all
#'     ranks from 1 to the number of rows. Row and column names must
#'     be defined. Row names must match names used in pathways
#'     (unrecognized names will not be used).
#' @param parallel If TRUE, gene sets will be processed in
#'     parallel. Requires a parallel backend.
#' @param existing What to do if PEPs for a given DB of gene sets are
#'     already present. Can be one of:
#'
#' \itemize{
#'
#'   \item{stop: }{the default, throws an error. This is the safest
#'   approach.}
#'
#'   \item{skip: }{the existing DB will be skipped. This is especially
#'   useful to build a repository incrementally by repeatedly using
#'   the same call to \code{buildPEPs}}.
#'
#'   \item{overwrite: }{all the existing PEPs for the DB will be replaced by
#'   new PEPs.}
#'
#'   \item{append: }{the new PEPs will be appended to the existing
#'   ones. Useful to update the repository.}
#'
#' }
#'
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
#' ## [15:45:06] Storing pathway data for DB: C3_TFT
#' ## [15:45:06] Storing pathway data for DB: C3_MIR
#' ## [15:45:06] Storing pathway data for DB: C4_CGN
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
#' unlink(repo_path, TRUE)
#'
#' @export
buildPEPs <- function(rp, geps, parallel=FALSE, existing="stop")
{
    okargs <- c("stop", "overwrite", "skip", "append")
    
    if(length(existing)>1 || (! existing %in% okargs))
        say(paste0("existing must be one of: ",
                   paste(okargs, collapse=", ")), TRUE)
    
    dbs <- getCollections(rp)

    for(i in 1:length(dbs))
    {        
        say(paste0("Working on DB: ", dbs[i],
                   " (", i, "/", length(dbs), ")" ))

        if(rp$has(dbs[i]))
            if(existing != "stop") {
                say(paste0("Existing DB found, action: ", existing))
            } else {
                say(paste0("DB already present, check parameter: existing"),
                    TRUE)
            }

        if(rp$has(dbs[i]) && existing == "skip")
            next
        
        thisdb <- rp$get(paste0(dbs[i], "_sets"))
        peps <- gep2pep(geps, thisdb, parallel)

        storePEPs(rp, dbs[i], peps, existing)

        cat("\n")        
    }
    
}


#' Export PertSEA or PathSEA results to XLS format
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
    if(! tolower(type) %in% c("pertsea", "modsea"))
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

#' Performs Perturbagen Set Enrichment Analysis (PSEA)
#' @inheritParams dummyFunction
#' @param pgset A vector of names of perturbagens. Corresponding PEPs
#'     must exist in all the pathway collections currently in
#'     \code{rp}.
#' @param bgset The background against which to compare
#'     \code{pgset}. If set to \code{all} (default), all the remaining
#'     PEPs will be used. Corresponding PEPs must exist in all the
#'     pathway collections currently in \code{rp}.
#' @param details If TRUE (default) details will be reported for each
#'     perturbagen in \code{pgset}.
#' @return A list of 2, by names "PertSEA" and "details". The
#'     "PertSEA" entry is a 2-columns matrix including ESs and
#'     p-values (see details) for each pathway database and
#'     perturbagen. The "details" entry reports the rank of each
#'     perturbagen in \code{pgset} for each pathway.
#' @details Perturbagen Set Enrichment Analysis (PSEA) can be seen as
#'     a Gene-SEA performed over rows (as opposed to columns) of a
#'     profile matrix. For each pathway, all perturbagens are
#'     ranked and the p-value of a Kolmogorov-Smirnov test is computed
#'     for the ranks assigned to perturbagens in \code{pgset} against
#'     the ranks assigned to perturbagens in \code{bgset}. A positive
#'     (negative) Enrichment Score (ES) is also computed to indicate
#'     if each pathway is UP- (DOWN-) regulated by \code{pgset} as
#'     compared to \code{bgset}. See reference.
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
#' ## [16:26:58] Working on DB: C3_TFT
#' ## [16:26:58] Row-ranking DB...
#' ## [16:26:58] Computing enrichments...
#' ## [16:26:58] done.
#' ## [16:26:58] Working on DB: C3_MIR
#' ## [16:26:58] Row-ranking DB...
#' ## [16:26:58] Computing enrichments...
#' ## [16:26:58] done.
#' ## [16:26:58] Working on DB: C4_CGN
#' ## [16:26:58] Row-ranking DB...
#' ## [16:26:58] Computing enrichments...
#' ## [16:26:58] done.
#'
#' psea$PertSEA$C3_TFT
#' ##                ES  PV
#' ## M5067   1.0000000 0.2
#' ## M2501   1.0000000 0.2
#' ## M3128  -0.6666667 0.6
#' ## M11607  0.6666667 0.6
#' ## M10817 -0.6666667 0.6
#' ## M16694 -0.6666667 0.6
#' ## M14463  0.6666667 0.6
#' ## M14686 -0.6666667 0.6
#' ## M12599  0.5000000 0.9
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
            say(paste0("The following DBs could not be found: ",
                       paste(off, collapse=", ")), TRUE)
    }

    if(details)
        thedetails <- list() else thedetails <- NULL
    
    res <- list()
    for(i in 1:length(dbs)) {        
        say(paste0("Working on DB: ", dbs[i]))

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
        
        say(paste0("Row-ranking DB..."))
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
#' @inheritParams dummyFunction
#' @param modsets A list where each entry must be a vector of pathway
#'     IDs and must be named like a pathway database. PSEA will be
#'     performed for each database separately.
#' @param bgsets Another list like \code{modsets}, representing the
#'     statistical background for each database. If set to "all" (the
#'     default), all pathways not in \code{modsets} will be used.
#' @param details If TRUE (default) details will be reported for each
#'     perturbagen in \code{pgset}.
#' @return A list of 2, by names "PathSEA" and "details". The "PathSEA"
#'     entry is a 2-columns matrix including ESs and p-values for each
#'     collection and perturbagen. The "details" entry reports the
#'     rank of each pathway in \code{modsets} for each perturbagen.
#' @details PathSEA is the analogous of Gene Set Enrichment Analysis
#'     (GSEA), but for pathways instead of single genes.
#' @examples
#' db <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
#' repo_path <- file.path(tempdir(), "gep2pepTemp")
#'
#' rp <- createRepository(repo_path, db)
#' geps <- readRDS(system.file("testgep.RDS", package="gep2pep"))
#' buildPEPs(rp, geps)
#'
#' modsets <- list(C3_TFT = c("M11607", "M10817", "M16694"),
#'                 C4_CGN = c("M19723", "M5038", "M13419", "M1094"))
#'
#' psea <- PathSEA(rp, modsets)
#' ## [15:35:29] Working on DB: C3_TFT
#' ## [15:35:29] Common pathway sets removed from bgset.
#' ## [15:35:29] Col-ranking DB...
#' ## [15:35:29] Computing enrichments...
#' ## [15:35:29] done.
#' ## [15:35:29] Working on DB: C4_CGN
#' ## [15:35:29] Common pathway sets removed from bgset.
#' ## [15:35:29] Col-ranking DB...
#' ## [15:35:29] Computing enrichments...
#' ## [15:35:29] done.
#'
#' psea$PathSEA$C3_TFT
#' ##                         ES        PV
#' ## (_)_mk_801       0.7142857 0.1666667
#' ## (_)_atenolol     0.7142857 0.1666667
#' ## (+)_isoprenaline 0.5714286 0.4000000
#' ## (+/_)_catechin   0.5714286 0.4000000
#' ## (+)_chelidonine  0.3809524 0.8500000
#'
#' unlink(repo_path, TRUE)
#'
#' @export
PathSEA <- function(rp_peps, modsets, bgsets="all", collections="all",
                    details=TRUE)
{
    dbs <- collections
    if(length(dbs) == 1 && dbs=="all") {
        dbs <- getCollections(rp_peps)
    } else {
        off <- setdiff(dbs, getCollections)
        if(length(off)>0)
            say(paste0("The following DBs could not be found: ",
                       paste(off, collapse=", ")), TRUE)
    }
    
    if(! all(names(modsets) %in% dbs))
        say("Names of modsets should match pathway DB names.", TRUE)

    if(details)
        thedetails <- list() else thedetails <- NULL     
    
    dbs <- names(modsets)
    res <- list()
    for(i in 1:length(modsets))
    {
        say(paste0("Working on DB: ", dbs[i]))
        gmd <- modsets[[dbs[i]]]

        if(length(bgsets) == 1 && bgsets=="all") {
            bgset <- names(rp_peps$get(paste0(dbs[i], "_sets")))
        } else bgset <- bgsets[[i]]

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

        say(paste0("Col-ranking DB..."))
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


storePEPs <- function(rp, db_id, peps, existing) {

    okargs <- c("overwrite", "append", "stop")
    if(length(existing)>1 || (! existing %in% okargs))
        say(paste0("existing must be one of: ",
                   paste(okargs, collapse=", ")), TRUE)
    
    replace <- FALSE
    
    if(rp$has(db_id)) {

        if(existing == "stop")
            say("PEP already exists", TRUE)

        if(existing == "append") {
            say("Merging PEPs...")
            curpep <- rp$get(db_id)
            peps[["ES"]] <- cbind(curpep[["ES"]], peps[["ES"]])
            peps[["PV"]] <- cbind(curpep[["PV"]], peps[["PV"]])
            replace <- TRUE
        }
        
        if(existing == "overwrite")
            replace <- TRUE
    }
    
    say("Storing pathway expression profiles...")        
    rp$put(peps, db_id,
           paste0("Pathway data for DB ", db_id,
                  ". It contains 2 matrices: 1 for enrichement scores ",
                  "(signed Kolmogorov Smirnov statistic) and one for ",
                  "the corresponding p values."),
           c("gep2pep", "pep"), replace=replace)

    curids <- getCollections(rp)
    
    rp$put(colnames(peps[[1]]), "perturbagens",
           "Names of perturbagens inducing expression profiles",
           c("gep2pep", "meta"), replace=TRUE)

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

    ESs <- peps[["ES"]][, rankingset]
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


findGenePathways <- function(rp, gene)
{
    ## gene <- intersect(intersect(db[[3]]$set, db[[4]]$set), db[[7]]$set)[1]
    
    dbs <- getCollections(rp)
    mods <- list()
    for(i in 1:length(dbs)) {
        db <- rp$get(paste0(dbs[i], "_sets"))
        w <- sapply(db, function(x) gene %in% x$set)
        mods[[dbs[i]]] <- db[w]
    }

    return(mods)        
}
