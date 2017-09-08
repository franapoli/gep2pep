#' gep2mep: creation and analysis of Module Expression Profiles
#'
#' @details
#' Gene-module expression profiles (MEPs) are based on the expression
#' of gene modules (sets of genes) as opposed to individual
#' genes. This package converts gene expression profiles to MEPs and
#' performs enrichment analysis of modules or perturbagens.
#' \code{gep2mep} creates a local repository of gene modules, which
#' can also be imported from the MSigDB database. The local repository
#' is in the \code{repo} format. When a gene expression profile (GEP)
#' is passed to \code{\link{buildMEPs}}, it refers to the existing
#' gene sets to convert the GEP to a MEP.
#'
#' One type of analysis that can be performed on MEPs and that is
#' directly supported by \code{gep2mep} is the Drug-Set Enrichment
#' Analysis (DSEA, see reference below). It finds gene modules that
#' are consistently dysregulated by a set of drugs, as opposed to a
#' background of other drugs. Of course MEPs may refer to
#' non-pharmacological perturbagens (genetic perturbations, disease
#' states, etc.) for analogous analyses. See \code{\link{PertSEA}}
#' function.
#'
#' A complementary approach is that of finding perturbagens that
#' consistently dysregulate a set of gene modules. This is the
#' module-based version of the Gene Set Enrichment Analysis
#' (GSEA). See \code{\link{ModSEA}}.
#'
#' @references Napolitano F. et al, Drug-set enrichment analysis: a
#'     novel tool to investigate drug mode of action. Bioinformatics
#'     32, 235–241 (2016).
#' 
#' @docType package
#' @name gep2mep-package
#' @author Francesco Napolitano \email{franapoli@@gmail.com}
#' @aliases gep2mep
NULL

## repo is for storage of repositories
#'@import repo
## XML is to import MSigDB data
#'@import XML
## foreach is for easy parallelization support
#'@import foreach
## utils is for txtProgressBar
#'@import utils
## stats is for ks.test
#'@import stats




#' Dummy function for parameter inheritance
#' @param rp A repository created by \code{buildEmptyDB}.
#' @param rp_meps A repository created with \code{buildEmptyDB}, and
#'     containing MEPs created with \code{buildMEPs}.
#' @return Nothing
dummyFunction <- function(rp) {}


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


#' Imports gene set information from an MsigDB XML file.
#'
#' @param fname Path to an XML file downloaded from MsigDB
#' @return A list of gene set entries (see details)
#' @details The format required by gep2mep for a database of gene sets
#'     entries, provided by this function, is a list where each item
#'     includes the following fields:
#' \itemize{
#' \item{setid: }{a unique identifier of the gene set}
#'
#' \item{setname: }{a descriptive name of the gene set}
#'
#' \item{db: }{an ID for the database this gene set belongs to (for
#' example "GO")}
#'
#' \item{subdb: }{an ID for the sub-database this gene set belongs to
#' (for example "BP")}
#'
#' \item{organism: }{name of the organism (for example "homo sapiens")}
#'
#' \item{desc: }{text description of the gene set (typically one short
#' sentence)}
#'
#' \item{desc_full: }{a long, detailed description of the gene set}
#'
#' \item{set: }{genes in the set, as a vector of (typically) gene
#' symbols}
#' }
#'
#' @examples
#' \dontrun{
#' 
#' ## To run this example, first obtain the MSigDB database in XML
#' ## format (see
#' ## http://software.broadinstitute.org/gsea/downloads.jsp). It
#' ## assumed that this file is available as "msigdb_v6.0.xml".
#' 
#' db <- importMsigDB.xml("msigdb_v6.0.xml")
#'
#' ## The db is now in an acceptable format to create a local
#' ## repository using buildEmptyDB
#'
#' length(db)
#' ## [1] 18643
#'
#' head(names(db))
#' ## [1] "M3128"  "M11607" "M12599" "M5067"  "M10817" "M16694"
#'
#' str(db[[1]], nchar.max=20)
#' ## List of 8
#' ##  $ setid    : chr "M3128"
#' ##  $ setname  : chr "AAANWWTGC_UNKNOWN"
#' ##  $ db       : chr "C3"
#' ##  $ subdb    : chr "TFT"
#' ##  $ organism : chr "Homo sapiens"
#' ##  $ desc     : chr "Genes having at lea"| __truncated__
#' ##  $ desc_full: chr "Comprehensive ident"| __truncated__
#' ##  $ set      : chr [1:193] "MEF2C" "ATP1B1" "RORA" "CITED2" ...
#'
#' }
#'
#' @export
importMsigDB.xml <- function(fname) {

    xml <- xmlTreeParse(fname, useInternalNodes=TRUE)
    sets <- xml["/MSIGDB/GENESET"]

    ids <- sapply(sets, function(x) xmlAttrs(x)[["SYSTEMATIC_NAME"]])
    msigDB <- data.frame(
        setid = ids,
        setname = sapply(sets, function(x) xmlAttrs(x)[["STANDARD_NAME"]]),
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

## #' Converts Gene Expression Profiles (GEPs) to Module-EPs (MEPs).
## #' @param geps A matrix of ranks where each row corresponds to a gene
## #'     and each column to a perturbagen. Each column must include all
## #'     ranks from 1 to the number of rows. Row and column names must
## #'     be defined.
## #' @param gmd A database of gene modules. See \code{importMsigDB.xml}
## #'     for the format.
## #' @param parallel If TRUE, gene sets will be processed in
## #'     parallel. Requires a parallel backend.
## #' @return A list of two matrices, one for Enrichment Scores and one
## #'     for p-values. Each entry (i,j) refers to gene set i and
## #'     perturbagen j.
## #' @seealso buildMEPs
## #' @export
gep2mep <- function(geps, gmd, parallel=FALSE) {

    pathw <- gmd
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
        sapply(pathw, "get", x="setid")
    colnames(ES) <- colnames(PV) <- colnames(genemat)
    
    return(list(ES=ES, PV=PV))
}


#' Returns the names of the module databases in a repository.
#' @inheritParams dummyFunction
#' @return Vector of names.
#' @details The names are obtained from the repository entries that
#'     are tagged with "gmd", and removing the "_gmd" suffix.
#' @examples
#'
#' db <- readRDS(system.file("testgmd.RDS", package="gep2mep"))
#' repo_path <- file.path(tempdir(), "gep2mepTemp")
#'
#' rp <- buildEmptyDB(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for DB: C3_TFT
#' ## [15:45:06] Storing pathway data for DB: C3_MIR
#' ## [15:45:06] Storing pathway data for DB: C4_CGN
#'
#' getDBlist(rp)
#' ## [1] "C3_TFT" "C3_MIR" "C4_CGN"
#'
#' unlink(repo_path, TRUE)
#'
#' @export
getDBlist <- function(rp)
{    
    w <- sapply(lapply(rp$entries(), get, x="tags"), `%in%`, x="gmd")
    res <- gsub("_gmd", "", names(w[w]))
    return(res)
}


#' Creates an empty DB of gene sets for converting GEPs.
#' @param path Path to an empty folder where the repository will be
#'     created.
#' @param gmd A list of gene modules (see \code{importMsigDB.xml}).
#' @param name Name of the repository. Defaults to \code{NULL} (a
#'     generic name will be given).
#' @param description Description of the repository. Defaults to \code{NULL}
#'     (a generic repository will be given).
#' @return An object of class \code{repo}.
#' @seealso buildMEPs, buildDBIDs
#' @examples
#'
#' db <- readRDS(system.file("testgmd.RDS", package="gep2mep"))
#' repo_path <- file.path(tempdir(), "gep2mepTemp")
#'
#' rp <- buildEmptyDB(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for DB: C3_TFT
#' ## [15:45:06] Storing pathway data for DB: C3_MIR
#' ## [15:45:06] Storing pathway data for DB: C4_CGN
#'
#' rp
#' ##         ID Dims     Size
#' ## C3_TFT_gmd   10 18.16 kB
#' ## C3_MIR_gmd   10 17.25 kB
#' ## C4_CGN_gmd   10   6.9 kB
#'
#' unlink(repo_path, TRUE)
#' @export
buildEmptyDB <- function(path, gmd, name=NULL, description=NULL)
{
    if(is.null(name))
        name <- "gep2mep database"
    if(is.null(description))
        description <- paste("This database contains pathway information",
                             "and pathway expression profiles created by",
                             "the gep2mep package.")
    
    if(file.exists(path)) {
        say("Can not create DB in existing folder", TRUE) 
    } else rp <- repo_open(path, TRUE)
    
    rp$project(name, description)
    
    db_ids <- buildDBids(gmd)
    subdbs <- unique(db_ids)

    for(i in 1:length(subdbs))
    {
        dbi <- subdbs[i]
        say(paste("Storing pathway data for DB:", dbi))
        rp$put(gmd[db_ids == dbi],
               paste0(subdbs[i], "_gmd"),
               paste("Pathway information for DB", subdbs[i]),
               c("gep2mep", "gmd"))
    }    
    
    ## rp$put(unique(db_ids), "DB list",
    ##        "IDs of pathway sub-databases included in this repository",
    ##        c("gep2mep", "meta"))
    
    return(rp)
}


storeMEPs <- function(rp, db_id, peps, existing) {

    okargs <- c("overwrite", "append", "stop")
    if(length(existing)>1 || (! existing %in% okargs))
        say(paste0("existing must be one of: ",
                   paste(okargs, collapse=", ")), TRUE)
    
    replace <- FALSE
    
    if(rp$has(db_id)) {

        if(existing == "stop")
            say("MEP already exists", TRUE)

        if(existing == "append") {
            say("Merging MEPs...")
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
           c("gep2mep", "mep"), replace=replace)

    curids <- getDBlist(rp)
    
    rp$put(colnames(peps[[1]]), "perturbagens",
           "Names of perturbagens inducing expression profiles",
           c("gep2mep", "meta"), replace=TRUE)

    say("Done.")
}



#' Creates unique identifiers for gene module databases.
#' @param gmds A gene module database in the same format created by
#'     \code{importMsigDB.xml}.
#' @return A vector of identifiers, one per gene module, with the
#'     format: "db_subdb".
#' @seealso importMsigDB.xml
#' @examples
#' db <- readRDS(system.file("testgmd.RDS", package="gep2mep"))
#' ids <- buildDBids(db)
#'
#' length(ids)
#' ## [1] 30
#'
#' unique(ids)
#' ## [1] "C3_TFT" "C3_MIR" "C4_CGN"
#'
#' @export
buildDBids <- function(gmds) {
    dbs <- sapply(gmds, get, x="db")
    subdbs <- sapply(gmds, get, x="subdb")
    subdbs[subdbs==""] <- dbs[subdbs==""]
    db_ids <- paste(dbs, subdbs, sep="_")
    return(db_ids)
}

#' Build MEPs using an existing repository and stores them in it.
#' @inheritParams dummyFunction
#' @param geps A matrix of ranks where each row corresponds to a gene
#'     and each column to a perturbagen. Each column must include all
#'     ranks from 1 to the number of rows. Row and column names must
#'     be defined.
#' @param parallel If TRUE, gene sets will be processed in
#'     parallel. Requires a parallel backend.
#' @param existing What to do if MEPs for a given DB of gene sets are
#'     already present. Can be one of:
#'
#' + stop: the default, throws an error. This is the safest approach.
#'
#' + skip: the existing DB will be skipped. This is especially useful
#' to build a repository incrementally by repeatedly using the same
#' call to \code{buildMEPs}.
#'
#' + overwrite: all the existing MEPs for the DB will be replaced by
#' new MEPs.
#'
#' + append: the new MEPs will be appended to the existing
#' ones. Useful to update the repository.
#'
#' @return Nothing. The computed MEPs will be available in the
#'     repository.
#' @seealso buildMEPs
#' @examples
#' db <- readRDS(system.file("testgmd.RDS", package="gep2mep"))
#' repo_path <- file.path(tempdir(), "gep2mepTemp")
#'
#' rp <- buildEmptyDB(repo_path, db)
#' ## Repo root created.
#' ## Repo created.
#' ## [15:45:06] Storing pathway data for DB: C3_TFT
#' ## [15:45:06] Storing pathway data for DB: C3_MIR
#' ## [15:45:06] Storing pathway data for DB: C4_CGN
#'
#' rp
#' ##         ID Dims     Size
#' ## C3_TFT_gmd   10 18.16 kB
#' ## C3_MIR_gmd   10 17.25 kB
#' ## C4_CGN_gmd   10   6.9 kB
#'
#' ## Loading sample gene expression profiles
#' geps <- readRDS(system.file("testgep.RDS", package="gep2mep"))
#'
#' geps[1:3,1:3]
#' ##       (+)_chelidonine (+)_isoprenaline (+/_)_catechin
#' ## AKT3             2050             2703          10435
#' ## MED6             8495            10020            984
#' ## NR2E3            9021             2769          11127
#'
#' buildMEPs(rp, geps)
#'
#' rp
#' ##           ID Dims     Size
#' ##   C3_TFT_gmd   10 18.16 kB
#' ##   C3_MIR_gmd   10 17.25 kB
#' ##   C4_CGN_gmd   10   6.9 kB
#' ##       C3_TFT    2  1.07 kB
#' ## perturbagens    5    114 B
#' ##       C3_MIR    2  1.07 kB
#' ##       C4_CGN    2  1.04 kB
#' unlink(repo_path, TRUE)
#'
#' @export
buildMEPs <- function(rp, geps, parallel=FALSE, existing="stop")
{
    okargs <- c("stop", "overwrite", "skip", "append")
    
    if(length(existing)>1 || (! existing %in% okargs))
        say(paste0("existing must be one of: ",
                   paste(okargs, collapse=", ")), TRUE)
    
    dbs <- getDBlist(rp)

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
        
        thisdb <- rp$get(paste0(dbs[i], "_gmd"))
        peps <- gep2mep(geps, thisdb, parallel)

        storeMEPs(rp, dbs[i], peps, existing)

        cat("\n")        
    }
    
}

## #' Shows repository statistics
## #' @param rp A gep2mep repository repository.
## #' @return Nothing.
## #' @examples
## #'
## #' db <- readRDS(system.file("testgmd.RDS", package="gep2mep"))
## #' repo_path <- file.path(tempdir(), "gep2mepTemp")
## #'
## #' rp <- buildEmptyDB(repo_path, db)
## #'
## #' dbStats(rp)
## #' ## Number of perturbagens: 5
## #' ## Databases and their sizes:
## #' ##    C3_TFT:    10
## #' ##    C3_MIR:    10
## #' ##    C4_CGN:    10
## #'
## #' unlink(repo_path, TRUE)
## #'
## #' @export
## dbStats <- function(rp) {
##     pert <- length(rp$get("perturbagens"))
##     dbs <- getDBlist(rp)

##     message(paste("Number of perturbagens:", pert))
##     message("Databases and their sizes:")
##     for(i in 1:length(dbs))
##         message(paste0("   ", dbs[i], ":\t",
##                    length(rp$get(paste0(dbs[i], "_gmd")))))
## }


## newMEPs <- function(rp, geps, gmd, id, parallel=F, overwrite=F)
## {
##     n <- length(subdbs)
##     peps <- gep2pep(geps, thisdb, parallel)
##     addDBtoRepo(rp, id, peps, gmd)
## }



## mergeRepoDBs <- function(path1, path2, path3, overwrite=F)
## {
##     quitDiffDBs <- function()
##         say("The 2 DBs must have the same sub-DBs and pathways", T)

##     if(!file.exists(path1))
##         say(paste("Path does not exist:", path1), T)
##     rp1 <- repo_open(path1)
##     if(!file.exists(path2))
##         say(paste("Path does not exist:", path2), T)
##     rp2 <- repo_open(path2)
##     rp3 <- buildEmptyRepoDB(path3, overwrite)

##     dbs1 <- rp1$get("DB list")
##     dbs2 <- rp2$get("DB list")
##     if(!setequal(dbs1, dbs2))
##         quitDiffDBs()

##     pert1 <- rp1$get("perturbagens")
##     pert2 <- rp2$get("perturbagens")

##     say(paste("Found", length(pert1), "perturbagens in", path1))
##     say(paste("Found", length(pert2), "perturbagens in", path2))
##     say(paste("Found", length(intersect(pert1, pert2)),
##                   "perturbagens in common, which will be",
##                   "treated as different."))

##     for(i in 1:length(dbs))
##     {
##         say(paste("Merging sub-DB:", dbs[i]))
##         prl1 <- rp1$get(dbs[i])
##         pw1 <- rp1$get(paste0(dbs[i], "_gmd"))
##         prl2 <- rp2$get(dbs[i])
##         pw2 <- rp2$get(paste0(dbs[i], "_gmd"))
##         if(!identical(pw1, pw2))
##             quitDiffDBs()
##         prl3 <- cbind(prl1, prl2)
##         addDBtoRepo(rp, dbs[i], cbind(prl1, prl2), pw1)        
##     }

## }


rankMEPsByCols <- function(peps, rankingset="all")
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


rankMEPsByRows <- function(peps, rankingset="all")
{
    if(length(rankingset) == 1 && rankingset == "all")
        rankingset <- 1:ncol(peps[["ES"]])

    ESs <- peps[["ES"]][, rankingset]
    x <- t(apply(-ESs, 1, rank, ties.method = "random", na.last="keep"))
    return(x)
}


#' Export PertSEA or ModSEA results to XLS format
#' @inheritParams dummyFunction
#' @param results The output of \code{PertSEA} or \code{ModSEA}.
#' @param outname Name of the XLS file to be created.
#' @return Nothing.
#' @seealso PertSEA, ModSEA
#' @examples
#'
#' db <- readRDS(system.file("testgmd.RDS", package="gep2mep"))
#' repo_path <- file.path(tempdir(), "gep2mepTemp")
#'
#' rp <- buildEmptyDB(repo_path, db)
#' geps <- readRDS(system.file("testgep.RDS", package="gep2mep"))
#' buildMEPs(rp, geps)
#'
#' pgset <- c("(+)_chelidonine", "(+/_)_catechin")
#' psea <- PertSEA(rp, pgset)
#'
#' exportSEA(rp, psea)
#'
#' unlink(repo_path, TRUE)
#'
#' @export
exportSEA <- function(rp, results, outname=NULL)
{
    type <- names(results)[[1]]
    if(! tolower(type) %in% c("pertsea", "modsea"))
        say("type must be on of: PertSEA, ModSEA", T)

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

attachInfo <- function(rp, results)
{
    type <- names(results)[[1]]
    dbs <- names(results[[type]])
    newres <- list()
    for(i in 1:length(results[[type]])) {
        db <- rp$get(paste0(dbs[i], "_gmd"))
        resmat <- results[[type]][[i]]

        if(type == "PertSEA") {
            modIDs <- rownames(resmat)
            newres[[i]] <- cbind(
                Module = sapply(db[modIDs], get, x="setname"),
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
                modnames <- sapply(db[modIDs], get, x="setname")
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

#' Performs Perturbagen Set Enrichment Analysis
#' @inheritParams dummyFunction
#' @param pgset A vector of names of perturbagens. Corresponding MEPs
#'     must exist in all the gene module databases currently in
#'     \code{rp}.
#' @param bgset The background against which to compare
#'     \code{pgset}. If set to \code{all} (default), all the remaining
#'     MEPs will be used.
#' @param dbs A subset of database names as returned by
#'     \code{getDBlist}. If set to "all" (default), all of them will
#'     be used.
#' @param details If TRUE (default) details will be reported for each
#'     perturbagen in \code{pgset}.
#' @return A list of 2, by names "PertSEA" and "details". The
#'     "PertSEA" entry is a 2-columns matrix including ESs and
#'     p-values (see details) for each gene module database and
#'     perturbagen. The "details" entry reports the rank of each
#'     perturbagen in \code{pgset} for each gene module.
#' @details Perturbagen Set Enrichment Analysis (PSEA) can be seen as
#'     a Gene-SEA performed over rows (as opposed to columns) of a
#'     profile matrix. For each gene module, all perturbagens are
#'     ranked and the p-value of a Kolmogorov-Smirnov test is computed
#'     for the ranks assigned to perturbagens in \code{pgset} against
#'     the ranks assigned to perturbagens in \code{bgset}. A positive
#'     (negative) Enrichment Score (ES) is also computed to indicate
#'     if each gene module is UP- (DOWN-) regulated by \code{pgset} as
#'     compared to \code{bgset}.
#' @examples
#' db <- readRDS(system.file("testgmd.RDS", package="gep2mep"))
#' repo_path <- file.path(tempdir(), "gep2mepTemp")
#'
#' rp <- buildEmptyDB(repo_path, db)
#' geps <- readRDS(system.file("testgep.RDS", package="gep2mep"))
#' buildMEPs(rp, geps)
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
PertSEA <- function(rp_meps, pgset, bgset="all", dbs="all", details=TRUE)
{
    if(length(dbs) == 1 && dbs=="all") {
        dbs <- getDBlist(rp_meps)
    } else {
        off <- setdiff(dbs, getDBlist)
        if(length(off)>0)
            say(paste0("The following DBs could not be found: ",
                       paste(off, collapse=", ")), TRUE)
    }

    if(length(bgset) == 1 && bgset=="all")
        bgset <- rp_meps$get("perturbagens")
    if(length(intersect(pgset, bgset))>0) {
        bgset <- setdiff(bgset, pgset)
        say("Common perturbagens removed from bgset.")
    }
    rankingset <- c(bgset, pgset)

    if(details)
        thedetails <- list() else thedetails <- NULL
    
    res <- list()
    for(i in 1:length(dbs)) {
        say(paste0("Working on DB: ", dbs[i]))
        say(paste0("Row-ranking DB..."))
        ranked <- rankMEPsByRows(rp_meps$get(dbs[i]), rankingset)
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


findGeneModules <- function(rp, gene)
{
    ## testing
    ## gene <- intersect(intersect(db[[3]]$set, db[[4]]$set), db[[7]]$set)[1]
    
    dbs <- getDBlist(rp)
    mods <- list()
    for(i in 1:length(dbs)) {
        db <- rp$get(paste0(dbs[i], "_gmd"))
        w <- sapply(db, function(x) gene %in% x$set)
        mods[[dbs[i]]] <- db[w]
    }

    return(mods)        
}



#' Performs Module Set Enrichment Analysis
#'
#' @inheritParams dummyFunction
#' @param modsets A list. Each entry must be a vector of module IDs
#'     and must be named like a module database. MSEA will be
#'     performed for each database separately.
#' @param bgsets A list with the same format as \code{modsets},
#'     representing the statistical background for each database. If
#'     set to "all" (the default), all modules not in \code{modsets}
#'     will be used.
#' @return A list of 2, by names "ModSEA" and "details". The "ModSEA"
#'     entry is a 2-columns matrix including ESs and p-values for each
#'     collection and perturbagen. The "details" entry reports the
#'     rank of each module in \code{modsets} for each perturbagen.
#' @details ModSEA is the analogous of Gene Set Enrichment Analysis
#'     (GSEA), but for gene modules instead of single genes.
#' @examples
#' db <- readRDS(system.file("testgmd.RDS", package="gep2mep"))
#' repo_path <- file.path(tempdir(), "gep2mepTemp")
#'
#' rp <- buildEmptyDB(repo_path, db)
#' geps <- readRDS(system.file("testgep.RDS", package="gep2mep"))
#' buildMEPs(rp, geps)
#'
#' modsets <- list(C3_TFT = c("M11607", "M10817", "M16694"),
#'                 C4_CGN = c("M19723", "M5038", "M13419", "M1094"))
#'
#' msea <- ModSEA(rp, modsets)
#' ## [15:35:29] Working on DB: C3_TFT
#' ## [15:35:29] Common module sets removed from bgset.
#' ## [15:35:29] Col-ranking DB...
#' ## [15:35:29] Computing enrichments...
#' ## [15:35:29] done.
#' ## [15:35:29] Working on DB: C4_CGN
#' ## [15:35:29] Common module sets removed from bgset.
#' ## [15:35:29] Col-ranking DB...
#' ## [15:35:29] Computing enrichments...
#' ## [15:35:29] done.
#'
#' msea$ModSEA$C3_TFT
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
ModSEA <- function(rp_meps, modsets, bgsets="all", details=T)
{
    
    dbs <- getDBlist(rp_meps)
    if(! all(names(modsets) %in% dbs))
        say("Names of modsets should match gene module DB names.", T)

    if(details)
        thedetails <- list() else thedetails <- NULL     
    
    dbs <- names(modsets)
    res <- list()
    for(i in 1:length(modsets))
    {
        say(paste0("Working on DB: ", dbs[i]))
        gmd <- modsets[[dbs[i]]]

        if(length(bgsets) == 1 && bgsets=="all") {
            bgset <- names(rp_meps$get(paste0(dbs[i], "_gmd")))
        } else bgset <- bgsets[[i]]

        if(length(intersect(gmd, bgset)) > 0) {
            bgset <- setdiff(bgset, gmd)
            say("Common module sets removed from bgset.")
        }
        rankingset <- c(gmd, bgset)

        peps <- rp_meps$get(dbs[i])
        notok <- rankingset[rankingset %in% rownames(peps)]
        if(length(notok)>0)
            say(paste0("Module set ids not found in ", dbs[i], ": ",
                       paste(notok, collapse=", ")), T)

        say(paste0("Col-ranking DB..."))
        ranked <- rankMEPsByCols(peps, rankingset)
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

    return(list(ModSEA=res, details=thedetails))
}


gsea <- function(S, ranks_list, check=FALSE, alternative = "two.sided")
                                        #                ,leadedge = FALSE)
{
    S <- S[!(is.na(S))]
    S1 <- ranks_list[S]
    S2 <- ranks_list[-S]

    if(length(S1)<1 || length(S2)<1 || all(is.na(S1)) || all(is.na(S2)))
        return(list(ES=NA, p=NA))

    ks <- ks.test.2(S1, S2, alternative=alternative, maxCombSize=10^10)

    ## lead=NA
    ## if(leadedge){
    ##     sSm <- sort(Sm)
    ##     if(ks$ES > 0)
    ##         lead=sSm[sSm<=ks$edge] else lead=sSm[sSm>=ks$edge]
    ## }

    ## return(list(ES=ks$ES, p=ks$"p.value", edge=ks$edge, lead=lead));
    return(list(ES=ks$ES, p=ks$"p.value", edge=ks$edge));
}


