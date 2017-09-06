#' gep2mep: working with Module Expression Profiles
#'
#' @details TODO
#' 
#' @docType package
#' @name gep2mep-package
#' @author Francesco Napolitano \email{franapoli@@gmail.com}
#' @aliases gep2mep-package
NULL


#'@import repo
## XML is to import MSigDB data
#'@import XML 
#'@import foreach
## utils is for txtProgressBar
#'@import utils
#'@import stats


ks.sign <- function (x, y)
{
    n.x <- as.double(length(x))
    n.y <- length(y)

    w <- c(x, y)        
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))

    return(sign(z[which.max(abs(z))]))
}



## S1 <- sample(1:100, 10);
## S2 <- (1:100)[-S1];
## system.time(for(i in 1:1000) {
##                 x <- sign(ks.test.2(S1, S2)$ES);
## }
## )
## system.time(for(i in 1:1000) {
##                 x <- sign(ks.sign(S1, S2));
## }
## )

## for(i in 1:10000) {
## S1 <- sample(1:100, 10);
## S2 <- (1:100)[-S1];
##     if(sign(ks.sign(S1, S2))!=sign(ks.test.2(S1, S2)$ES))
##         stop();
## }

ks.test.2 <- function(x, y, ...) {
    ks <- ks.test(x, y, ...)
    ks[["ES"]] <- unname(ks.sign(x, y)*ks$statistic)
    return(ks)
}


repoExists <- function(path)
{
    return(file.exists(file.path(path, "R_repo.RDS")))
}

say <- function(txt, stopping=F) {
    msg <- paste0("[",
                  format(Sys.time(), format="%H:%M:%S"),
                  "] ",
                  txt)
    
    if(stopping)
        stop(msg, call.=F) else message(msg)
}



#' Imports gene set information from an MsigDB XML file.
#'
#' @param fname Path to an XML file downloaded from MsigDB
#' @return A list of gene set entries (see details)
#' @details The format required by gep2mep for a database of gene sets
#'     entries, provided by this function, is a list where each item
#'     includes the following fields:
#'
#' + setid: a unique identifier of the gene set
#'
#' + setname: a descriptive name of the gene set
#'
#' + db: an ID for the database this gene set belongs to (for
#' example "GO")
#'
#' + subdb: an ID for the sub-database this gene set belongs to (for
#' example "BP")
#'
#' + organism: name of the organism (for example "homo sapiens")
#'
#' + desc: text description of the gene set (typically one short
#' sentence)
#'
#' + desc_full: a long, detailed description of the gene set
#'
#' + set: genes in the set, as a vector of (typically) gene symbols
#' 
#' @export
importMsigDB.xml <- function(fname) {

    xml <- xmlTreeParse(fname, useInternalNodes=T)
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
        stringsAsFactors=F
    )


    msigDB <- apply(msigDB, 1, as.list)
    msigDB <- lapply(msigDB, function(x) {x$set <- strsplit(x$set, ",")[[1]]; x})
    names(msigDB) <- ids
    
    return(msigDB)    
}

#' Converts Gene Expression Profiles (GEPs) to Module-EPs (MEPs).
#' @param geps A matrix of ranks where each row corresponds to a gene
#'     and each column to a perturbagen. Each column must include all
#'     ranks from 1 to the number of rows. Row and column names must
#'     be defined.
#' @param gmd A database of gene modules. See \code{importMsigDB.xml}
#'     for the format.
#' @param parallel If TRUE, gene sets will be processed in
#'     parallel. Requires a parallel backend.
#' @return A list of two matrices, one for Enrichment Scores and one
#'     for p-values. Each entry (i,j) refers to gene set i and
#'     perturbagen j.
#' @seealso buildMEPs
#' @export
gep2mep <- function(geps, gmd, parallel=F) {

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
            gsea(where, genematj, F)
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
#' @param rp A repository created by \code{buildEmptyDB}.
#' @return Vector of names.
#' @details The names are obtained from the repository entries that
#'     are tagged with "gmd" and removing the "_gmd" suffix.
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
       say("Can not create DB in existing folder", T) 
    } else rp <- repo_open(path, T)
    
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
                   paste(okargs, collapse=", ")), T)
    
    replace <- F
    
    if(rp$has(db_id)) {

        if(existing == "stop")
            say("MEP already exists", T)

        if(existing == "append") {
            say("Merging MEPs...")
            curpep <- rp$get(db_id)
            peps[["ES"]] <- cbind(curpep[["ES"]], peps[["ES"]])
            peps[["PV"]] <- cbind(curpep[["PV"]], peps[["PV"]])
            replace <- T
        }
        
        if(existing == "overwrite")
            replace <- T
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
           c("gep2mep", "meta"), replace=T)

    say("Done.")
}



#' Creates unique identifiers for gene module databases.
#' @param gmds A gene module database in the same format created by
#'     \code{importMsigDB.xml}.
#' @return A vector of identifiers created as "DB_SUBDB".
#' @seealso importMsigDB.xml
#' @export
buildDBids <- function(gmds) {
    dbs <- sapply(gmds, get, x="db")
    subdbs <- sapply(gmds, get, x="subdb")
    subdbs[subdbs==""] <- dbs[subdbs==""]
    db_ids <- paste(dbs, subdbs, sep="_")
    return(db_ids)
}

#' Build MEPs using an existing repository and stores them in it.
#' @param rp A gep2mep repository (see \code{buildEmptyDB}).
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
#' @export
buildMEPs <- function(rp, geps, parallel=F, existing="stop")
{
    okargs <- c("stop", "overwrite", "skip", "append")
    
    if(length(existing)>1 || (! existing %in% okargs))
        say(paste0("existing must be one of: ",
                   paste(okargs, collapse=", ")), T)
    
    dbs <- getDBlist(rp)

    for(i in 1:length(dbs))
    {        
        say(paste0("Working on DB: ", dbs[i],
                   " (", i, "/", length(dbs), ")" ))

        if(rp$has(dbs[i]))
            if(existing != "stop") {
                say(paste0("Existing DB found, action: ", existing))
            } else say(paste0("DB already present, check parameter: existing"), T)

        if(rp$has(dbs[i]) && existing == "skip")
            next
        
        thisdb <- rp$get(paste0(dbs[i], "_gmd"))
        peps <- gep2mep(geps, thisdb, parallel)

        storeMEPs(rp, dbs[i], peps, existing)

        cat("\n")        
    }
    
}

#' Shows repository statistics
#' @param path Path to the root directory of an existing repository.
#' @return Nothing, used for side effects.
#' @export
dbStats <- function(path="MEPDB") {
    rp <- repo_open(path)
    pert <- length(rp$get("perturbagens"))
    dbs <- getDBlist(rp)

    cat(paste("Number of perturbagens:", pert, "\n"))
    cat("Databases and their sizes:\n")
    for(i in 1:length(dbs))
        cat(paste0("   ", dbs[i], ":\t", length(rp$get(paste0(dbs[i], "_gmd"))), "\n"))
}


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
##                   "perturbagens in common, which will be treated as different."))
    
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
        pos[is.na(pos)] <- F
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


#' Export PertSEA results to XLS format
#' @param rp A repository created by \code{buildEmptyDB}.
#' @param results The output of \code{PertSEA}.
#' @param outname Name of the XLS file to be created.
#' @return Nothing.
#@' @export
exportPertSEA <- function(rp, results, outname="PertSEA.xls")
{
        if (requireNamespace("WriteXLS", quietly = TRUE)) {
            sheets <- attachModInfo(rp, results)
            names(sheets) <- gsub(":","_",names(sheets))
            WriteXLS::WriteXLS(sheets, outname, AutoFilter=T, BoldHeaderRow=T, FreezeRow=1)
        } else {
            stop("The suggested package WriteXLS is not installed.")
        }              
}

attachModInfo <- function(rp, results)
{
    dbs <- names(results[["PertSEA"]])
    newres <- list()
    for(i in 1:length(results[["PertSEA"]])) {
        db <- rp$get(paste0(dbs[i], "_gmd"))
        modnms <- rownames(results[["PertSEA"]][[i]])
        newres[[i]] <- cbind(
          Module = sapply(db[modnms], get, x="setname"),
          Description = sapply(db[modnms], get, x="desc"),
          results[["PertSEA"]][[i]],
          results[["details"]][[i]]
        )
    }
    names(newres) <- names(results[["PertSEA"]])
    return(newres)
}

#@' @export
PertSEA <- function(rp, pgset, bgset="all", details=T)
{
    dbs <- getDBlist(rp)
    if(length(bgset) == 1 && bgset=="all")
        bgset <- rp$get("perturbagens")
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
        ranked <- rankMEPsByRows(rp$get(dbs[i]), rankingset)
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

#@' @export
ModSEA <- function(rp, gmds, bgsets="all")
{
    dbs <- getDBlist(rp)
    if(! all(names(gmds) %in% dbs))
        say("Names of gmds should match gene module DB names.", T)

    dbs <- names(gmds)
    res <- list()
    for(i in 1:length(gmds))
    {
        say(paste0("Working on DB: ", dbs[i]))
        gmd <- gmds[[dbs[i]]]
        
        if(length(bgsets) == 1 && bgsets=="all") {
            bgset <- names(rp$get(paste0(dbs[i], "_gmd")))
        } else bgset <- bgsets[[i]]

        if(length(intersect(gmd, bgset)) > 0) {
            bgset <- setdiff(bgset, gmd)
            say("Common module sets removed from bgset.")
        }
        rankingset <- c(gmd, bgset)

        peps <- rp$get(dbs[i])
        notok <- rankingset[rankingset %in% rownames(peps)]
        if(length(notok)>0)
            say(paste0("Module set ids not found in ", dbs[i], ": ",
                       paste(notok, collapse=", ")), T)
        
        say(paste0("Col-ranking DB..."))
        ranked <- rankMEPsByCols(peps, rankingset)
        say(paste0("Computing enrichments..."))
        
        ks <- apply(ranked, 2, function(col) ks.test.2(col[gmd], col[bgset]))
        res[[dbs[i]]]$ES <- sapply(ks, get, x="ES")
        res[[dbs[i]]]$PV <- sapply(ks, get, x="p.value")
        say("done.")        
    }

    return(res)
}


gsea <- function(S, ranks_list, check=F, alternative = "two.sided")
#                ,leadedge = F)
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

 # return(list(ES=ks$ES, p=ks$"p.value", edge=ks$edge, lead=lead));
    return(list(ES=ks$ES, p=ks$"p.value", edge=ks$edge));
}


