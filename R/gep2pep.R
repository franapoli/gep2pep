
#'@import repo
#'@import XML
#'@import foreach


ks.test.2 <- function (x, y, signs=rep(1, length(x)+length(y)), ...,
                       alternative = c("two.sided", "less", "greater"),
                       exact = NULL, maxCombSize=10000, cumsum.return=F)
{
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 1L) 
        stop("not enough 'x' data")
    PVAL <- NULL
    if (is.numeric(y)) {
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        y <- y[!is.na(y)]
        n.x <- as.double(n)
        n.y <- length(y)
        if (n.y < 1L) 
            stop("not enough 'y' data")
        if (is.null(exact)) {
            exact <- (n.x * n.y < maxCombSize)
            if(!exact)
                warning(paste("P-value not computed exactly because",
                              "of combined sample size"))
        }
        METHOD <- "Two-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        n <- n.x * n.y/(n.x + n.y)
        w <- c(x, y)
        
        z <- cumsum(ifelse((order(w) <= n.x)*signs, 1/n.x, -1/n.y))
        if (length(unique(w)) < (n.x + n.y)) {
            if (exact) {
                warning("cannot compute exact p-value with ties")
                exact <- FALSE
            }
            else warning("p-value will be approximate in the presence of ties")
            z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
            TIES <- TRUE
        }
        STATISTIC <- switch(alternative, two.sided = max(abs(z)), 
            greater = max(z), less = -min(z))

        edge <- which.max(abs(z))
        ES <- z[edge]
        
        nm_alternative <- switch(alternative, two.sided = "two-sided", 
            less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y")
        if (exact && (alternative == "two.sided") && !TIES) 
            PVAL <- 1 - .Call(stats:::C_pSmirnov2x, STATISTIC, n.x, n.y)
    }
    else {
        stop("The case of is.character(y)=T has not been implemented")
        if (is.character(y)) 
            y <- get(y, mode = "function", envir = parent.frame())
        if (!is.function(y)) 
            stop("'y' must be numeric or a function or a string naming a valid function")
        METHOD <- "One-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        if (length(unique(x)) < n) {
            warning("ties should not be present for the Kolmogorov-Smirnov test")
            TIES <- TRUE
        }
        if (is.null(exact)) 
            exact <- (n < 100) && !TIES
        x <- y(sort(x), ...) - (0:(n - 1))/n
        STATISTIC <- switch(alternative, two.sided = max(c(x, 
            1/n - x)), greater = max(1/n - x), less = max(x))
        if (exact) {
            PVAL <- 1 - if (alternative == "two.sided")
                result = tryCatch({
                .C(C_pkolmogorov2x, p = as.double(STATISTIC), 
                  as.integer(n), PACKAGE = "stats")$p
                }, warning = function(w) {
                    warning(w)
                }, error = function(e) {
                    .Call(C_pKolmogorov2x, STATISTIC, n)
                }, finally = {
                })

            else {
                pkolmogorov1x <- function(x, n) {
                  if (x <= 0) 
                    return(0)
                  if (x >= 1) 
                    return(1)
                  j <- seq.int(from = 0, to = floor(n * (1 - 
                    x)))
                  1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 - 
                    x - j/n) + (j - 1) * log(x + j/n)))
                }
                pkolmogorov1x(STATISTIC, n)
            }
        }
        nm_alternative <- switch(alternative, two.sided = "two-sided", 
            less = "the CDF of x lies below the null hypothesis", 
            greater = "the CDF of x lies above the null hypothesis")
    }
    names(STATISTIC) <- switch(alternative, two.sided = "D", 
        greater = "D^+", less = "D^-")
    if (is.null(PVAL)) {
        pkstwo <- function(x, tol = 1e-06) {
            if (is.numeric(x)) 
                x <- as.double(x)
            else stop("argument 'x' must be numeric")
            p <- rep(0, length(x))
            p[is.na(x)] <- NA
            IND <- which(!is.na(x) & (x > 0))
            if (length(IND))
                p[IND] <- tryCatch({
                    tryRes <- .C(stats:::C_pkstwo, length(x[IND]), p = x[IND], 
                             as.double(tol), PACKAGE = "stats")$p
                }, warning = function(w) {
                    warning(w)
                }, error = function(e) {
                    tryRes <- .Call(stats:::C_pKS2, p = x[IND], tol)
                }, finally = {
                })
            p
        }
        PVAL <- ifelse(alternative == "two.sided", 1 - pkstwo(sqrt(n) * 
            STATISTIC), exp(-2 * n * STATISTIC^2))
    }
    PVAL <- min(1, max(0, PVAL))
    RVAL <- list(statistic = STATISTIC, p.value = PVAL,
                 alternative = nm_alternative, method = METHOD,
                 data.name = DNAME, ES = ES, edge = edge)
    if(cumsum.return)
        csum = RVAL <- c(RVAL, list(csum = z))

    class(RVAL) <- "htest"
    return(RVAL)
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
#' @details The format required by gep2pep for a database of gene sets
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

#' Converts Gene Expression Profiles (GEPs) to Pathway-EPs.
#' @param geps A matrix of ranks where each row corresponds to a gene
#'     and each column to a perturbagen. Each column must include all
#'     ranks from 1 to the number of rows. Row and column names must
#'     be defined.
#' @param pathw A database of gene sets. See \code{importMsigDB.xml}
#'     for the format.
#' @param parallel If TRUE, gene sets will be processed in
#'     parallel. Requires a parallel backend.
#' @return A list of two matrices, one for Enrichment Scores and one
#'     for p-values. Each entry (i,j) refers to gene set i and
#'     perturbagen j.
#' @seealso buildPEPs
#' @export
gep2pep <- function(geps, pathw, parallel=F) {
    
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
#' @param name Name of the repository. Defaults to \code {NULL} (a
#'     generic name will be given).
#' @param description Description of the repository. Defaults to \code {NULL}
#'     (a generic repository will be given).
#' @return An object of class \code{repo}.
#' @seealso buildPEPs
#' @export
buildEmptyDB <- function(path, gmd, name=NULL, description=NULL)
{
    if(is.null(name))
        name <- "gep2pep database"
    if(is.null(description))
        description <- paste("This database contains pathway information",
                             "and pathway expression profiles created by",
                             "the gep2pep package.")
    
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
               c("gep2pep", "gmd"))
    }    
    
    ## rp$put(unique(db_ids), "DB list",
    ##        "IDs of pathway sub-databases included in this repository",
    ##        c("gep2pep", "meta"))
        
    return(rp)
}


storePEPs <- function(rp, db_id, peps, existing) {

    okargs <- c("overwrite", "append", "stop")
    if(length(existing)>1 || (! existing %in% okargs))
        say(paste0("existing must be one of: ",
                   paste(okargs, collapse=", ")), T)
    
    replace <- F
    
    if(rp$has(db_id)) {

        if(existing == "stop")
            say("PEP already exists", T)

        if(existing == "append") {
            say("Merging peps...")
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
           c("gep2pep", "pep"), replace=replace)

    curids <- getDBlist(rp)
    
    rp$put(colnames(peps[[1]]), "perturbagens",
           "Names of perturbagens inducing expression profiles",
           c("gep2pep", "meta"), replace=T)
}


## addFromGEPs <- function(geps, parallel=F, path="PEPDB")
## {    
##     rp <- repo_open(path)
##     dbs <- rp$get("DB list")
##     for(i in 1:length(dbs)) {
##         say(paste0("Working on DB: ", dbs[i], " (", i, "/", length(dbs), ")" ))
##         thisdb <- rp$get(paste0(dbs[i], "_gmd"))
##         peps <- gep2pep(geps, thisdb, parallel)
##         curpep <- rp$get(dbs[i])
##         rp$set(dbs[i], obj=list(ESmat=cbind(curpep[["ESmat"]], peps[["ESmat"]]),
##                                 PVmat=cbind(curpep[["PVmat"]], peps[["PVmat"]])))
##     }
##     rp$set("perturbagens",
##            obj=c(rp$get("perturbagens"),
##                  colnames(peps[[1]])))
## }

#' @export
buildDBids <- function(pathdb) {
    dbs <- sapply(pathdb, get, x="db")
    subdbs <- sapply(pathdb, get, x="subdb")
    subdbs[subdbs==""] <- dbs[subdbs==""]
    db_ids <- paste(dbs, subdbs, sep="_")
    return(db_ids)
}

#' Build PEPs using an existing repository and stores them in it.
#' @param rp A gep2pep repository (see \code{buildEmptyDB}).
#' @param geps A matrix of ranks where each row corresponds to a gene
#'     and each column to a perturbagen. Each column must include all
#'     ranks from 1 to the number of rows. Row and column names must
#'     be defined.
#' @param parallel If TRUE, gene sets will be processed in
#'     parallel. Requires a parallel backend.
#' @param existing What to do if PEPs for a given DB of gene sets are
#'     already present. Can be one of:
#'
#' + stop: the default, throws an error. This is the safest approach.
#'
#' + skip: the existing DB will be skipped. This is especially useful
#' to build a repository incrementally by repeatedly using the same
#' call to \code{buildPEPs}.
#'
#' + overwrite: all the existing PEPs for the DB will be replaced by
#' new PEPs.
#'
#' + append: the new PEPs will be appended to the existing
#' ones. Useful to update the repository.
#'
#' @return Nothing. The computed PEPs will be available in the
#'     repository.
#' @seealso buildPEPs
#' @export
buildPEPs <- function(rp, geps, parallel=F, existing="stop")
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
        peps <- gep2pep(geps, thisdb, parallel)

        storePEPs(rp, dbs[i], peps, existing)

        cat("\n")        
    }
    
}

#' @param path Path to the root directory of an existing repository.
#' @return Nothing, used for side effects.
#' @export
dbStats <- function(path="PEPDB") {
    rp <- repo_open(path)
    pert <- length(rp$get("perturbagens"))
    dbs <- getDBlist(rp)

    cat(paste("Number of perturbagens:", pert, "\n"))
    cat("Databases and their sizes:\n")
    for(i in 1:length(dbs))
        cat(paste0("   ", dbs[i], ":\t", length(rp$get(paste0(dbs[i], "_gmd"))), "\n"))
}


## newPEPs <- function(rp, geps, gmd, id, parallel=F, overwrite=F)
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
#' @export
rankPEPsByCols <- function(peps, rankingset="all")
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

#' @export
rankPEPsByRows <- function(peps, rankingset="all")
{
    if(length(rankingset) == 1 && rankingset == "all")
        rankingset <- 1:ncol(peps[["ES"]])

    ESs <- peps[["ES"]][, rankingset]
    x <- t(apply(-ESs, 1, rank, ties.method = "random", na.last="keep"))
    return(x)
}

#@' @export
PertSEA <- function(rp, pgset, bgset="all")
{
    dbs <- getDBlist(rp)
    if(length(bgset) == 1 && bgset=="all")
        bgset <- rp$get("perturbagens")
    if(length(intersect(pgset, bgset))>0) {
        bgset <- setdiff(bgset, pgset)
        say("Common perturbagens removed from bgset.")
    }
    rankingset <- c(bgset, pgset)

    res <- list()
    for(i in 1:length(dbs)) {
        say(paste0("Working on DB: ", dbs[i]))
        say(paste0("Row-ranking DB..."))
        ranked <- rankPEPsByRows(rp$get(dbs[i]), rankingset)
        say(paste0("Computing enrichments..."))
        ks <- apply(ranked, 1, function(row) ks.test.2(row[pgset], row[bgset]))
        res[[dbs[i]]]$ES <- sapply(ks, get, x="ES")
        res[[dbs[i]]]$PV <- sapply(ks, get, x="p.value")
        say("done.")
    }
    
    return(res)
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
        ranked <- rankPEPsByCols(peps, rankingset)
        say(paste0("Computing enrichments..."))
        
        ks <- apply(ranked, 2, function(col) ks.test.2(col[gmd], col[bgset]))
        res[[dbs[i]]]$ES <- sapply(ks, get, x="ES")
        res[[dbs[i]]]$PV <- sapply(ks, get, x="p.value")
        say("done.")        
    }

    return(res)
}


gsea <- function(S, ranks_list, check=F, alternative = "two.sided",
                 leadedge = F)
{
    S <- S[!(is.na(S))]
    S1 <- ranks_list[S]
    S2 <- ranks_list[-S]

    if(length(S1)<1 || length(S2)<1 || all(is.na(S1)) || all(is.na(S2)))
      return(list(ES=NA, p=NA))

    ks <- ks.test.2(S1, S2, alternative=alternative, maxCombSize=10^10)

    lead=NA
    if(leadedge){
        sSm <- sort(Sm)
        if(ks$ES > 0)
            lead=sSm[sSm<=ks$edge] else lead=sSm[sSm>=ks$edge]
    }

    return(list(ES=ks$ES, p=ks$"p.value", edge=ks$edge, lead=lead));
}


if(F) {
## real data
    source("PEP.R")
    db <- importMsigDB.xml("data/in/msigdb_v6.0.xml")
    library(repo)
    rp <- repo_open()
    prls <- rp$get("Gmantra_PRL")
    library(doParallel)
    registerDoParallel(22)
    library(pbarETA)
    library(foreach)

    ## system("rm -r PEPDB")
    rp <- buildEmptyDB("PEPDB", db)
    buildPEPs(rp, prls, parallel=T, existing="stop")
}

if(F) {
## test data   
    source("PEP.R")
    
    library(foreach)
    system("rm -r TEST")
    testpws <- readRDS("testpws.RDS")    
    rp <- buildEmptyDB("TEST", testpws)
    testprl <- readRDS("testprl.RDS")
    buildPEPs(rp, testprl, parallel=F, existing="stop")
    
    ## CODE TO RESTART FROM A DB
    rp2 <- repo_open("PEPDB")
    dbids <- buildDBids(db)
    todo <- setdiff(unique(dbids), getDBlist(rp2))
    buildRepoDB(prls, db[dbids %in% todo], T, overwrite=T)

    testprl <- prls[1:500, 1:5]
    subdbs <- sapply(db, get, x="subdb")
    dbs <- sapply(db, get, x="db")
    subdbs[subdbs==""] <- dbs[subdbs==""]
    testpws <- do.call(c, lapply(unique(subdbs)[1:3], function(x) db[subdbs==x][1:10]))

    saveRDS(testprl, "testprl.RDS")
    saveRDS(testpws, "testpws.RDS")


    buildRepoDB(testprl[,1:3], testpws, T, overwrite=T)
    addGeps(testprl[,4:5], T)
    rp <- repo_open("PEPDB")
    pgset <- rp$get("perturbagens")[1:3]
    bgset <- "all"
    peps <- gep2pep(testprl, testpws)
    PGsea(rp, pgset)
}



