

## Workflow:
## library(GSEABase)
## library(devtools)
## library(testthat)
## load_all()

loadPEPs <- gep2pep:::.loadPEPs

dbfolder <- file.path(tempdir(), "gep2pepDB")

clear_test_repo <- function(suffix=NULL) {
    folder <- paste0(dbfolder, suffix)
    if(file.exists(folder))
        unlink(folder, TRUE)
}

create_test_repo <- function(suffix=NULL) {
    folder <- paste0(dbfolder, suffix)
    clear_test_repo(suffix)
    return(
        suppressMessages(
            createRepository(folder, testpws)
        )
    )    
}

testgep <- loadSampleGEP()
testpws <- loadSamplePWS()
testpws_old <- gep2pep:::convertFromGSetClass(testpws)

rp <- create_test_repo()
dbs <- makeCollectionIDs(testpws)
expected_dbs <- c("c3_TFT", "c3_MIR", "c4_CGN")


context("gep2pep: gep format checks")

randmat <- matrix(runif(9), 3)
nm_randmat <- randmat
rownames(nm_randmat) <- colnames(nm_randmat) <- 1:3
okmat <- apply(nm_randmat,2,rank)
okmatNA <- okmat; okmatNA[okmatNA==2] <- NA
okmatRep <- okmat; okmatRep[2:3, 1] <- 1
test_that("check geps", {
    expect_error(checkGEPsFormat(randmat))
    expect_error(checkGEPsFormat(nm_randmat))
    expect_silent(checkGEPsFormat(okmatNA))
    expect_error(checkGEPsFormat(okmatRep))
    expect_silent(checkGEPsFormat(okmat))
})

context("gep2pep: creation of newdb")


test_that("new db creation", {
    expect_error(createRepository(".", testgep))
    expect_equal(length(dbs), length(testpws))
    expect_true(setequal(unique(dbs), expected_dbs))
    expect_true(setequal(getCollections(rp), expected_dbs))
    expect_true(rp$has(paste0(expected_dbs[1], "_sets")))
    expect_true(rp$has(paste0(expected_dbs[2], "_sets")))
    expect_true(rp$has(paste0(expected_dbs[3], "_sets")))
    expect_true(rp$has("gep2pep repository"))
    expect_equal(length(rp$entries()), 5)
    expect_equal(length(rp$get(paste0(expected_dbs[1], "_sets"))), 10)
    expect_equal(length(rp$get(paste0(expected_dbs[2], "_sets"))), 10)
    expect_equal(length(rp$get(paste0(expected_dbs[3], "_sets"))), 10)
    expect_failure(expect_warning(suppressMessages(checkRepository(rp))))
    expect_error(loadESmatrix(rp, "random name"))
    expect_error(loadESmatrix(rp, "c3_TFT"))
    expect_error(loadPVmatrix(rp, "random name"))
    expect_error(loadPVmatrix(rp, "c3_TFT"))
})

context("gep2pep: creation of SGE db")

addSingleGeneSets(rp, rownames(testgep))
test_that("SGE were created", {
    expect_true("SGE_sets" %in% names(rp$entries()))
})

context("gep2pep: creation of peps")

suppressMessages(
    buildPEPs(rp, testgep, progress_bar=FALSE,
              min_size=3)
)

test_that("build first PEPs", {
    expect_equal(length(rp$entries()), 10)
    expect_equal(length(dbs), length(testpws))
    expect_true(rp$has(expected_dbs[1]))
    expect_true(rp$has(expected_dbs[2]))
    expect_true(rp$has(expected_dbs[3]))
    expect_equal(names(rp$get(expected_dbs[1])), c("ES", "PV"))
    expect_equal(names(rp$get(expected_dbs[3])), c("ES", "PV"))
    expect_equal(nrow(rp$get(expected_dbs[1])[[1]]), 10)
    expect_equal(nrow(rp$get(expected_dbs[1])[[2]]), 10)
    expect_equal(nrow(rp$get(expected_dbs[3])[[1]]), 10)
    expect_equal(nrow(rp$get(expected_dbs[3])[[2]]), 10)
    expect_equal(ncol(rp$get(expected_dbs[1])[[1]]), ncol(testgep))
    expect_equal(ncol(rp$get(expected_dbs[1])[[2]]), ncol(testgep))
    expect_equal(ncol(rp$get(expected_dbs[3])[[1]]), ncol(testgep))
    expect_equal(ncol(rp$get(expected_dbs[3])[[2]]), ncol(testgep))
    expect_failure(expect_warning(suppressMessages(checkRepository(rp))))
    expect_error(loadESmatrix(rp, "random name"))
    expect_equal(loadESmatrix(rp, "c3_TFT"), rp$get("c3_TFT")$ES)
    expect_error(loadPVmatrix(rp, "random name"))
    expect_equal(loadPVmatrix(rp, "c3_TFT"), rp$get("c3_TFT")$PV)
})



res <- list()
for(i in 1:3) {
    testi <- sample(1:length(testpws),1)
    testj <- sample(1:ncol(testgep),1)
    set <- testpws_old[[testi]]$set
    id <- testpws_old[[testi]]$id
    tomatch <- intersect(rownames(testgep), set)
    inset <- testgep[match(tomatch, rownames(testgep)), testj]
    if(length(tomatch) >= 3) {
        ks <- ks.test.2(inset, (1:nrow(testgep))[-inset], maxCombSize=10^10)
    } else ks <- list(ES=as.numeric(NA), p.value=as.numeric(NA))
    dbi <- dbs[testi]
    res[[i]] <- list(id=id, testj=testj, ks=ks, dbi=dbi)
}

s <- sample(1:ncol(testgep), 1)
g <- sample(1:nrow(testgep), 1)
ksSGE <- ks.test.2(testgep[g,s], (1:nrow(testgep))[-testgep[g,s]])

test_that("KS statistics", {
    i <- 1
    id <- res[[i]]$id; testj <- res[[i]]$testj; ks <- res[[i]]$ks; dbi <- res[[i]]$dbi
    expect_equal(loadPEPs(rp, dbi)$ES[id, testj], ks$ES)
    expect_equal(loadPEPs(rp, dbi)$PV[id, testj], ks$p.value)
    i <- 2
    id <- res[[i]]$id; testj <- res[[i]]$testj; ks <- res[[i]]$ks; dbi <- res[[i]]$dbi
    expect_equal(loadPEPs(rp, dbi)$ES[id, testj], ks$ES)
    expect_equal(loadPEPs(rp, dbi)$PV[id, testj], ks$p.value)
    i <- 3
    id <- res[[i]]$id; testj <- res[[i]]$testj; ks <- res[[i]]$ks; dbi <- res[[i]]$dbi
    expect_equal(loadPEPs(rp, dbi)$ES[id, testj], ks$ES)
    expect_equal(loadPEPs(rp, dbi)$PV[id, testj], ks$p.value)
    expect_equal(ksSGE$ES, loadESmatrix(rp, "SGE")[g,s])
    expect_equal(ksSGE$p.value, loadPVmatrix(rp, "SGE")[g,s])
})

context("gep2pep: adding existing peps")

oldTFT <- loadPEPs(rp, "c3_TFT")
test_that("Adding PEPs", {
    expect_warning(
        suppressMessages(buildPEPs(rp, testgep[, 1:3], progress_bar=FALSE))
    )
    expect_failure(expect_warning(suppressMessages(checkRepository(rp))))    
})

untouchedTFT <- oldTFT

subs <- c(2,4,5)
smallTFT <- list(ES=oldTFT$ES[, subs],
                 PV=oldTFT$PV[, subs])
## the following will create conflicts with
## the "perturbagens" item
rp$set("c3_TFT", smallTFT)

rebuiltTFT <- loadPEPs(rp, "c3_TFT")
test_that("Adding PEPs", {
    expect_warning(
        suppressMessages(buildPEPs(rp, testgep[, 1:3], progress_bar=FALSE))
    )
    expect_equal(untouchedTFT, oldTFT)
    expect_equal(rebuiltTFT$ES, oldTFT$ES[,colnames(rebuiltTFT$ES)])
    expect_equal(rebuiltTFT$PV, oldTFT$PV[,colnames(rebuiltTFT$PV)])
    expect_warning(suppressMessages(checkRepository(rp)))
})

rp$set("c3_TFT", oldTFT)


context("gep2pep: adding peps one by one")

rp2 <- create_test_repo("2")
for(i in 1:ncol(testgep))
    suppressMessages(
        buildPEPs(rp2, testgep[,i,drop=FALSE], progress_bar=FALSE)
    )

rptft <- rp$get("c3_TFT")
rptft2 <- rp2$get("c3_TFT")
eseq <- rptft$ES == rptft2$ES
pveq <- rptft$PV == rptft2$PV
test_that("adding one by one", {
    expect_true(all(eseq[!is.na(eseq)]))
    expect_true(all(pveq[!is.na(pveq)]))
    expect_true(all(colnames(rptft) == colnames(rptft2)))
    expect_true(all(rownames(rptft) == rownames(rptft2)))    
    expect_failure(expect_warning(suppressMessages(checkRepository(rp2))))
})



context("gep2pep: Row and Col ranking")

peps1 <- rp$get(expected_dbs[1])
peps3 <- rp$get(expected_dbs[3])

es1 <- peps1$ES
es3 <- peps3$ES
RowRanked1 <- rankPEPsByRows.ES(peps1)
nas1 <- which(is.na(RowRanked1[,1]))
RowRanked3 <- rankPEPsByRows.ES(peps3)
nas3 <- which(is.na(RowRanked3[,3]))

test_that("Row ranking", {
    expect_true(all(apply(RowRanked1[-nas1,], 1, setequal, 1:5)))
    expect_true(all(apply(RowRanked3[-nas3,], 1, setequal, 1:5)))
    expect_equal(RowRanked1[1,1], 5)
    expect_equal(RowRanked1[4,3], 1)
    expect_equal(RowRanked3[5,4], 5)
    expect_equal(RowRanked3[8,3], 1)
})


ColRanked1 <- rankPEPsByCols.SPV(peps1)
ColRanked3 <- rankPEPsByCols.SPV(peps3)

randj <- sample(ncol(ColRanked3),1)
PVs <- peps1$PV[,randj]
ESs <- peps1$ES[,randj]
if(any(ESs<0)) {
    PVs[ESs>=0] <- Inf
    lastid <- which.min(PVs)
} else lastid <- which.max(PVs)

test_that("Column ranking with SPV", {
    expect_true(all(apply(ColRanked1[-nas1,], 2, setequal, 1:9)))
    expect_true(all(apply(ColRanked3[-nas3,], 2, setequal, 1:9)))
    expect_equal(ColRanked1[2,1], 1)
    expect_equal(ColRanked1[1,1], 9)
    expect_equal(ColRanked3[8,3], 1)
    expect_equal(ColRanked3[10,3], 9)
    expect_equal(ColRanked3[4,5], 9)
    expect_equal(ColRanked1[lastid, randj], 9)
})


ColRanked3 <- rankPEPsByCols.NES(peps3)
manual <- t(scale(t(peps3$ES)))
stopifnot(all(apply(manual[-nas3,],1,sd)-1<10^-15))
manualR <- apply(-manual, 2, rank, na.last="keep")

test_that("Column ranking with NES", {
    expect_true(all(ColRanked3[-nas3,]==manualR[-nas3,]))
    expect_true(all(is.na(manualR[nas3,])))
})


context("gep2pep: CondSEA")
pgset <- c("(+)_chelidonine",  "(+/_)_catechin")
suppressMessages(
    res <- CondSEA(rp, pgset)
)
randi <- sample(1:length(testpws), 1)
pwsid <- testpws_old[[randi]]$id
randDB <- dbs[randi]
ranked <- rankPEPsByRows.ES(loadPEPs(rp, randDB))
inset <- ranked[pwsid, pgset]
outset <- ranked[pwsid, setdiff(colnames(ranked), pgset)]
if(length(inset[!is.na(inset)])>0 &&
   length(outset[!is.na(outset)])>0) {
    ks <- ks.test.2(inset, outset, maxCombSize=10^10)
} else {
    ks <- list(ES=as.numeric(NA), p.value=as.numeric(NA))
}
test_that("CondSEA", {
    expect_equal(getDetails(res, "c3_TFT"), res$details[["c3_TFT"]])   
    expect_equal(getResults(res, "c3_TFT"), res$CondSEA[["c3_TFT"]])
    expect_equal(unname(res$CondSEA[[randDB]][pwsid, "ES"]),
                 ks$ES)
    expect_equal(unname(res$CondSEA[[randDB]][pwsid, "PV"]),
                 ks$p.value)
})


context("gep2pep: PathSEA")

db1 <- expected_dbs[1]
db3 <- expected_dbs[3]
pws1 <- sapply(testpws[makeCollectionIDs(testpws)==db1][c(2,5,6,9)], setName)
pws3 <- sapply(testpws[makeCollectionIDs(testpws)==db3][c(1,3,10)], setName)
res <- suppressMessages(PathSEA(rp, testpws[c(pws1, pws3)]))
setids1 <- sapply(testpws[pws1], setIdentifier)
setids3 <- sapply(testpws[pws3], setIdentifier)
randj1 <- sample(1:ncol(testgep), 1)
ranked <- rankPEPsByCols.SPV(rp$get(db1))
peps <- rp$get(db1)
inset <- ranked[setids1, randj1]
outset <- ranked[setdiff(rownames(ranked), setids1), randj1]
inset <- inset[!is.na(inset)]
outset <- outset[!is.na(outset)]
if(length(inset)>0 &&
   length(outset)>0) {
    ks1 <- ks.test.2(inset, outset, maxCombSize=10^10)
} else {
    ks1 <- list(ES=as.numeric(NA), p.value=as.numeric(NA))
}
randj3 <- sample(1:ncol(testgep), 1)
ranked <- rankPEPsByCols.SPV(rp$get(db3))
peps <- rp$get(db3)
inset <- ranked[setids3, randj3]
outset <- ranked[setdiff(rownames(ranked), setids3), randj3]
inset <- inset[!is.na(inset)]
outset <- outset[!is.na(outset)]
if(length(inset)>0 &&
   length(outset)>0) {
    ks3 <- ks.test.2(inset, outset, maxCombSize=10^10)
} else {
    ks3 <- list(ES=as.numeric(NA), p.value=as.numeric(NA))
}
name1 <- colnames(testgep)[randj1]
name3 <- colnames(testgep)[randj3]
test_that("PathSEA", {
    expect_equal(getDetails(res, "c3_TFT"), res$details[["c3_TFT"]])
    expect_equal(getResults(res, "c3_TFT"), res$PathSEA[["c3_TFT"]])
    expect_equal(unname(res[["PathSEA"]][[db1]][name1, "ES"]),
                 ks1$ES)
    expect_equal(unname(res[["PathSEA"]][[db1]][name1, "PV"]),
                 ks1$p.value)
    expect_equal(unname(res[["PathSEA"]][[db3]][name3, "ES"]),
                 ks3$ES)
    expect_equal(unname(res[["PathSEA"]][[db3]][name3, "PV"]),
                 ks3$p.value)
})

subset <- c("(+)_isoprenaline", "(_)_atenolol")
res2 <- suppressMessages(PathSEA(rp, testpws[c(pws1, pws3)],
                                 subset=subset))
test_that("PathSEA on subset of PEPs", {
    expect_equal(getDetails(res, "c3_TFT")[,subset],
                 getDetails(res2, "c3_TFT")[,subset])
})

## A gene that is found in at least 3 pathways:
gene <- intersect(intersect(geneIds(testpws[[3]]), geneIds(testpws[[4]])),
                  geneIds(testpws[[7]]))[1]
subpw <- gene2pathways(rp, gene)
## it is actually found also in a 4th:
extrapw <- setdiff(sapply(subpw, setIdentifier),
                   sapply(testpws[c(3,4,7)], setIdentifier))

test_that("gene2pathways", {
    expect_true(length(subpw) >= 3)
    expect_true(all(sapply(testpws[c(3,4,7)], setIdentifier)
                    %in% sapply(subpw, setIdentifier)))
    expect_true(gene %in%
                geneIds(testpws[sapply(testpws, setIdentifier) == extrapw])
                [[1]])
})



context("gep2pep: missing pathway points")

oldpeps <- rp$get("c3_TFT")

peps <- oldpeps
nas <- sample(1:length(peps$ES), 5)
peps$ES[nas] <- peps$PV[nas] <- NA
rp$set("c3_TFT", obj=peps)

pgset <- sample(colnames(ranked), 2)
bgset <- sample(setdiff(colnames(ranked), pgset), 2)

ranked <- t(apply(-peps$ES[,c(pgset,bgset)], 1,
                        rank, na.last="keep"))

suppressMessages(
    rankedByFunc <- .getPEPRanks(rp, "c3_TFT", c(pgset,bgset),
                                 subcols=c(pgset,bgset),
                                 usecache=FALSE,
                                 rankingFun=rankPEPsByRows.ES)
)

tests <- rankedByFunc==ranked
test_that("ranking with random missing points", {
    expect_true(all(tests[!is.na(tests)]))
})

ksES <- ksPV <- vector("numeric")
for(i in 1:nrow(ranked)) {
    inset <- ranked[i, pgset]
    inset <- inset[!is.na(inset)]
    outset <- ranked[i, bgset]
    outset <- outset[!is.na(outset)]
    if(length(inset)>0 &&
       length(outset)>0) {
        ks <- ks.test.2(inset, outset, maxCombSize=10^10)
    } else {
        ks <- list(ES=as.numeric(NA), p.value=as.numeric(NA))
    }
    ksES[[i]] <- ks$ES
    ksPV[[i]] <- ks$p.value
}

suppressMessages(
    csea <- CondSEA(rp, pgset, bgset,
                    "c3_TFT", sortoutput=FALSE)$CondSEA$c3_TFT
)

test_that("condsea with missing pathways", {
    expect_equal(ksES, csea$ES)
    expect_equal(ksPV, csea$PV)
})

rp$set("c3_TFT", obj=oldpeps)




context("gep2pep: XLS export")

pgset <- c("(+)_isoprenaline", "(_)_atenolol")
res1 <- suppressMessages(CondSEA(rp, pgset, details=FALSE))
res2 <- suppressMessages(CondSEA(rp, pgset, details=FALSE))
out1 <- tempfile()
out2 <- tempfile()
exportSEA(rp, res1, out1)
exportSEA(rp, res2, out2)

test_that("excel files produced", {
    expect_true(file.exists(out1))
    expect_true(file.exists(out2))
})
