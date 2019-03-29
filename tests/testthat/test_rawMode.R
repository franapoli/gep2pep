

## ## Workflow:
## library(GSEABase)
## library(devtools)
## library(testthat)
## library(digest)
## load_all()

context("rawMode: creation of RAW peps")

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
## adding problematic characters
rp$set("c3_MIR_sets", newname="c3_M:I*R_sets")

dbs <- makeCollectionIDs(testpws)
nms <- colnames(testgep)
suppressMessages(buildPEPs(rp, testgep, progress_bar=FALSE))
suppressMessages(origCondSEA <- CondSEA(rp, nms[c(1,3)], nms[1:4],
                                        usecache=TRUE))



suppressMessages(
    buildPEPs(rp, testgep[,1:2], progress_bar=FALSE,
              rawmode_id=1)
)
suppressMessages(
    buildPEPs(rp, testgep[,3:5], progress_bar=FALSE,
              rawmode_id=2)
)

outfiles1 <- gsub("[^[:alnum:]]", "_", getCollections(rp))
outfiles2 <- outfiles1
outfiles1 <- paste0(outfiles1, "#1.RDS")
outfiles2 <- paste0(outfiles2, "#2.RDS")
outfiles <- c(outfiles1, outfiles2)
outdir <- file.path(rp$root(), "raw")
f1 <- readRDS(paste0(file.path(outdir, outfiles[1])))

test_that("build hdf5 PEPs", {
    expect_true(all(sapply(outfiles, `%in%`, list.files(outdir))))
    expect_equal(f1$ES[,1], rp$get("c3_TFT")$ES[,1])
    expect_equal(f1$PV[,1], rp$get("c3_TFT")$PV[,1])
    expect_equal(f1$ES[,2], rp$get("c3_TFT")$ES[,2])
    expect_equal(f1$PV[,2], rp$get("c3_TFT")$PV[,2])
})


colls <- getCollections(rp)

oldpep2 <- rp$get(colls[2])
rp$rm(tags="pep", force=TRUE)
suppressMessages(
    importFromRawMode(rp)
    )

pep2 <- loadPEPs(rp, colls[2])
w <- is.na(pep2$ES[,1])

subset <- c(3,5)
subsetC <- colnames(pep2$ES)[subset]

test_that("check hdf5 PEPs", {
    expect_equal(oldpep2$ES[!w,], pep2$ES[!w,])
    expect_equal(oldpep2$PV[!w,], pep2$PV[!w,])
    expect_true(all(is.na(oldpep2$PV[w,])))
    expect_equal(loadPEPs(rp, colls[2], subset)$ES[!w,],
                 oldpep2$ES[!w, subset])
    expect_equal(loadPEPs(rp, colls[2], subsetC)$PV[!w,],
                 oldpep2$PV[!w, subsetC])
})


context("rawMode: Merging in raw mode")

newroot <- file.path(tempdir(), "mergedRAW")
if(file.exists(newroot))
    unlink(newroot, TRUE)

nms <- colnames(testgep)
pgA <- sample(nms, sample(1:2, 1))
pgB <- sample(setdiff(nms, pgA), sample(1:2, 1))
pgC <- setdiff(nms, c(pgA,pgB))
mergestr <- list(B=pgB, A=pgA, C=pgC)


suppressMessages(
    createMergedRepository(rp$root(), newroot, mergestr,
                           progressBar=FALSE)
)

suppressMessages(
    rpM <- openRepository(newroot)
    )

peps <- loadPEPs(rp, "c3_M:I*R")

chis <- function(ps) {
    S=-2*sum(log(ps))
    cump <- pchisq(S, 2*length(ps))
    return(cump)
}

manualMergedES <- matrix(NA, 10, 3)
manualMergedPV <- matrix(NA, 10, 3)
for(i in 1:length(mergestr)) {
    manualMergedES[,i] <- apply(peps$ES[,mergestr[[i]],drop=FALSE], 1, mean)
    manualMergedPV[,i] <- apply(peps$PV[,mergestr[[i]],drop=FALSE], 1, chis)
    }
colnames(manualMergedES) <-
    colnames(manualMergedPV) <-
    c("B", "A", "C")
rownames(manualMergedES) <-
    rownames(manualMergedPV) <-
    rownames(peps$ES)


mergedPEPs <- loadPEPs(rpM, "c3_M:I*R")
test_that("merged repository in HDF5", {
    expect_equal(length(rp$entries()), length(rpM$entries())+3-1)
    expect_equal(colnames(mergedPEPs$ES), c("B", "A", "C"))
    expect_equal(colnames(mergedPEPs$PV), c("B", "A", "C"))
    for(i in 1:3) {
        expect_equal(mergedPEPs$ES[,i], manualMergedES[,i])
        expect_equal(mergedPEPs$PV[,i], manualMergedPV[,i])
    }
})



context("rawMode: Ranking in raw mode")

peps1 <- loadPEPs(rp, getCollections(rp)[1])
peps3 <- loadPEPs(rp, getCollections(rp)[3])

es1 <- peps1$ES
es3 <- peps3$ES
RowRanked1 <- rankPEPsByRows.ES(peps1)
nas1 <- which(is.na(RowRanked1[,1]))
RowRanked3 <- rankPEPsByRows.ES(peps3)
nas3 <- which(is.na(RowRanked3[,3]))

test_that("Row ES-ranking in raw mode", {
    expect_true(all(apply(RowRanked1[-nas1,], 1, setequal, 1:5)))
    expect_true(all(apply(RowRanked3[-nas3,], 1, setequal, 1:5)))
    expect_equal(RowRanked1[1,1], 5)
    expect_equal(RowRanked1[4,3], 1)
    expect_equal(RowRanked3[5,4], 5)
    expect_equal(RowRanked3[8,3], 1)
})

ColRanked1 <- rankPEPsByCols.ES(peps1)
nas1 <- which(is.na(ColRanked1[,1]))
ColRanked3 <- rankPEPsByCols.ES(peps3)
nas3 <- which(is.na(ColRanked3[,3]))

test_that("Row ES-ranking in raw mode", {
    expect_true(all(apply(ColRanked1[-nas1,], 2, setequal, 1:9)))
    expect_true(all(apply(ColRanked3[-nas3,], 2, setequal, 1:9)))
})



subrank <- suppressMessages(
    .getPEPRanks(rp, getCollections(rp)[1],
                 nms[c(1,3,4)],
                 subcols=nms[c(1,4)],
                 rankingFun = rankPEPsByRows.ES,
                 usecache=TRUE)
)
nas <- !is.na(subrank[,1])
manualranks <- t(apply(-peps1$ES[, c(1,3,4)], 1,
                       rank, na.last="keep"))[,c(1,3)]
dig <- digest(list(getCollections(rp)[1],
                   sort(nms[c(1,3,4)]), "rankPEPsByRows.ES"))

test_that("Getting rank subsets", {
    expect_equal(subrank[nas,], manualranks[nas,])
    expect_true(rp$has(dig))
})

suppressMessages(
    newCondSEA <- CondSEA(rp, nms[1:2], nms[1:4], usecache=TRUE)
)
test_that("CondSEA in raw mode", {
    identical(origCondSEA, newCondSEA)
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

test_that("Column ranking with SPV in raw mode", {
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

test_that("Column ranking with NES in raw mode", {
    expect_true(all(ColRanked3[-nas3,]==manualR[-nas3,]))
    expect_true(all(is.na(manualR[nas3,])))
})



## context("gep2pep: adding peps one by one")

## rp2 <- create_test_repo("2")
## for(i in 1:2)
##     suppressMessages(
##         buildPEPs(rp2, testgep[,i,drop=F], progress_bar=FALSE,
##                   rawmode_id = i)
##     )
## importFromRawMode(rp2)

## for(i in 3:5)
##     suppressMessages(
##         buildPEPs(rp2, testgep[,i,drop=F], progress_bar=FALSE,
##                   rawmode_id = i)
##     )
## importFromRawMode(rp2)

## rptft <- rp$get("c3_TFT")
## rptft2 <- rp2$get("c3_TFT")
## eseq <- rptft$ES == rptft2$ES
## pveq <- rptft$PV == rptft2$PV
## test_that("adding one by one in raw mode", {
##     expect_true(all(eseq[!is.na(eseq)]))
##     expect_true(all(pveq[!is.na(pveq)]))
##     expect_true(all(colnames(rptft) == colnames(rptft2)))
##     expect_true(all(rownames(rptft) == rownames(rptft2)))    
##     expect_failure(expect_warning(suppressMessages(checkRepository(rp2))))
## })
