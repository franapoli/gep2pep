

## ## Workflow:
## library(GSEABase)
## library(devtools)
## library(testthat)
## load_all()

loadPEPs <- gep2pep::.loadPEPs

dbfolder <- file.path(tempdir(), "gep2pepDB")

clear_test_repo <- function(suffix=NULL) {
    folder <- paste0(dbfolder, suffix)
    if(file.exists(folder))
        unlink(folder, T)
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
testpws <- as.CategorizedCollection(
    loadSamplePWS()
)
testpws_old <- gep2pep:::convertFromGSetClass(testpws)

rp <- create_test_repo()
dbs <- makeCollectionIDs(testpws)
expected_dbs <- c("c3_TFT", "c3_MIR", "c4_CGN")


suppressMessages(buildPEPs(rp, testgep, progress_bar=FALSE))

context("creation of RAW peps")

suppressMessages(
    buildPEPs(rp, testgep[,1:2], progress_bar=FALSE,
              rawmode_id=1)
)
suppressMessages(
    buildPEPs(rp, testgep[,3:5], progress_bar=FALSE,
              rawmode_id=2)
)

outfiles1 <- paste0(getCollections(rp), "#1.RDS")
outfiles2 <- paste0(getCollections(rp), "#2.RDS")
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
rp$rm(tags="pep", force=T)
importFromRawMode(rp)

pep2 <- loadPEPs(rp, colls[2])
w <- is.na(pep2$ES[,1])
    
test_that("check hdf5 PEPss", {
    expect_true(all(oldpep2$ES[!w,]==pep2$ES[!w,]))
    expect_true(all(oldpep2$PV[!w,]==pep2$PV[!w,]))
    expect_true(all(is.na(oldpep2$PV[w,])))
    expect_true(all(rownames(oldpep2$ES)==rownames(pep2$ES)))
    expect_true(all(rownames(oldpep2$PV)==rownames(pep2$PV)))
    expect_true(all(colnames(oldpep2$ES)==colnames(pep2$ES)))
    expect_true(all(colnames(oldpep2$PV)==colnames(pep2$PV)))
})

