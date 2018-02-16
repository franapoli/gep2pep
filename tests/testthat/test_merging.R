## ## Workflow:
## library(GSEABase)
## library(devtools)
## library(testthat)
## load_all()

loadPEPs <- gep2pep:::.loadPEPs
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

testpws <- loadSamplePWS()


rp <- create_test_repo()

testgep <- loadSampleGEP()
suppressMessages(buildPEPs(rp, testgep, progress_bar=FALSE))


context("merging: Merging")

peps <- loadPEPs(rp, "c3_MIR")

nms <- colnames(testgep)
pgA <- sample(nms, sample(1:2, 1))
pgB <- sample(setdiff(nms, pgA), sample(1:2, 1))
pgC <- setdiff(nms, c(pgA,pgB))
mergestr <- list(B=pgB, A=pgA, C=pgC)

mergedPEPs <- suppressMessages(
    .mergePEPs(rp, "c3_MIR", mergestr,
               progressBar=F)
    )

chis <- function(ps) {
    S=-2*sum(log(ps))
    cump <- pchisq(S, 2*length(ps))
    return(cump)
}

manualMergedES <- matrix(NA, 10, 3)
manualMergedPV <- matrix(NA, 10, 3)
for(i in 1:length(mergestr)) {
    manualMergedES[,i] <- apply(peps$ES[,mergestr[[i]],drop=F], 1, mean)
    manualMergedPV[,i] <- apply(peps$PV[,mergestr[[i]],drop=F], 1, chis)
    }
colnames(manualMergedES) <-
    colnames(manualMergedPV) <-
    c("B", "A", "C")
rownames(manualMergedES) <-
    rownames(manualMergedPV) <-
    rownames(peps$ES)

test_that("merging function", {
    expect_equal(colnames(mergedPEPs$ES), c("B", "A", "C"))
    expect_equal(colnames(mergedPEPs$PV), c("B", "A", "C"))
    for(i in 1:3) {
        expect_equal(mergedPEPs$ES[,i], manualMergedES[,i])
        expect_equal(mergedPEPs$PV[,i], manualMergedPV[,i])
    }
})



newroot <- file.path(tempdir(), "merged")
if(file.exists(newroot))
    unlink(newroot, T)
suppressMessages(
    createMergedRepository(rp$root(), newroot, mergestr,
                           progressBar=FALSE)
)
suppressMessages(
    rpM <- openRepository(newroot)
    )

mergedPEPs <- rpM$get("c3_MIR")
test_that("merged repository", {
    expect_equal(length(rp$entries()), length(rpM$entries())-1)
    expect_equal(colnames(mergedPEPs$ES), c("B", "A", "C"))
    expect_equal(colnames(mergedPEPs$PV), c("B", "A", "C"))
    for(i in 1:3) {
        expect_equal(mergedPEPs$ES[,i], manualMergedES[,i])
        expect_equal(mergedPEPs$PV[,i], manualMergedPV[,i])
    }
})



if(file.exists(newroot))
    unlink(newroot, T)
suppressMessages(

    createMergedRepository(rp$root(), newroot, mergestr,
                           progressBar=FALSE,
                           mergeFunc="CondSEA",
                           collections="c3_MIR")

)

suppressMessages({
    unmerged <- rankPEPsByRows.ES(loadPEPs(rp, "c3_MIR"))
    merged <- loadPEPs(openRepository(newroot), "c3_MIR")
})

wrow <- sample(nrow(unmerged), 1)
wmer <- sample(1:length(mergestr), 1)
pgset <- mergestr[[wmer]]
inset <- unmerged[wrow, pgset]
inset <- inset[!is.na(inset)]
outset <- (1:5)[-inset]
if(length(inset)>0 && length(outset)>0) {
    ks <- ks.test.2(inset, outset)
    } else ks <- list(ES=as.numeric(NA), p.value=as.numeric(NA))
test_that("CondSEA merging", {
    expect_equal(length(rp$entries()), length(rpM$entries())-1)
    expect_equal(colnames(merged$ES), c("B", "A", "C"))
    expect_equal(colnames(merged$PV), c("B", "A", "C"))
    expect_equal(ks$ES, merged$ES[wrow, wmer])
    expect_equal(ks$p.value, merged$PV[wrow, wmer])
})
