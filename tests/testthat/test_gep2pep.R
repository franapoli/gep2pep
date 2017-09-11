
context("gep2pep")

testgep <- readRDS(system.file("testgep.RDS", package="gep2pep"))
testgep <- testgep[rownames(testgep) != "NA",]
testgep <- apply(testgep, 2, rank)
testpws <- readRDS(system.file("testgmd.RDS", package="gep2pep"))
dbfolder <- file.path(tempdir(), "gep2pepDB")
if(file.exists(dbfolder))
    unlink(dbfolder, T)
rp <- createRepository(dbfolder, testpws)


dbs <- makeCollectionIDs(testpws)
expected_dbs <- c("C3_TFT", "C3_MIR", "C4_CGN")

test_that("new db creation", {
  expect_equal(length(dbs), length(testpws))
  expect_true(setequal(unique(dbs), expected_dbs))
  expect_true(setequal(getCollections(rp), expected_dbs))
  expect_true(rp$has(paste0(expected_dbs[1], "_sets")))
  expect_true(rp$has(paste0(expected_dbs[2], "_sets")))
  expect_true(rp$has(paste0(expected_dbs[3], "_sets")))
  expect_true(rp$has("gep2pep database"))
  expect_equal(length(rp$entries()), 4)
  expect_equal(length(rp$get(paste0(expected_dbs[1], "_sets"))), 10)
  expect_equal(length(rp$get(paste0(expected_dbs[2], "_sets"))), 10)
  expect_equal(length(rp$get(paste0(expected_dbs[3], "_sets"))), 10)
})

buildPEPs(rp, testgep)

test_that("build first PEPs", {
  expect_equal(length(rp$entries()), 8)
  expect_equal(length(dbs), length(testpws))
  expect_true(rp$has(expected_dbs[1]))
  expect_true(rp$has(expected_dbs[2]))
  expect_true(rp$has(expected_dbs[3]))
  expect_true(rp$has("perturbagens"))

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
})


res <- list()
for(i in 1:3) {
  testi <- sample(1:length(testpws),1)
  testj <- sample(1:ncol(testgep),1)
  set <- testpws[[testi]]$set
  id <- testpws[[testi]]$id
  tomatch <- intersect(rownames(testgep), set)
  inset <- testgep[match(tomatch, rownames(testgep)), testj]
  ks <- ks.test.2(inset, (1:nrow(testgep))[-inset], maxCombSize=10^10)
  dbi <- dbs[testi]
  res[[i]] <- list(id=id, testj=testj, ks=ks, dbi=dbi)
}

test_that("KS statistics", {
  i <- 1
  id <- res[[i]]$id; testj <- res[[i]]$testj; ks <- res[[i]]$ks; dbi <- res[[i]]$dbi
  expect_equal(rp$get(dbi)$ES[id, testj], ks$ES)
  expect_equal(rp$get(dbi)$PV[id, testj], ks$p.value)

  i <- 2
  id <- res[[i]]$id; testj <- res[[i]]$testj; ks <- res[[i]]$ks; dbi <- res[[i]]$dbi
  expect_equal(rp$get(dbi)$ES[id, testj], ks$ES)
  expect_equal(rp$get(dbi)$PV[id, testj], ks$p.value)

  i <- 3
  id <- res[[i]]$id; testj <- res[[i]]$testj; ks <- res[[i]]$ks; dbi <- res[[i]]$dbi
  expect_equal(rp$get(dbi)$ES[id, testj], ks$ES)
  expect_equal(rp$get(dbi)$PV[id, testj], ks$p.value)
})



peps1 <- rp$get(expected_dbs[1])
peps3 <- rp$get(expected_dbs[3])

es1 <- peps1$ES
es3 <- peps3$ES
RowRanked1 <- rankPEPsByRows(peps1)
RowRanked3 <- rankPEPsByRows(peps3)

test_that("Row ranking", {
    expect_true(all(apply(RowRanked1, 1, setequal, 1:5)))
    expect_true(all(apply(RowRanked3, 1, setequal, 1:5)))
    expect_equal(RowRanked1[1,1], 5)
    expect_equal(RowRanked1[4,3], 1)
    expect_equal(RowRanked3[5,4], 5)
    expect_equal(RowRanked3[8,3], 1)
})

ColRanked1 <- rankPEPsByCols(peps1)
ColRanked3 <- rankPEPsByCols(peps3)

randj <- sample(ncol(ColRanked3),1)
PVs <- peps1$PV[,randj]
ESs <- peps1$ES[,randj]
if(any(ESs<0)) {
    PVs[ESs>=0] <- Inf
    lastid <- which.min(PVs)
} else lastid <- which.max(PVs)

test_that("Row ranking", {
    expect_true(all(apply(ColRanked1, 2, setequal, 1:10)))
    expect_true(all(apply(ColRanked3, 2, setequal, 1:10)))

    expect_equal(ColRanked1[2,1], 1)
    expect_equal(ColRanked1[1,1], 10)
    expect_equal(ColRanked3[8,3], 1)
    expect_equal(ColRanked3[10,3], 10)
    expect_equal(ColRanked3[4,5], 10)

    expect_equal(ColRanked1[lastid, randj], 10)    
})


pgset <- c("(+)_chelidonine",  "(+/_)_catechin")
res <- PertSEA(rp, pgset)
randi <- sample(1:length(testpws), 1)
pwsid <- testpws[[randi]]$id
randDB <- dbs[randi]
ranked <- rankPEPsByRows(rp$get(randDB))
inset <- ranked[pwsid, pgset]
outset <- ranked[pwsid, setdiff(colnames(ranked), pgset)]
ks <- ks.test.2(inset, outset)

test_that("PertSEA", {
    expect_equal(unname(res[["PertSEA"]][[randDB]][pwsid, "ES"]),
                 ks$ES)
    expect_equal(unname(res[["PertSEA"]][[randDB]][pwsid, "PV"]),
                 ks$p.value)
})



db1 <- expected_dbs[1];
db3 <- expected_dbs[3]
pws1 <- names(rp$get(paste0(db1, "_sets")))[c(2,5,6,9)]
pws3 <- names(rp$get(paste0(db3, "_sets")))[c(1,3,10)]
subpws <- list(pws1, pws3)
names(subpws) <- c(db1, db3)
res <- PathSEA(rp, subpws)

randj1 <- sample(1:ncol(testgep), 1)
ranked <- rankPEPsByCols(rp$get(db1))
peps <- rp$get(db1)
inset <- ranked[pws1, randj1]
outset <- ranked[setdiff(rownames(ranked), pws1), randj1]
ks1 <- ks.test.2(inset, outset)

randj3 <- sample(1:ncol(testgep), 1)
ranked <- rankPEPsByCols(rp$get(db3))
peps <- rp$get(db3)
inset <- ranked[pws3, randj3]
outset <- ranked[setdiff(rownames(ranked), pws3), randj3]
ks3 <- ks.test.2(inset, outset)

name1 <- colnames(testgep)[randj1]
name3 <- colnames(testgep)[randj3]

test_that("PathSEA", {
    expect_equal(unname(res[["PathSEA"]][[db1]][name1, "ES"]),
                 ks1$ES)
    expect_equal(unname(res[["PathSEA"]][[db1]][name1, "PV"]),
                 ks1$p.value)

    expect_equal(unname(res[["PathSEA"]][[db3]][name3, "ES"]),
                 ks3$ES)
    expect_equal(unname(res[["PathSEA"]][[db3]][name3, "PV"]),
                 ks3$p.value)    
})

