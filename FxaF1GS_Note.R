# 22-Dec-2020
# This R script was written for 
# "Genomic selection for F1 hybrid breeding in strawberry (Fragaria x ananassa)."
# correspondence: yame.repos@gmail.com

################################################
###    Genetic Relationship (Figure 2C)      ###
#############################################

n.core = 20

# loading data and functions
load("FxaF1GS.RData")
attach(Data)
source("FxaF1GS.R")

# simulating 5,460 possible hybrid genome
id = colnames(geno)[-(1:3)]
ix = 1:length(id)
pair = expand.grid(id, id)
ix = expand.grid(ix, ix)
pair = pair[which(ix$Var1 < ix$Var2),]
ph = simulateHybrid(pair, geno, n.core = n.core)
ix.ph = rep(1, ncol(ph))
rm(id, ix) ; gc()

# simulating 275 test F1 hybrid genome
pair = pheno[is.element(pheno$type, "trainF1"), is.element(colnames(pheno), c("Seed", "Pollen"))]
th = simulateHybrid(pair, geno, n.core = n.core)
ix.th = rep(2, ncol(th))
rm(pair) ; gc()

# integrating genotype information
X = t(cbind(ph, th, geno[,-(1:3)]))
col.ix = pch.ix = c(ix.ph, ix.th, rep(3, ncol(geno) - 3))
rm(ph, th, ix.ph, ix.th) ; gc()

# calculating additive relationship matrix
A = myAmat(X, n.core = n.core)

# performing PCA
pca = prcomp(A)

# plotting PCA results
col.ix[col.ix == 1] = "grey"
col.ix[col.ix == 2] = "black"
col.ix[col.ix == 3] = "black"
pch.ix[pch.ix == 3] = 21
pch.ix[pch.ix == 2] = 3
pch.ix[pch.ix == 1] = 3
plot(pca$x[, 1], pca$x[, 2], col = col.ix, pch = pch.ix, bg = "white")


################################################
###    Heritability in 105-ILs (Table 1)     ###
#############################################

n.core = 20

# loading data and functions
load("FxaF1GS.RData")
attach(Data)
source("FxaF1GS.R")

# associating genotype to phenotype
pheno = pheno[is.element(pheno$type, "trainP"), -c(2, 3, 4)]
use.id = intersect(colnames(geno), pheno[,1])
geno = cbind(geno[, 1:3], geno[, is.element(colnames(geno), use.id)])
pheno = pheno[is.element(pheno[, 1], use.id), ]
X = t(geno[, -(1:3)])
Z = makeZ(pheno[, 1], rownames(X))
X = Z %*% X
rm(use.id, Z) ; gc()

# calculating additive and dominant relationship matrices
A = myAmat(X, n.core = n.core)
D = myDmat(X, n.core = n.core)

# estimating heritability
res = matrix(NA, nrow = ncol(pheno) - 1, ncol = 2)
colnames(res) = c("A", "AD")
rownames(res) = colnames(pheno)[-1]
for (i in 2:ncol(pheno)) {
  y = pheno[, i]
  use = !is.na(y)
  Ai = A[use, use]
  Di = D[use, use]
  soln = BGLR::BGLR(y = y[use], 
                    ETA = list(list(K=Ai, model='RKHS')), 
                    nIter = 15000, burnIn = 5000, 
                    verbose = FALSE, saveAt = "A")
  varU = scan('AETA_1_varU.dat')
  varE = scan('AvarE.dat')
  h2a = varU / (varU + varE)
  soln = BGLR::BGLR(y = y[use], 
                    ETA = list(list(K = Ai, model='RKHS'), list(K = Di, model='RKHS')), 
                    nIter = 15000, burnIn = 5000, 
                    verbose = FALSE, saveAt = "AD")
  varU1 = scan('ADETA_1_varU.dat')
  varU2 = scan('ADETA_2_varU.dat')
  varE = scan('ADvarE.dat')
  h2ad = (varU1 + varU2) / (varU1 + varU2 + varE)
  res[i-1, ] = c(mean(h2a), mean(h2ad))
  system("rm *varU.dat")
  system("rm *varE.dat")
  system("rm *mu.dat")
}#i
rm(h2a, h2ad, i, use, varU, varU1, varU2, varE, y, X, soln, Ai, Di, A, D) ; gc()
print(res)
write.csv(res, "105-ILs_h2.csv")


################################################
###    Heritability in 275-F1s (Table 1)     ###
#############################################

n.core = 20

# loading data and functions
load("FxaF1GS.RData")
attach(Data)
source("FxaF1GS.R")

# simulating 275 test F1 hybrid genome
pair = pheno[is.element(pheno$type, "trainF1"), is.element(colnames(pheno), c("Seed", "Pollen"))]
X = t(simulateHybrid(pair, geno, n.core = n.core))
rm(pair) ; gc()

# associating genotype to phenotype
pheno = pheno[is.element(pheno$type, "trainF1"), -c(2, 3, 4)]
use.id = intersect(rownames(X), pheno[,1])
pheno = pheno[is.element(pheno[, 1], use.id), ]
X = X[is.element(rownames(X), use.id), ]
Z = makeZ(pheno[, 1], rownames(X))
X = Z %*% X

# calculating additive and dominant relationship matrices
A = myAmat(X, n.core = n.core)
D = myDmat(X, n.core = n.core)
rm(use.id, Z) ; gc()

# estimating heritability
res = matrix(NA, nrow = ncol(pheno) - 1, ncol = 2)
colnames(res) = c("A", "AD")
rownames(res) = colnames(pheno)[-1]
for (i in 2:ncol(pheno)) {
  y = pheno[, i]
  use = !is.na(y)
  Ai = A[use, use]
  Di = D[use, use]
  soln = BGLR::BGLR(y = y[use], 
                    ETA = list(list(K=Ai, model='RKHS')), 
                    nIter = 15000, burnIn = 5000, 
                    verbose = FALSE, saveAt = "A")
  varU = scan('AETA_1_varU.dat')
  varE = scan('AvarE.dat')
  h2a = varU / (varU + varE)
  soln = BGLR::BGLR(y = y[use], 
                    ETA = list(list(K = Ai, model='RKHS'), list(K = Di, model='RKHS')), 
                    nIter = 15000, burnIn = 5000, 
                    verbose = FALSE, saveAt = "AD")
  varU1 = scan('ADETA_1_varU.dat')
  varU2 = scan('ADETA_2_varU.dat')
  varE = scan('ADvarE.dat')
  h2ad = (varU1 + varU2) / (varU1 + varU2 + varE)
  res[i-1, ] = c(mean(h2a), mean(h2ad))
  system("rm *varU.dat")
  system("rm *varE.dat")
  system("rm *mu.dat")
}#i
rm(h2a, h2ad, i, use, varU, varU1, varU2, varE, y, X, soln, Ai, Di, A, D) ; gc()
print(res)
write.csv(res, "275-F1s_h2.csv")


################################################
###  Cross-validation in 105-ILs (Table 2)   ###
#############################################

n.core = 20

# loading data and functions
load("FxaF1GS.RData")
attach(Data)
source("FxaF1GS.R")

# cross-validation
pheno = pheno[is.element(pheno$type, "trainP"), -(2:4)]
traits = colnames(pheno)[-1]
for (trait in traits) {
  y = pheno[[trait]]
  names(y) = pheno[,1]
  
  n.fold = 2
  n.cycles = 100
  
  y = y[!is.na(y)]
  use = intersect(names(y), colnames(geno))
  geno = cbind(geno[, 1:3], geno[, is.element(colnames(geno), use)])
  y = y[is.element(names(y), use)]
  
  Z = makeZ(names(y), colnames(geno)[-(1:3)])
  
  X = t(geno[, -(1:3)])
  H = matrix(0, nrow = nrow(X), ncol = ncol(X))
  H[X == 0] = 1
  XH = cbind(X, H)
  
  X = Z %*% X
  XH = Z %*% XH
  
  A = myAmat(X, n.core = n.core)
  D = myDmat(X, n.core = n.core)
  G = dist(X, method = "euclidean")
  
  A = Z %*% A %*% t(Z)
  D = Z %*% D %*% t(Z)
  
  hp.BL = matrix(c(1, 0.001, 1, 0.01, 1, 0.1, 1, 1, 1, 5), nc = 2, byrow = TRUE)
  hp.BB = matrix(c(5, 1, 0.001, 5, 1, 0.01, 5, 1, 0.1, 5, 1, 0.5), nc = 3, byrow = TRUE)
  
  dat = data.frame(pheno = y, geno = rownames(X))
  
  n = length(y)
  
  Res = matrix(NA, nrow = n.cycles, ncol = 8)
  colnames(Res) = c("A", "AD", "BB", "BBd", "BL", "BLd", "G", "RF")
  
  for (cyc in 1:n.cycles) {
    set.seed(cyc)
    
    ids = dat$geno
    size.matrix = ceiling(n / n.fold) * n.fold
    add.NAs = rep(NA, size.matrix - n)
    shuffle = sample(1:n)
    Partition = matrix(c(shuffle, add.NAs), ncol=n.fold, byrow=TRUE)
    
    res = data.frame(resA = numeric(n.fold),
                     resAD = numeric(n.fold),
                     resBB = numeric(n.fold),
                     resBBd = numeric(n.fold),
                     resBL = numeric(n.fold),
                     resBLd = numeric(n.fold),
                     resG = numeric(n.fold),
                     resRF = numeric(n.fold))
    for (i in 1:n.fold) {
      pred = Partition[, i]
      pred = pred[!is.na(pred)]
      dati = dat
      dati$pheno[pred] = NA
      soln = suppressWarnings(BGLR::BGLR(y = dati$pheno, 
                                         ETA = list(list(K = A, model='RKHS')), 
                                         burnIn = 5000, nIter = 15000, 
                                         verbose = FALSE, saveAt = paste0(trait, cyc)))
      res$resA[i] = cor(soln$yHat[pred], dat$pheno[pred])
      soln = suppressWarnings(BGLR::BGLR(y = dati$pheno, 
                                         ETA = list(list(K = A, model='RKHS'), list(K = D, model='RKHS')), 
                                         burnIn = 5000, nIter = 15000, 
                                         verbose = FALSE, saveAt = paste0(trait, cyc)))
      res$resAD[i] = cor(soln$yHat[pred], dat$pheno[pred])
      
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.BB, Function = "tuning")
      MSE = soln$MSE
      hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:4, drop = FALSE])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.USE, Function = "fitting")
      yhat = X %*% soln$Beta
      res$resBB[i] = cor(yhat[pred], dat$pheno[pred])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.BB, Function = "tuning")
      MSE = soln$MSE
      hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:4, drop = FALSE])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = XH, Method = "BayesB", Hyperparameters = hp.USE, Function = "fitting")
      yhat = XH %*% soln$Beta
      res$resBBd[i] = cor(yhat[pred], dat$pheno[pred])
      
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.BL, Function = "tuning")
      MSE = soln$MSE
      hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:3, drop = FALSE])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.USE, Function = "fitting")
      yhat = X %*% soln$Beta
      res$resBL[i] = cor(yhat[pred], dat$pheno[pred])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.BL, Function = "tuning")
      MSE = soln$MSE
      hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:3, drop = FALSE])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = XH, Method = "BL", Hyperparameters = hp.USE, Function = "fitting")
      yhat = XH %*% soln$Beta
      res$resBLd[i] = cor(yhat[pred], dat$pheno[pred])
      
      soln = rrBLUP::kin.blup(data = dati, geno = "geno", pheno = "pheno", K = G, GAUSS = TRUE)
      res$resG[i] = cor(soln$pred[pred], dat$pheno[pred])
      
      rf = randomForest::randomForest(x = X[-pred,], y = dat$pheno[-pred])
      yhat = predict(rf, X[pred,])
      res$resRF[i] = cor(yhat, dat$pheno[pred])
    }#i
    Res[cyc, ] = apply(res, 2, mean, na.rm = TRUE)
  }#cyc
  csv.name = paste0("105-ILs_", trait, "_CV.cs")
  write.csv(res, csv.name)
}#trait


################################################
###  Cross-validation in 275-F1s (Table 2)   ###
#############################################

n.core = 20

# loading data and functions
load("FxaF1GS.RData")
attach(Data)
source("FxaF1GS.R")

# associating genotype to phenotype
pair = pheno[is.element(pheno$type, "trainF1"), is.element(colnames(pheno), c("Seed", "Pollen"))]
X = simulateHybrid(pair, geno, n.core = n.core)
geno = cbind(geno[, 1:3], X)

# cross-validation
pheno = pheno[is.element(pheno$type, "trainF1"), -(2:4)]
traits = colnames(pheno)[-1]
for (trait in traits) {
  y = pheno[[trait]]
  names(y) = pheno[,1]
  
  n.fold = 2
  n.cycles = 100
  
  y = y[!is.na(y)]
  use = intersect(names(y), colnames(geno))
  geno = cbind(geno[, 1:3], geno[, is.element(colnames(geno), use)])
  y = y[is.element(names(y), use)]
  
  Z = makeZ(names(y), colnames(geno)[-(1:3)])
  
  X = t(geno[, -(1:3)])
  H = matrix(0, nrow = nrow(X), ncol = ncol(X))
  H[X == 0] = 1
  XH = cbind(X, H)
  
  X = Z %*% X
  XH = Z %*% XH
  
  A = myAmat(X, n.core = n.core)
  D = myDmat(X, n.core = n.core)
  G = dist(X, method = "euclidean")
  
  A = Z %*% A %*% t(Z)
  D = Z %*% D %*% t(Z)
  
  hp.BL = matrix(c(1, 0.001, 1, 0.01, 1, 0.1, 1, 1, 1, 5), nc = 2, byrow = TRUE)
  hp.BB = matrix(c(5, 1, 0.001, 5, 1, 0.01, 5, 1, 0.1, 5, 1, 0.5), nc = 3, byrow = TRUE)
  
  dat = data.frame(pheno = y, geno = rownames(X))
  
  n = length(y)
  
  Res = matrix(NA, nrow = n.cycles, ncol = 8)
  colnames(Res) = c("A", "AD", "BB", "BBd", "BL", "BLd", "G", "RF")
  
  for (cyc in 1:n.cycles) {
    set.seed(cyc)
    
    ids = dat$geno
    size.matrix = ceiling(n / n.fold) * n.fold
    add.NAs = rep(NA, size.matrix - n)
    shuffle = sample(1:n)
    Partition = matrix(c(shuffle, add.NAs), ncol=n.fold, byrow=TRUE)
    
    res = data.frame(resA = numeric(n.fold),
                     resAD = numeric(n.fold),
                     resBB = numeric(n.fold),
                     resBBd = numeric(n.fold),
                     resBL = numeric(n.fold),
                     resBLd = numeric(n.fold),
                     resG = numeric(n.fold),
                     resRF = numeric(n.fold))
    for (i in 1:n.fold) {
      pred = Partition[, i]
      pred = pred[!is.na(pred)]
      dati = dat
      dati$pheno[pred] = NA
      soln = suppressWarnings(BGLR::BGLR(y = dati$pheno, 
                                         ETA = list(list(K = A, model='RKHS')), 
                                         burnIn = 5000, nIter = 15000, 
                                         verbose = FALSE, saveAt = paste0(trait, cyc)))
      res$resA[i] = cor(soln$yHat[pred], dat$pheno[pred])
      soln = suppressWarnings(BGLR::BGLR(y = dati$pheno, 
                                         ETA = list(list(K = A, model='RKHS'), list(K = D, model='RKHS')), 
                                         burnIn = 5000, nIter = 15000, 
                                         verbose = FALSE, saveAt = paste0(trait, cyc)))
      res$resAD[i] = cor(soln$yHat[pred], dat$pheno[pred])
      
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.BB, Function = "tuning")
      MSE = soln$MSE
      hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:4, drop = FALSE])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.USE, Function = "fitting")
      yhat = X %*% soln$Beta
      res$resBB[i] = cor(yhat[pred], dat$pheno[pred])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.BB, Function = "tuning")
      MSE = soln$MSE
      hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:4, drop = FALSE])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = XH, Method = "BayesB", Hyperparameters = hp.USE, Function = "fitting")
      yhat = XH %*% soln$Beta
      res$resBBd[i] = cor(yhat[pred], dat$pheno[pred])
      
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.BL, Function = "tuning")
      MSE = soln$MSE
      hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:3, drop = FALSE])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.USE, Function = "fitting")
      yhat = X %*% soln$Beta
      res$resBL[i] = cor(yhat[pred], dat$pheno[pred])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.BL, Function = "tuning")
      MSE = soln$MSE
      hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:3, drop = FALSE])
      soln = VIGoR::vigor(Pheno = dati$pheno, Geno = XH, Method = "BL", Hyperparameters = hp.USE, Function = "fitting")
      yhat = XH %*% soln$Beta
      res$resBLd[i] = cor(yhat[pred], dat$pheno[pred])
      
      soln = rrBLUP::kin.blup(data = dati, geno = "geno", pheno = "pheno", K = G, GAUSS = TRUE)
      res$resG[i] = cor(soln$pred[pred], dat$pheno[pred])
      
      rf = randomForest::randomForest(x = X[-pred,], y = dat$pheno[-pred])
      yhat = predict(rf, X[pred,])
      res$resRF[i] = cor(yhat, dat$pheno[pred])
    }#i
    Res[cyc, ] = apply(res, 2, mean, na.rm = TRUE)
  }#cyc
  csv.name = paste0("275-F1s_", trait, "_CV.cs")
  write.csv(res, csv.name)
}#trait



################################################
###  Across-population prediction (Table 3)  ###
#############################################

n.core = 20

# loading data and functions
load("FxaF1GS.RData")
attach(Data)
source("FxaF1GS.R")

# simulating 275 test F1 hybrid genome
pair = pheno[is.element(pheno$type, "trainF1"), is.element(colnames(pheno), c("Seed", "Pollen"))]
th = simulateHybrid(pair, geno, n.core = n.core)
geno = cbind(geno, th)
rm(pair) ; gc()

# preparing formats
pheno = pheno[is.element(pheno$type, c("trainP", "trainF1")), -(2:3)]
use = intersect(pheno[,1], colnames(geno))
geno = cbind(geno[, 1:3], geno[, is.element(colnames(geno), use)])
pheno = pheno[is.element(pheno[,1], use), ]
rm(use) ; gc()

X = t(geno[, -(1:3)])
H = matrix(0, nrow = nrow(X), ncol = ncol(X))
H[X == 0] = 1
XH = cbind(X, H)
rm(H) ; gc()

Z = makeZ(pheno[,1], rownames(X))
X = Z %*% X
XH = Z %*% XH
rm(Z) ; gc()

A = myAmat(X, n.core = n.core)
D = myDmat(X, n.core = n.core)
G = dist(X, method = "euclidean")

traits = colnames(pheno)[-(1:2)]

hp.BL = matrix(c(1, 0.001, 1, 0.01, 1, 0.1, 1, 1, 1, 5), nc = 2, byrow = TRUE)
hp.BB = matrix(c(5, 1, 0.001, 5, 1, 0.01, 5, 1, 0.1, 5, 1, 0.5), nc = 3, byrow = TRUE)

# prediction
for (p in 1:5) {
  trait = traits[p]
  print(trait)
  
  dat = data.frame(pheno = pheno[[trait]], geno = pheno[, 1])
  
  res = data.frame(resA = numeric(2),
                   resAD = numeric(2),
                   resBB = numeric(2),
                   resBBd = numeric(2),
                   resBL = numeric(2),
                   resBLd = numeric(2),
                   resG = numeric(2),
                   resRF = numeric(2))
  for (i in 1:2) {
    pred = NULL
    if (i == 1) pred = is.element(pheno$type, "trainF1")
    if (i == 2) pred = is.element(pheno$type, "trainP") 
    dati = dat
    dati$pheno[pred] = NA
    soln = suppressWarnings(BGLR::BGLR(y = dati$pheno, 
                                       ETA = list(list(K = A, model='RKHS')), 
                                       burnIn = 5000, nIter = 15000, 
                                       verbose = FALSE, saveAt = paste0(trait, "_cv_")))
    #system("rm *.dat")
    res$resA[i] = cor(soln$yHat[pred], dat$pheno[pred], use = "pairwise.complete.obs")
    soln = suppressWarnings(BGLR::BGLR(y = dati$pheno, 
                                       ETA = list(list(K = A, model='RKHS'), list(K = D, model='RKHS')), 
                                       burnIn = 5000, nIter = 15000, 
                                       verbose = FALSE, saveAt = paste0(trait, "_cv_")))
    #system("rm *.dat")
    res$resAD[i] = cor(soln$yHat[pred], dat$pheno[pred], use = "pairwise.complete.obs")
    
    soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.BB, Function = "tuning")
    MSE = soln$MSE
    hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:4, drop = FALSE])
    soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.USE, Function = "fitting")
    yhat = X %*% soln$Beta
    res$resBB[i] = cor(yhat[pred], dat$pheno[pred], use = "pairwise.complete.obs")
    soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BayesB", Hyperparameters = hp.BB, Function = "tuning")
    MSE = soln$MSE
    hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:4, drop = FALSE])
    soln = VIGoR::vigor(Pheno = dati$pheno, Geno = XH, Method = "BayesB", Hyperparameters = hp.USE, Function = "fitting")
    yhat = XH %*% soln$Beta
    res$resBBd[i] = cor(yhat[pred], dat$pheno[pred], use = "pairwise.complete.obs")
    
    soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.BL, Function = "tuning")
    MSE = soln$MSE
    hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:3, drop = FALSE])
    soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.USE, Function = "fitting")
    yhat = X %*% soln$Beta
    res$resBL[i] = cor(yhat[pred], dat$pheno[pred], use = "pairwise.complete.obs")
    soln = VIGoR::vigor(Pheno = dati$pheno, Geno = X, Method = "BL", Hyperparameters = hp.BL, Function = "tuning")
    MSE = soln$MSE
    hp.USE = as.matrix(MSE[MSE$MSE == min(MSE$MSE), 2:3, drop = FALSE])
    soln = VIGoR::vigor(Pheno = dati$pheno, Geno = XH, Method = "BL", Hyperparameters = hp.USE, Function = "fitting")
    yhat = XH %*% soln$Beta
    res$resBLd[i] = cor(yhat[pred], dat$pheno[pred], use = "pairwise.complete.obs")
    
    soln = rrBLUP::kin.blup(data = dati, geno = "geno", pheno = "pheno", K = G, GAUSS = TRUE)
    res$resG[i] = cor(soln$pred[pred], dat$pheno[pred], use = "pairwise.complete.obs")
    
    rf = randomForest::randomForest(x = X[!is.na(dati$pheno),], y = dat$pheno[!is.na(dati$pheno)])
    yhat = predict(rf, X[pred,])
    res$resRF[i] = cor(yhat, dat$pheno[pred], use = "pairwise.complete.obs")
  }#i
  rownames(res) = c("ILtoF1", "F1toIL")
  csv.name = paste0("AcrossPop_", trait, ".csv")
  write.csv(res, csv.name)
}#p


################################################
###       Genomic selection (Figure 4)       ###
#############################################

n.core = 20

# loading data and functions
load("FxaF1GS.RData")
attach(Data)
source("FxaF1GS.R")

# simulating GS F1 hybrids
pair = pheno[is.element(pheno$type, "testF1"), is.element(colnames(pheno), c("Seed", "Pollen"))]
gsh = simulateHybrid(pair, geno, n.core = 1)
colnames(gsh) = paste0(colnames(gsh), ".test")

# simulating 275 test F1 hybrid genome
pair = pheno[is.element(pheno$type, "trainF1"), is.element(colnames(pheno), c("Seed", "Pollen"))]
th = simulateHybrid(pair, geno, n.core = n.core)

# simulating 5,460 possible hybrid genome
id = colnames(geno)[-(1:3)]
ix = 1:length(id)
pair = expand.grid(id, id)
ix = expand.grid(ix, ix)
pair = pair[which(ix$Var1 < ix$Var2),]
ph = simulateHybrid(pair, geno, n.core = n.core)
ph = ph[, !is.element(colnames(ph), colnames(th))]
rm(id, ix) ; gc()

# arranging data format
geno = cbind(geno, gsh, th, ph)
rm(pair, gsh, th, ph) ; gc()

pheno = pheno[!is.element(pheno$type, "testF1"), -(2:3)]
pheno = pheno[is.element(pheno[, 1], colnames(geno)), ]

addIDs = colnames(geno)[-(1:3)]
addIDs = addIDs[!is.element(addIDs, pheno[, 1])]
addMat = as.data.frame(matrix(NA, nrow = length(addIDs), ncol = ncol(pheno)))
colnames(addMat) = colnames(pheno)
addMat$ID = addIDs
pheno = rbind(pheno, addMat)
rm(addMat, addIDs) ; gc()

# associating genotype to phenotype
X = t(geno[, -(1:3)])
pix = pheno
pix[] = NA
for (i in 1:nrow(X)) pix[i,] = pheno[pheno$ID == rownames(X)[i],]
pheno = pix
pheno$type[is.na(pheno$type)] = "trainF1"
rm(pix) ; gc()

# calculating additive and dominant relationship matrices
A = myAmat(X, n.core = n.core)
D = myDmat(X, n.core = n.core)

# constructing GS models and predicting phenotypes
traits = colnames(pheno)[-(1:2)]
for (p in 1:5) {
  trait = traits[p]
  print(trait)
  soln = suppressWarnings(BGLR::BGLR(y = pheno[[trait]], 
                                     ETA = list(list(K = A, model = 'RKHS'), 
                                                list(K = D, model = 'RKHS'), 
                                                list(~ factor(type), 
                                                     data = pheno, model = 'FIXED')), 
                                     burnIn = 5000, nIter = 15000, 
                                     verbose = FALSE, saveAt = paste0(trait, "_gs_")))
  gs = data.frame(gid = rownames(X),
                  y = soln$y,
                  yHat = soln$yHat)
  save(gs, file = paste0(trait, "_gs.RData"))
}#p
system("rm *.dat")
rm(A, D, soln, X, i) ; gc()

# arranging results
files = dir()
files = files[grep("_gs.RData", files)]
cn = c()
ct = 1
for (fi in files) {
  load(fi)
  if (ct == 1) {
    Yhat = matrix(NA, nrow = nrow(gs), ncol = length(files))
    rownames(Yhat) = gs$gid
  }
  cni = unlist(strsplit(fi, "_gs.RData"))
  cn = c(cn, cni)
  Yhat[, ct] = gs$yHat
  ct = ct + 1
}#fi
colnames(Yhat) = cn
Yhat = as.data.frame(Yhat)
rm(gs, cn, cni, ct, fi, files) ; gc()

# plotting results
gid = rownames(Yhat)
col.ix = rep("grey", nrow(Yhat))
col.ix[grep(".test", gid)] = "black"
pch.ix = rep(3, nrow(Yhat))
ix = rep("o", nrow(Yhat))
dat = data.frame(gid = gid,
                 col = col.ix,
                 pch = pch.ix,
                 ix = ix, 
                 FruitHardness = Yhat$FruitHardness,
                 PericarpColor = Yhat$PericarpColor,
                 PetioleLength = Yhat$PetioleLength, 
                 Brix = Yhat$Brix, 
                 LeafArea = Yhat$LeafArea)
dat = dat[rev(order(dat$col)),]

LFH =c("A42-04xC46-04.test", "A22-01xB08-01.test", "B08-01xC42-02.test", "B12-06xC42-02.test", "B30-09xC42-02.test", "B04-12xC47-11.test")
HFH = c("A48-01xA51-10.test", "A48-01xA51-02.test")
HPC = c("B26-09xC50-03.test", "B05-11xC50-03.test", "B16-02xC50-03.test", "B12-06xC50-01.test")
LPC = c("A49-08xB24-12.test", "A51-02xB24-05.test", "A53-08xB24-05.test", "A54-03xB24-12.test", "A54-03xA60-18.test", "A48-01xA64-03.test", "A55-08xB23-08.test")
IM = c("A51-02xC50-03.test", "A48-01xC44-08.test")

dat$pch[is.element(dat$gid, LFH)] = 23
dat$pch[is.element(dat$gid, HPC)] = 19
dat$pch[is.element(dat$gid, HFH)] = 18
dat$pch[is.element(dat$gid, LPC)] = 21
dat$pch[is.element(dat$gid, IM)] = 17

dat$ix[dat$pch == 23] = "Low-FH"
dat$ix[dat$pch == 19] = "High-PC"
dat$ix[dat$pch == 18] = "Hihg-FH"
dat$ix[dat$pch == 21] = "Low-PC"
dat$ix[dat$pch == 17] = "IM"

# Figure 4A
plot(dat$FruitHardness, dat$PericarpColor, col = dat$col, pch = dat$pch, 
     bg = "white", xlab = "Fruit hardness", ylab = "Pericarp color",
     main = "Figure 4A")

# Figure 4B
obs = Data$pheno[Data$pheno$type == "testF1", ]
pre = dat[!is.element(dat$ix, "o"), ]
use = intersect(pre$gid, obs$ID)
obs = obs[is.element(obs$ID, use), ]
pre = pre[is.element(pre$gid, use), ]
obs = obs[sort(order(obs$ID)), ]
pix = pre
for (i in 1:nrow(pix)) pix[i, ] = pre[pre$gid == obs$ID[i], ]
plot(obs$FruitHardness, obs$PericarpColor, pch = pix$pch, col = pix$col,
     bg = "white", xlab = "Fruit hardness", ylab = "Pericarp color", main = "Figure 4B")

# Figure 4C
plot(pix$PetioleLength, obs$PetioleLength, pch = pix$pch, col = pix$col,
     bg = "white", xlab = "Predicted values", ylab = "Observed values", main = "Petiole length")
abline(lm(obs$PetioleLength ~ pix$PetioleLength))

plot(pix$LeafArea, obs$LeafArea, pch = pix$pch, col = pix$col,
     bg = "white", xlab = "Predicted values", ylab = "Observed values",main = "Leaf area")
abline(lm(obs$LeafArea ~ pix$LeafArea))

plot(pix$Brix, obs$Brix, pch = pix$pch, col = pix$col,
     bg = "white", xlab = "Predicted values", ylab = "Observed values", main = "Brix")
abline(lm(obs$Brix ~ pix$Brix))

plot(pix$FruitHardness, obs$FruitHardness, pch = pix$pch, col = pix$col,
     bg = "white", xlab = "Predicted values", ylab = "Observed values", main = "Fruit hardness")
abline(lm(obs$FruitHardness ~ pix$FruitHardness))

plot(pix$PericarpColor, obs$PericarpColor, pch = pix$pch, col = pix$col,
     bg = "white", xlab = "Predicted values", ylab = "Observed values", main = "Pericarp color")
abline(lm(obs$PericarpColor ~ pix$PericarpColor))


