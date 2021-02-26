# 26-Feb-2021
# This R script was written for 
# "Genomic selection for F1 hybrid breeding in strawberry (Fragaria x ananassa)."
# https://doi.org/10.3389/fpls.2021.645111
# This is a source file for functions used in FxaF1GS_Note.R.
# correspondence: yame.repos@gmail.com

makeZ = function(row.id, col.id) {
  Z = matrix(0, nrow=length(row.id), ncol=length(col.id))
  rownames(Z) = row.id
  colnames(Z) = col.id
  Z[cbind(1:nrow(Z),match(row.id,col.id))] = 1
  return(Z)
}#makeZ()
myAmat = function(X, n.core = 1) {
  X = X + 1
  n = nrow(X)
  m = ncol(X)
  v1 = matrix(1, n, 1)
  if (n.core > 1) {
    it = split(1:m, factor(cut(1:m, n.core, labels = FALSE)))
    resit = parallel::mclapply(it, function(markers) { apply(X[, markers], 2, mean) /2 }, mc.cores = n.core)
    p = unlist(resit)
  } else {
    p = apply(X, 2, mean) / 2
  }
  q = 1 - p
  var.A = 2 * mean(p * q)
  Mp = tcrossprod(v1, matrix(p, m, 1))
  W = X - 2 * Mp
  A = tcrossprod(W) / var.A / m
  return(A)
}#myAmat()
myDmat = function(X, n.core = 1) {
  X = X + 1
  n = nrow(X)
  m = ncol(X)
  H = matrix(0, n, m)
  H[X == 1] = 1
  rownames(H) = rownames(X)
  v1 = matrix(1, n, 1)
  if (n.core > 1) {
    it = split(1:m, factor(cut(1:m, n.core, labels = FALSE)))
    resit = parallel::mclapply(it, function(markers) { apply(X[, markers], 2, mean) /2 }, mc.cores = n.core)
    p = unlist(resit)
  } else {
    p = apply(X, 2, mean) / 2
  }
  q = 1 - p
  pq = p * q
  Mpq = tcrossprod(v1, matrix(p * q, m, 1))
  H = H - 2 * Mpq
  D = tcrossprod(H)
  var.D = 2 * mean(pq) * (1 - 2 * mean(pq))
  W = H - 2 * Mpq
  D = tcrossprod(W) / var.D / m
  return(D)
}#myDmat()
simulateHybrid = function(Pair, geno, n.core = 1) 
{
  simHyb = function(Pair) {
    g = matrix(NA, nrow = nrow(geno), ncol = nrow(Pair))
    rownames(g) = geno[,1]
    cn = rep(NA, ncol(g))
    for (i in 1:nrow(Pair)) {
      pair = sort(as.character(unlist(Pair[i, ])))
      gi = geno[, is.element(colnames(geno), pair), drop = FALSE]
      if (ncol(gi) < 2) next
      gi = gi + 1
      gi = apply(gi, 1, sum)
      gi[gi == 0] = -1
      gi[gi == 1] = 0
      gi[gi == 2] = 0
      gi[gi == 3] = 0
      gi[gi == 4] = 1
      g[, i] = gi
      cn[i] = paste(pair, collapse = "x")
    }#i
    colnames(g) = cn
    g = g[, !is.na(colnames(g)), drop = FALSE]
    return(g)
  }#simHyb()
  if (n.core > 1) {
    m = nrow(Pair)
    it = split(1:m, factor(cut(1:m, n.core, labels = FALSE)))
    resit = parallel::mclapply(it, function(ix) { simHyb(Pair[ix, ]) }, mc.cores = n.core)
    N = unlist(lapply(resit, ncol))
    m = sum(N)
    G = matrix(NA, nrow = nrow(geno), ncol = m)
    cn = rep(NA, m)
    for (i in 1:length(N)) {
      sp = min(which(is.na(G[1,])))
      resi = resit[[i]]
      G[, sp:(sp + ncol(resi) - 1)] = resi
      cn[sp:(sp + ncol(resi) - 1)] = colnames(resi)
    }#i
    colnames(G) = cn
  } else {
    G = simHyb(Pair)
  }
  return(G)
}#simulateHybrid()
