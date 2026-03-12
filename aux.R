# ===========================================================
#   Contains several auxiliary functions to GSS
#   Miguel Perez-Enciso
#   mperezenciso@gmail.com
# ===========================================================


#---------------------------
# Calculate MAF for each SNP
#---------------------------
maf = function(X) {
  freq = colMeans(X, na.rm = TRUE) / 2
  freq = 0.5-abs(0.5-freq)
  return(freq)
}


# ============================================================
# Numerator Relationship Matrix via the Tabular Method
# ============================================================
# The tabular method builds A row-by-row using the rules:
#   - a_ii = 1 + 0.5 a[sire_i, dam_i]  (F_i = inbreeding coefficient of i)
#   - a_ij = 0.5 * (a[sire_i, j] + a[dam_i, j])  for j < i
# Inds must be sorted oldest to youngest (ancestors first)
# ============================================================
#-------------------
doA <- function(ped) {
#-------------------
  # pedigree: data.frame with columns [id, sire, dam]
  n <- nrow(ped)
  id <- ped[,1]
  # Initialize A matrix
  A <- diag(1, n)
  sire = ped[,2]
  dam = ped[,3]
  
  for (i in 1:n) {
    for (j in 1:i) {
      if (i == j) {
        # Diagonal: a_ii = 1 + F_i = 1 + 0.5 * a[sire_i, dam_i]
        if (sire[i] != 0 && dam[i] != 0) A[i, i] <- A[i,i] + 0.5 * A[sire[i], dam[i]]
        
      } else {
        # Off-diagonal: a_ij = 0.5 * (a[sire_j, i] + a[dam_j, i])
        if (sire[i] != 0) A[j,i] <- 0.5 * (A[sire[i],j])
        if (dam[i] != 0)  A[j,i] <- A[j,i] + 0.5 * (A[dam[i],j])
        A[i,j] <- A[j,i]  # symmetric
      }
    }
  
  }
  return(A)
}

# ============================================================
# Inverse of Numerator Relationship Matrix via Henderson's rules
# ============================================================
# No inbreeding
# Inds must be sorted oldest to youngest (ancestors first)
# ============================================================
#----------------------
doAInv <- function(ped) {
#----------------------
  # ped is a data frame or matrix with id, father and mother ids, 0 for unknown parent
  w <- c(1, -0.5, -0.5)
  res <- c(2, 4/3, 1)
  nind <- nrow(ped)
  
  # Convert ped to matrix if a data frame
  if (is.data.frame(ped)) { ped <- as.matrix(ped) }
  
  # a dummy 0 row and 0 column is created to facilitate addressing positions
  AI <- matrix(0, nrow = nind + 1, ncol = nind + 1)
  
  s <- integer(nind)
  s[ped[,2] == 0] <- s[ped[,2] == 0] + 1
  s[ped[,3] == 0] <- s[ped[,3] == 0] + 1
  
  # ULL with id indices
  for (i in 1:nind) {
    for (k1 in 1:3) {
      for (k2 in 1:3) {
        AI[ped[i, k1]+1, ped[i, k2]+1] = AI[ped[i, k1]+1, ped[i, k2]+1] + 
                                         w[k1] * w[k2] * res[s[i]+1]
      }
    }
  }
  
  # remove row 0 and col 0
  return(AI[2:(nind + 1), 2:(nind + 1)])
}

# ----------------------------------
#  GRM function (G de los Campos)
# ----------------------------------
GRM=function(X) {
  X=scale(X,center=TRUE,scale=FALSE)
  G=tcrossprod(X)
  G=G/mean(diag(G))
  return(G)
}

# =====================================
# BLUP using G matrix
#  y: phenotypes
#  whichna: ids with missing phenotypes
#  h2: heritability
# RETURN
#  uhat: predicted y
#  pev: predicted error variance
#  rel: reliability
# =====================================
#----------------------------------
BLUP <- function(G, y, whichna, h2) {
#----------------------------------
  n = ncol(G)
  va = h2 * var(as.numeric((y)))
  Ginv = chol2inv(chol(G))
  LHS = diag(n)
  diag(LHS)[whichna] = 0
  LHS = LHS + Ginv * (1-h2)/h2
  LHSinv = (chol2inv(chol(LHS)))
  
  # predictive ability
  y[whichna] = 0
  uhat = (LHSinv %*% y)

  pev = diag(LHSinv)
  rel = 1 - pev / (diag(G)*va)
  return(list(uhat=uhat, pev=pev, rel=rel))
}
