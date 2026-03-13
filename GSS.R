# Load AlphaSimR package

library(AlphaSimR)
library(MASS)
library(ggplot2)

setwd('/Users/miguel/Dropbox/0_Research/2026_GSStep/src')

# load auxiliary functions
source('aux.R')

# ============================================================
# alphasimR wrapper
# ============================================================
# simulates nGen cycles of nGen discrete generations with 
# random crossings 
#  - h2: heritability
#  - nFounder: n founders
#  - nOffspring: per cross
#  - nChr: no chromosomes
#  - nQTL per chr
#  - nSNP per chr
# ============================================================
#------------------------------------------------------------------------
runAlphasim <- function(h2, nFounder, nOffspring, nChr, nQTL, nSNP, nGen) {
#------------------------------------------------------------------------
   #--> Founder population
   founderPop = runMacs(nInd = nFounder,  
                        nChr = nChr,   
                        segSites = nSNP)
   
   # ULL: SP is now global variable
   SP <<- SimParam$new(founderPop)
   SP$addTraitA(nQtlPerChr = nQTL)   # add QTL
   SP$setVarE(h2 = h2)               # sets error variance to achieve desired h²

   #--> mating cycles
   pops = list()
   pops[[1]] = newPop(founderPop) 
   for (t in 1:nGen) {
      # Random mating 
      pops[[t+1]] <- randCross(pops[[t]], 
                               nCrosses = nFounder,
                               nProgeny = nOffspring,
                               simParam = SP)
   }
   # merge all generations
   pop = Reduce(c, pops)

   #--> Assign phenotypes 
   pop = setPheno(pop, h2 = h2, simParam = SP)
   y = pop@pheno

   #--> Get all segregating SNPs + QTL
   g = pullSegSiteGeno(pop, simParam = SP)
   
   #--> Get all segregating SNPs + QTL
   q = pullQtlGeno(pop, simParam = SP) 
   
   #--> Get QTL genotypes
   tbv = pop@gv

   #--> get ped, add generation column
   ped = data.frame(pop@iid, as.integer(pop@father), as.integer(pop@mother))
   ped <- cbind( ped, c(rep(0,nFounder), sort(rep(1:nGen,nFounder*nOffspring))) )

   # returns pedigree, phenotypes, sequence, qtl genotypes, true breeding value
   return(list(ped=ped, y=y, g=g, q=q, tbv=tbv))
}


# ======================================================
# Computes H matrix for any number of hierarchical GRMs
# ======================================================
#------------------------
doHgss <- function(Glist) {
#------------------------
  # Input is a list of molecular relationship matrices in increasing marker density
  # Output is the covariance H matrix
  # WARNING: G list is in decreasing size order
  #          G[[1]] supposed to be NRM but could be the lowest density marker GRM
  #                 if pedigree unavailable
  # All individuals in G[i+1] are to be in G[i] and in same order
  # All markers in G[i] are supposed to be in G[i+1]
  nb = length(Glist)     # n G's
  n  = ncol(Glist[[1]])  # ninds
  H  = matrix(0, nrow=n, ncol=n)
  
  n0 = ncol(Glist[[nb]])
  # This is the span of inds for highest density genotyped
  H[(n-n0+1):n, (n-n0+1):n] = Glist[[nb]]
  for (i in (nb-1):1) {
    G  = Glist[[i]]   # G of current block
    n1 = ncol(G)      # coordinates wrt to current G
    b1 = 1:(n1-n0)    # individuals with current array ungenotyped with former denser array
    b0 = (n1-n0+1):n1 # previous individuals genotyped with denser array
    h1 = b1 + (n-n1)  # coordinates in global H: 1...nind
    h0 = b0 + (n-n1)
    Gi = chol2inv(chol(G[b0,b0]))
    H[h1,h1] = G[b1,b1] - 
               t(G[b0,b1]) %*% Gi %*% G[b0,b1] + 
               t(G[b0,b1]) %*% Gi %*% H[h0,h0] %*% Gi %*% G[b0,b1]
    H[h1,h0] = t(G[b0,b1]) %*% Gi %*% H[h0,h0]
    H[h0,h1] = t(H[h1,h0])
    n0 = n1
  }
  return(H)
}

#-----------------------
doHinv = function(Glist) {
#-----------------------
  # return H inverse using algorithm in Legarra et al.
  # input is list of G matrices in increasing marker density
  nb = length(Glist)
  n  = ncol(Glist[[1]]) # total nind
  n1 = ncol(Glist[[2]])  # nind of lowest density array
  b0 = (n-n1+1):n  # inds genotyped 
  
  A = Glist[[1]] # NRM
  
  # Hinv is initialized to NRM inverse (brute force simply for demonstration purposes)
  Hinv = chol2inv(chol(A))
  
  # Only H of matrices except last one is needed
  H = doHgss(Glist[-1])
  H = chol2inv(chol(H)) - chol2inv(chol(A[b0,b0]))
  Hinv[b0,b0] = Hinv[b0,b0] + H
  return(Hinv)
}

# ============================================================
# script begins here
# ============================================================

# for reproducibility
set.seed(42)

# The actual script begins here
# Parameters
h2 = 0.5     # heritability
nFounder = 1000  # n founder individuals
nOff = 1     # n offspring / cross
nChr = 20    # n chromosomes
nQTL = 5     # nQTL / chr
nSNP = 1000  # num SNP / chr
nGen = 3     # no. of generations

h2 = 0.5     # heritability
nFounder = 100  # n founder individuals
nOff = 1     # n offspring / cross
nChr = 20    # n chromosomes
nQTL = 5     # nQTL / chr
nSNP = 1000  # num SNP / chr
nGen = 3     # no. of generations

# Simulation program 
out <- runAlphasim(h2, nFounder=nFounder, nOffspring=nOff, nChr = nChr, nQTL=nQTL, nSNP = nSNP, nGen = nGen)
g <- out$g  # sequence data (whole population)
y <- out$y # phenotype
tbv <- out$tbv # true breeding value
ped <- out$ped # pedigree

nmkr <- ncol(g)
nind <- nrow(g)

hist(maf(g), main = 'Allele frequency distribution')
hist(colMeans(g/2), main = 'Unfolded site frequency spectrum')
heatmap(cor(g[,1:100]), Rowv=NA, Colv=NA, main = 'LD map, first 100 markers')
plot(out$tbv, out$y, xlab='TBV', ylab='y')

library(reshape2)
# Plot 1: Unfolded site frequency spectrum
df1 <- data.frame(allele_freq = colMeans(g/2))

p1 <- ggplot(df1, aes(x = allele_freq)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Unfolded Site Frequency Spectrum",
       x = "Allele Frequency", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot 2: LD map for first 100 markers
# Calculate correlation matrix for first 100 markers
cor_matrix <- cor(g[,1:min(100, ncol(g))])

# Convert to long format for ggplot
cor_melt <- melt(cor_matrix, varnames = c("Marker1", "Marker2"), value.name = "Correlation")

p2 <- ggplot(cor_melt, aes(x = Marker1, y = Marker2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Corr.") +
  theme_minimal() +
  labs(title = "LD Map - First 100 Markers", x = "Marker", y = "Marker") 

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)


#--> blocks defined
# b=nblock: highest marker density, b=1: lowest density or pedigree
# block definitions, inds in each block
# WARNING: ids must be in same order in each G
b_ids = list(1:nind, 
             (nFounder*nOff+nFounder+1):nind, # array genotyping starts at third generation inds
             (nind-nFounder*nOff+1):nind) # last generation is sequenced
nblock = length(b_ids)

#--> chip markers: contains marker ids included in each block array, 0 for pedigree
# - first block has no markers, only pedigree
# - second block ids genotyped with random 100 markers
# - third block genotyped with all markers (sequence)
mk_ids = list(0, sample(1:nmkr, 1000), 1:nmkr)

#--> obtain relationship matrices
Glist = list()
for (b in 1:nblock) {
  # numerator relationship matrix
  if (length(mk_ids[[b]])==1) {
    Glist[[b]] = doA(ped) # marker GRM for given inds and mkrs
  } else {
    Glist[[b]] = GRM(g[b_ids[[b]], mk_ids[[b]]])
    diag(Glist[[b]]) = diag(Glist[[b]]) + 0.5 # to stabilize
  }
}

H    <- doHgss(Glist)
Hinv <- doHinv(Glist)

# Proof that Hinv is the inverse of H
I = H %*% Hinv
print(diag(I))
print(round(I[1:10,1:10],6))

tst <- b_ids[[nblock]] # ids to be predicted those with highest density

# Predictive ability in individuals genotyped with highest density markers
Hgss   <- doHgss(Glist)
blup_gss <- BLUP(Hgss, y, tst, h2)
yHat <- blup_gss$uhat[tst]
print(cor(y[tst], yHat))

# Predictive ability if same individuals were genotyped with SNP array
H = doHgss(Glist[-nblock]) # computes H without sequence (standard SS)
blup_h = BLUP(H, y, tst, h2)
yHat <- blup_h$uhat[tst]

# only pedigree
A = doA(ped)
blup_a = BLUP(A, y, tst, h2)
yHat = blup_a$uhat[tst]

## all sequenced
G=GRM(g)
blup_s = BLUP(G, y, tst, h2)
yHat <- blup_s$uhat[tst]


## plot A elements
ids = sample(1:(nind*(nind-1)/2), 5000)

library(ggplot2)
library(cowplot)
library(gridExtra)

# Relationship matrices
# Calculate correlations
c1 <- cor(G[lower.tri(G)][ids], A[lower.tri(H)][ids])
c2 <- cor(G[lower.tri(G)][ids], Hgss[lower.tri(H)][ids])
c3 <- cor(G[lower.tri(G)][ids], H[lower.tri(H)][ids])
c4 <- cor(Hgss[lower.tri(G)][ids], H[lower.tri(H)][ids])

# Create individual plots with labels in the title
p1 <- ggplot(data.frame(x = G[lower.tri(G)][ids], y = A[lower.tri(H)][ids]), 
             aes(x = x, y = y)) +
  geom_point(alpha = 0.5, color = "firebrick", size = 1.5) +
  labs(x = "S elements", y = "A elements", 
       title = paste0("A) Seq vs A, cor = ", round(c1, 3))) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))

p2 <- ggplot(data.frame(x = G[lower.tri(G)][ids], y = H[lower.tri(H)][ids]), 
             aes(x = x, y = y)) +
  geom_point(alpha = 0.5, color = "firebrick", size = 1.5) +
  labs(x = "S elements", y = "Hss elements", 
       title = paste0("B) Seq vs Hss, cor = ", round(c3, 3))) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))

p3 <- ggplot(data.frame(x = G[lower.tri(G)][ids], y = Hgss[lower.tri(H)][ids]), 
             aes(x = x, y = y)) +
  geom_point(alpha = 0.5, color = "firebrick", size = 1.5) +
  labs(x = "S elements", y = "Hgss elements", 
       title = paste0("C) Seq vs Hgss, cor = ", round(c2, 3))) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))

p4 <- ggplot(data.frame(x = Hgss[lower.tri(G)][ids], y = H[lower.tri(H)][ids]), 
             aes(x = x, y = y)) +
  geom_point(alpha = 0.5, color = "firebrick", size = 1.5) +
  labs(x = "Hgss elements", y = "Hss elements", 
       title = paste0("D) Hgss vs Hss, cor = ", round(c4, 3))) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))

# Arrange plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Corrs between true and predicted breeding values
c1 <- cor(tbv[tst], blup_a$uhat[tst])
c2 <- cor(tbv[tst], blup_s$uhat[tst])
c3 <- cor(tbv[tst], blup_gss$uhat[tst])
c4 <- cor(tbv[tst], blup_h$uhat[tst])

p1 <- ggplot(data.frame(x = tbv[tst], y = blup_a$uhat[tst]), 
             aes(x = x, y = y)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
  labs(x = "TBV", y = "Pedigree", 
       title = paste0("Pedigree: cor = ", round(c1, 3))) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"))

p2 <- ggplot(data.frame(x = tbv[tst], y = blup_s$uhat[tst]), 
             aes(x = x, y = y)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
  labs(x = "TBV", y = "Sequence", 
       title = paste0("Sequence: cor = ", round(c2, 3))) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"))

p3 <- ggplot(data.frame(x = tbv[tst], y = blup_gss$uhat[tst]), 
             aes(x = x, y = y)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
  labs(x = "TBV", y = "Generalized Single Step", 
       title = paste0("GSS: cor = ", round(c3, 3))) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"))


p4 <- ggplot(data.frame(x = tbv[tst], y = blup_h$uhat[tst]), 
             aes(x = x, y = y)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
  labs(x = "TBV", y = "Single Step", 
       title = paste0("SS: cor = ", round(c4, 3))) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"))
grid.arrange(p1, p2, p3, p4, ncol = 2)



