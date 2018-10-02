rm(list = ls())
library(data.table)
#######################
# Import example data #
#######################
library(agridat)
dat <- john.alpha

library(lme4)
library(lmerTest)

g.ran <- lmer(data    = dat,
            formula = yield ~ rep + (1|gen) + (1|rep:block))

# manually reconstruct mixed model equation for this specific example
# to obtain var-cov-matrix for BLUPs of gen effect.
vc   <- as.data.table(VarCorr(g.ran)) # extract estimated variance components (vc)

# R = varcov-matrix for error term
n    <- length(summary(g.ran)$residuals) # numer of observations
vc.e <- vc[grp=="Residual", vcov]      # error vc
R    <- diag(n)*vc.e                   # R matrix = I * vc.e

# G = varcov-matrx for all random effects
  # varcov-matrix for genotypic effect
  n.g  <- summary(g.ran)$ngrps["gen"]    # number of genotypes
  vc.g <- vc[grp=="gen", vcov]         # genotypic vc
  G.g  <- diag(n.g)*vc.g               # gen part of G matrix = I * vc.g
  
  # varcov-matrix for incomplete block effect
  n.b  <- summary(g.ran)$ngrps["rep:block"] # number of incomplete blocks
  vc.b <- vc[grp=="rep:block", vcov]      # incomplete block vc
  G.b  <- diag(n.b)*vc.b                  # incomplete block part of G matrix = I * vc.b
  
G <- bdiag(G.g, G.b) # G is blockdiagonal with G.g and G.b

# Design Matrices
X <- as.matrix(getME(g.ran, "X")) # Design matrix fixed effects
Z <- as.matrix(getME(g.ran, "Z")) # Design matrix random effects

# Mixed Model Equation (HENDERSON 1986; SEARLE et al. 2006)
C11 <- t(X) %*% solve(R) %*% X
C12 <- t(X) %*% solve(R) %*% Z
C21 <- t(Z) %*% solve(R) %*% X
C22 <- t(Z) %*% solve(R) %*% Z + solve(G) 

C <- as.matrix(rbind(cbind(C11, C12),  # Combine components into one matrix C
                     cbind(C21, C22)))

# Mixed Model Equation Solutions 
C.inv <- solve(C)                                # Inverse of C
PEV.g <- C.inv[levels(dat$gen), levels(dat$gen)] # subset of C.inv that refers to genotypic BLUPs

# Mean variance of BLUP-difference from PEV matrix of genotypic BLUPs
one        <- t(t(rep(1, n.g)))                # vector of 1s
P.mu       <- diag(n.g, n.g) - one %*% t(one)  # P.mu = matrix that centers for overall-mean
vdBLUP.sum <- psych::tr(P.mu %*% PEV.g)        # sum of all variance of differences = trace of P.mu*PEV.g
vdBLUP.avg <- vdBLUP.sum * (2/(n.g*(n.g-1)))   # mean variance of BLUP-difference = divide sum by number of genotype pairs

#############
# H2 Cullis #
#############
H2Cullis <- 1 - (vdBLUP.avg / 2 / vc.g)
H2Cullis #0.8091336
