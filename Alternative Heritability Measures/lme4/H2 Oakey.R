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
C22.g <- C.inv[levels(dat$gen), levels(dat$gen)] # subset of C.inv that refers to genotypic BLUPs


ED    <- diag(n.g)-(solve(G.g)%*%C22.g)       # [see p. 813 bottom left in Oakey (2006)]
eM    <- eigen(ED)                            # obtain eigenvalues

############
# H2 Oakey #
############
# main method [see eq. (7) in Oakey (2006)]
H2Oakey <- sum(eM$values)/(n.g-1) 
H2Oakey # 0.8091336

# approximate method [see p. 813 top right in Oakey (2006)]
H2Oakey.approx <- 1 - psych::tr( as.matrix(solve(G.g)%*%C22.g / n.g ))
H2Oakey.approx # 0.7754197
