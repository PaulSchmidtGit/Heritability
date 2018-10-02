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
# to obtain var-cov-matrix for BLUPs (of gen effect).
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
F <- G[diag(G)==vc.g,]
D <- G[diag(G)==vc.g, diag(G)==vc.g]; all(D == G.g)

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
PEV   <- C.inv[-c(1:3), -c(1:3)]                 # subset of C.inv that refers to all BLUPS

# Gamma
M       <- G - PEV
inv.G   <- solve(G)
Q       <- F %*% inv.G %*% M %*% inv.G %*% t(F)
Omega   <- rbind(cbind(D,Q), cbind(Q,Q))
svdout  <- svd(Omega)
Gamma   <- svdout$u %*% diag(sqrt(svdout$d))

##############
# Simulation #
##############
n.sim  <- 10000
h2     <- list()
R      <- list()

for (i in 1:n.sim){
  z       <- rnorm(n=(2*n.g), mean=0, sd=1)
  w       <- Gamma %*% z
  g.hat   <- w[1      :   n.g ]
  g.true  <- w[(n.g+1):(2*n.g)]
  g.hat.s <- sort(g.hat, decreasing=T)
  #v       <- data.table(g.hat, g.true, ranks);  setorder(v, ranks)
  selmean <- list() 
  for (j in 1:n.g){
    selmean[[j]] <- mean(g.hat.s[1:j])
  }
  R[[i]]  <- unlist(selmean)
  h2[[i]] <- (t(g.hat) %*% g.true)**2 / (t(g.true) %*% g.true %*% t(g.hat) %*% g.hat)
}

##########
# H2 Sim #
##########
H2Sim <- mean(unlist(h2))
H2Sim # 0.7719582