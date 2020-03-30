library(agridat)
library(dplyr)
library(lme4) 
library(lmerTest) 
library(psych)
library(purrr)

### get example data
dat <- john.alpha

### fit model
# random genotype effect
g.ran <- lmer(data    = dat,
              formula = yield ~ rep + (1|gen) + (1|rep:block))

### handle model estimates
# to my knowledge, lme4 does not offer a function to
# extract variance-covariance-matrices for BLUPs (a.k.a. prediction error variance [PEV] matrix).
# therefore, I here manually reconstruct mixed model equation for this specific example.
# notice that this solution therefore only works for this specific model!

vc <- g.ran %>% VarCorr %>% as_tibble # extract estimated variance components (vc)

# R = varcov-matrix for error term
n <- g.ran %>% summary %>% pluck(residuals) %>% length # numer of observations
vc_e <- vc %>% filter(grp=="Residual") %>% pull(vcov)  # error vc
R    <- diag(n)*vc_e                                   # R matrix = I_n * vc_e

# G = varcov-matrx for all random effects
# subset of G regarding genotypic effects
n_g  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("gen") # number of genotypes
vc_g <- vc %>% filter(grp=="gen") %>% pull(vcov)              # genotypic vc
G_g  <- diag(n_g)*vc_g                                        # gen part of G matrix = I * vc.g

# subset of G regarding incomplete block effects
n_b  <- g.ran %>% summary %>% pluck("ngrps") %>% pluck("rep:block") # number of incomplete blocks
vc_b <- vc %>% filter(grp=="rep:block") %>% pull(vcov)              # incomplete block vc
G_b  <- diag(n_b)*vc_b                                              # incomplete block part of G matrix = I * vc.b

G <- bdiag(G_g, G_b) # G is blockdiagonal with G.g and G.b in this example

# Design Matrices
X <- g.ran %>% getME("X") %>% as.matrix # Design matrix fixed effects
Z <- g.ran %>% getME("Z") %>% as.matrix # Design matrix random effects

# Mixed Model Equation (HENDERSON 1986; SEARLE et al. 2006)
C11 <- t(X) %*% solve(R) %*% X
C12 <- t(X) %*% solve(R) %*% Z
C21 <- t(Z) %*% solve(R) %*% X
C22 <- t(Z) %*% solve(R) %*% Z + solve(G) 

C <- rbind(cbind(C11, C12),  
           cbind(C21, C22)) %>% as.matrix # Combine components into one matrix C

# Mixed Model Equation Solutions 
C_inv <- C %>% solve                             # Inverse of C
C22_g <- C_inv[levels(dat$gen), levels(dat$gen)] # subset of C.inv that refers to genotypic BLUPs


ED    <- diag(n_g) - (solve(G_g) %*% C22_g)       # [see p. 813 bottom left in Oakey (2006)]
eM    <- ED %>% eigen                             # obtain eigenvalues


### H2 Oakey
# main method [see eq. (7) in Oakey (2006)]
H2Oakey <- sum(eM$values)/(n_g-1) 
H2Oakey # 0.8091336

# approximate method [see p. 813 top right in Oakey (2006)]
H2Oakey.approx <- 1 - tr( as.matrix(solve(G_g) %*% C22_g / n_g ) )
H2Oakey.approx # 0.7754197
