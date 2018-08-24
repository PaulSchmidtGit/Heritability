rm(list = ls())
#######################
# Import example data #
#######################
library(agridat)
dat <- john.alpha

#############
# Fit model #
#############
library(sommer)
# Genotype as random effect
g.ran <- mmer2(fixed  = yield ~ rep,
               random =       ~ gen + rep:block,
               data   = dat)

##########################
# Handle model estimates #
##########################
n.g   <- g.ran$dimos[["gen"]][2]              # number of genotpyes
vc.g  <- as.numeric(g.ran$var.comp[["gen"]])  # genetic variance component
A     <- g.ran$ZETA[["gen"]]$K                # 
G.g   <- A*vc.g                               # 
C22.g <- g.ran$PEV.u.hat[["gen"]][[1]]        #
M     <- diag(n.g)-(solve(G.g)%*%C22.g)       # [see p. 813 bottom left in Oakey (2006)]
eM    <- eigen(M)                             # obtain eigenvalues

############
# H2 Oakey #
############
# main method [see eq. (7) in Oakey (2006)]
H2Oakey <- sum(eM$values)/(n.g-1) 
H2Oakey # 0.8091336

library(psych) # to compute trace of a matrix
# approximate method [see p. 813 top right in Oakey (2006)]
H2Oakey.approx <- 1 - tr( as.matrix(solve(G.g)%*%C22.g / n.g) )
H2Oakey.approx # 0.7754197
