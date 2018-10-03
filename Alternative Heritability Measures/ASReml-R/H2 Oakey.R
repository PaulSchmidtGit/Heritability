rm(list = ls())
#######################
# Import example data #
#######################
library(agridat)
dat <- john.alpha

#############
# Fit model #
#############
library(asreml)
# Genotype as random effect
g.ran <- asreml(fixed = yield ~ rep, 
                random=       ~ gen + rep:block, 
                data=dat)

##########################
# Handle model estimates #
##########################
n.g   <- as.numeric(g.ran$noeff["gen"]) # number of genotypes
vc.g  <- summary(g.ran)$varcomp['gen!gen.var','component'] # genetic variance component
G.g   <- diag(1,n.g)*vc.g               # note that this is manually created for simple diag structre
C22.g <- predict(g.ran, classify="gen", only="gen", vcov=TRUE)$pred$vcov
M     <- diag(n.g)-(solve(G.g)%*%C22.g) # [see p. 813 bottom left in Oakey (2006)]
eM    <- eigen(M)                       # obtain eigenvalues

############
# H2 Oakey #
############
# main method [see eq. (7) in Oakey (2006)]
H2Oakey <- sum(eM$values)/(n.g-1) 
H2Oakey # 0.8091336

library(psych) # to compute trace of a matrix
# approximate method [see p. 813 top right in Oakey (2006)]
H2Oakey.approx <- 1 - tr( solve(G.g)%*%C22.g / n.g ) 
H2Oakey.approx # 0.7754197

