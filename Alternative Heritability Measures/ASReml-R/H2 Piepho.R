rm(list = ls())
#######################
# Import example data #
#######################
library(agridat)
dat <- john.alpha

##############
# Fit models #
##############
library(asreml)
# Genotype as random effect
g.ran <- asreml(fixed = yield ~       rep, 
                random=       ~ gen + rep:block, 
                data=dat)
# Genotype as fixed effect
g.fix <- asreml(fixed = yield ~ gen + rep,
                random=       ~       rep:block, 
                data=dat)

##########################
# Handle model estimates #
##########################
# Genetic variance component
vc.g <- summary(g.ran)$varcomp['gen!gen.var','component']
vc.g #0.142902

# Mean variance of a difference of two genotypic BLUEs
vdBLUE.mat <- predict(g.fix, classify="gen", sed=TRUE)$pred$sed^2 # obtain squared s.e.d. matrix 
vdBLUE.avg <- mean(vdBLUE.mat[upper.tri(vdBLUE.mat, diag=FALSE)]) # take mean of upper triangle
vdBLUE.avg #0.07010875

#############
# H2 Piepho #
#############
H2.p <- vc.g/(vc.g + vdBLUE.avg/2)
H2.p #0.803017



