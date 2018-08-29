rm(list = ls())
#######################
# Import example data #
#######################
library(agridat)
dat <- john.alpha

##############
# Fit models #
##############
require(asreml)
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
gen.var <- summary(g.ran)$varcomp['gen!gen.var','component']
gen.var #0.142902

# Mean variance of a difference of two GBLUEs
vdBLUE.mat <- predict(g.fix, classify="gen", sed=TRUE)$pred$sed^2 # obtain squared s.e.d. matrix 
vdBLUE.avg <- mean(vdBLUE.mat[upper.tri(vdBLUE.mat, diag=FALSE)]) # take mean of upper triangle
vdBLUE.avg #0.07010875

#############
# H2 Piepho #
#############
H2.p <- gen.var/(gen.var + vdBLUE.avg/2)
H2.p #0.803017



