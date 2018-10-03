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
g.ran <- asreml(fixed = yield ~       rep, 
                random=       ~ gen + rep:block, 
                data=dat)

##########################
# Handle model estimates #
##########################
# Genetic variance component
vc.g <- summary(g.ran)$varcomp['gen!gen.var','component']
vc.g #0.142902

# Mean variance of a difference of two genotypic BLUPs
vdBLUP.mat <- predict(g.ran, classify="gen", only="gen", sed=TRUE)$pred$sed^2 # obtain squared s.e.d. matrix 
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
vdBLUP.avg #0.05455038

#############
# H2 Cullis #
#############
H2Cullis <- 1 - (vdBLUP.avg / 2 / vc.g)
H2Cullis #0.8091336
