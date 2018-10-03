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
vc.g     <- as.numeric(g.ran$var.comp[["gen"]]) # genetic variance component
n.g      <- g.ran$dimos[["gen"]][2]             # number of genotypes
C22.g    <- g.ran$PEV.u.hat[["gen"]][[1]]       #
trC22.g  <- psych::tr(as.matrix(C22.g))         # trace
vdBLUP.g <- 2/n.g*(trC22.g-(sum(C22.g)-trC22.g)/(n.g-1)) # Mean variance of a difference of two genotypic BLUPs

#############
# H2 Cullis #
#############
H2Cullis <- 1-(vdBLUP.g / 2 / vc.g)
H2Cullis #0.8091336