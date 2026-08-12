# =============================================================================
# ASReml-R 4.2.x version - COMMUNITY-CONTRIBUTED, NOT authored or tested by us
# =============================================================================
# Updated from the original ASReml-R 3.0 syntax to ASReml-R 4.2.x by GitHub
# user @mpeddu, submitted as pull request #5:
#   https://github.com/PaulSchmidtGit/Heritability/pull/5
#
# The original ASReml-R 3.0 version of this script is in the parent folder.
#
# NOTE: This code has NOT been run or verified by the repository maintainer,
# who does not have an ASReml licence. It is provided as-is; please review and
# test it yourself before relying on the results.
# =============================================================================

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
vc.g <- summary(g.ran)$varcomp['gen','component']
vc.g #0.142902

# Mean variance of a difference of two genotypic BLUPs
vdBLUP.mat <- predict(g.ran, classify="gen", only="gen", sed=TRUE)$sed^2 # obtain squared s.e.d. matrix 
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
vdBLUP.avg #0.05455038

#############
# H2 Cullis #
#############
H2Cullis <- 1 - (vdBLUP.avg / 2 / vc.g)
H2Cullis #0.8091336
