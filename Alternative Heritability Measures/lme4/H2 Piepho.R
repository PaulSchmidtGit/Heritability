rm(list = ls())
#######################
# Import example data #
#######################
library(agridat)
dat <- john.alpha

##############
# Fit models #
##############
library(lme4)
library(lmerTest)
library(emmeans)

# random genotype effect
g.ran <- lmer(data    = dat,
              formula = yield ~ rep + (1|gen) + (1|rep:block))

# fixed genotype effect
g.fix <- lmer(data    = dat,
              formula = yield ~ rep +    gen  + (1|rep:block))

##########################
# Handle model estimates #
##########################

# Genotypic variance component
vc.g <- as.data.table(VarCorr(g.ran))[grp=="gen", vcov]

# Obtaining adjusted means based on genotypic BLUEs
diffs.BLUE <- as.data.table(emmeans(g.fix, pairwise ~ gen)$contrasts) # get differences for all genotype pairs
vdBLUE.avg <- mean(diffs.BLUE$SE^2) # mean variance of a difference = mean squared standard error of a difference

#############
# H2 Piepho #
#############
H2.p <- vc.g/(vc.g + vdBLUE.avg/2)
H2.p #0.7966373