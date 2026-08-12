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

##############
# Fit models #
##############
library(asreml)

dat$Mu <- 1 #Create pseudo intercept to obtain estimate for Mu 

# Genotype as random effect
g.ran <- asreml(fixed = yield ~  -1 + Mu + rep, 
                random=       ~ gen + rep:block, 
                data=dat,
                ran.order = "user") #Force "gen" as first random effect in asreml object - makes BLUP extraction easier
# Genotype as fixed effect
g.fix <- asreml(fixed = yield ~ -1 + Mu + gen + rep,
                random=       ~ rep:block, 
                data=dat)

##########################
# Handle model estimates #
##########################
# gentoypic BLUEs
BLUEs <- predict(g.fix, classify="gen")$pvals[,c('gen','predicted.value')]

# gentoypic BLUPS
BLUPs <- predict(g.ran, classify="gen", only="gen")$pvals[,c('gen','predicted.value')]

# Overall mean in g.ran
Mu.ran <- predict(g.ran, classify="Mu")$pvals$predicted.value
Mu.ran #4.479517

# Combine BLUPs and BLUEs, obtain scaled BLUEs
Gpreds <- data.frame(gen  = BLUEs[,1], 
                     BLUP = BLUPs[,2], 
                     BLUE = BLUEs[,2], 
                     scaled.BLUE = BLUEs[,2]-Mu.ran)

################
# H2 BLUP~BLUE #
################
H2reg <- lm(data    = Gpreds,
            formula = BLUP ~ 0 + scaled.BLUE)$coefficients
H2reg #0.8178116

#############
# H2 sumdiv #
#############
H2sumdiv <- sum(abs(Gpreds$BLUP))/sum(abs(Gpreds$scaled.BLUE))
H2sumdiv #0.8205183
