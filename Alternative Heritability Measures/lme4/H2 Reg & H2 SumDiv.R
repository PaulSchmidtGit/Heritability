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

dat$Mu <- 1 #Create pseudo intercept to obtain estimate for Mu

# random genotype effect
g.ran <- lmer(data    = dat,
              formula = yield ~ 0 + Mu + rep + (1|gen) + (1|rep:block))

# fixed genotype effect
g.fix <- lmer(data    = dat,
              formula = yield ~          rep +    gen  + (1|rep:block))

##########################
# Handle model estimates #
##########################

# Obtaining genotypic BLUPs
BLUPs <- as.data.table(ranef(g.ran)$gen, keep.rownames=T)
names(BLUPs) <- c("gen", "BLUP")

# Obtaining estimated marginal means based on genotypic BLUEs
BLUEs <- as.data.table(emmeans(g.fix, pairwise ~ gen)$emmeans)[,c("gen", "emmean")] # get estimated marginal means for all genotype pairs

# Overall mean in g.ran
Mu.ran <- as.data.table(emmeans(g.ran, "Mu"))[,c("emmean")] # get estimated marginal overall mean for g.ran model

# Combine BLUPs and emmeans, obtain scaled emmeans
Gpreds <- data.frame(gen    = BLUEs[,1], 
                     BLUP   = BLUPs[,2], 
                     emmean = BLUEs[,2], 
                     scaled.emmean = BLUEs[,2]-as.numeric(Mu.ran))

################
# H2 BLUP~BLUE #
################
H2reg <- lm(data    = Gpreds,
            formula = BLUP ~ 0 + scaled.emmean)$coefficients
H2reg #0.8178116

############# 
# H2 sumdiv #
#############
H2sumdiv <- sum(abs(Gpreds$BLUP)) / sum(abs(Gpreds$scaled.emmean))
H2sumdiv #0.8205183
