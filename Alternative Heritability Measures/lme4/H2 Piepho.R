library(agridat)
library(emmeans)
library(lme4) 
library(tidyverse)

# get example data --------------------------------------------------------
dat <- agridat::john.alpha


# fit model ---------------------------------------------------------------
# random genotype effect
g_ran <- lmer(data    = dat,
              formula = yield ~ rep + (1|gen) + (1|rep:block))

# fixed genotype effect
g_fix <- lmer(data    = dat,
              formula = yield ~ rep +    gen  + (1|rep:block))


# handle model estimates --------------------------------------------------
# genotypic variance component
vc.g <- g_ran %>% 
  VarCorr %>% 
  as_tibble %>% 
  filter(grp=="gen") %>% 
  pull(vcov) # 0.1429021

# mean variance of a difference between genotypes
vdBLUE.avg <- g_fix %>% 
  emmeans(pairwise ~ gen) %>% 
  pluck("contrasts") %>% 
  as_tibble %>% 
  mutate(Var=SE^2) %>% 
  pull(Var) %>% 
  mean # 0.07295899


# H2 Piepho ---------------------------------------------------------------
H2.p <- vc.g/(vc.g + vdBLUE.avg/2)
H2.p # 0.7966375