library(agridat)
library(dplyr)
library(emmeans)
library(lme4) 
library(lmerTest) 
library(purrr)

### get example data
dat <- john.alpha

### fit models
# random genotype effect
g.ran <- lmer(data    = dat,
              formula = yield ~ rep + (1|gen) + (1|rep:block))

# fixed genotype effect
g.fix <- lmer(data    = dat,
              formula = yield ~ rep +    gen  + (1|rep:block))

### handle model estimates
# genotypic variance component
vc.g <- g.ran %>% 
  VarCorr %>% 
  as_tibble %>% 
  filter(grp=="gen") %>% 
  pull(vcov) # 0.1429021

# mean variance of a difference between genotypes
vdBLUE.avg <- g.fix %>% 
  emmeans(pairwise ~ gen) %>% 
  pluck("contrasts") %>% 
  as_tibble %>% 
  mutate(Var=SE^2) %>% 
  pull(Var) %>% 
  mean # vdBLUE.avg

### H2 Piepho
H2.p <- vc.g/(vc.g + vdBLUE.avg/2)
H2.p # 0.7966375