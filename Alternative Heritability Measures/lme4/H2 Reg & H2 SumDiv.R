library(agridat)
library(emmeans)
library(lme4) 
library(tidyverse)


# get example data --------------------------------------------------------
dat <- agridat::john.alpha %>% 
  mutate(Mu = 1) # Create dummy column for pseudo intercept in order to obtain estimate for Mu


# fit model ---------------------------------------------------------------
# random genotype effect
g_ran <- lmer(data    = dat,
              formula = yield ~ 0 + Mu + rep + (1|gen) + (1|rep:block)) # Default intercept set to 0, pseudo intercept Mu instead

# fixed genotype effect
g_fix <- lmer(data    = dat,
              formula = yield ~          rep +    gen  + (1|rep:block))


# handle model estimates --------------------------------------------------
# genotypic BLUPs
g_BLUPs <- g_ran %>% 
  ranef %>% as_tibble %>% 
  rename(BLUP=condval) %>%
  filter(grpvar=="gen") %>%
  mutate(gen = grp %>% as.character %>% as.factor) %>% 
  select(gen, BLUP)

# estimated marginal means (a.k.a. adjusted means) based on genotypic BLUEs
g_EMMs <- g_fix %>% 
  emmeans("gen") %>% 
  as_tibble %>% 
  select(gen, emmean)

# Overall mean in g_ran
Mu_ran <- g_ran %>% emmeans("Mu") %>% as_tibble %>% pull(emmean)

# Combine BLUPs and emmeans, compute scaled emmeans
Gpreds <- left_join(g_EMMs, g_BLUPs, by="gen") %>% 
  mutate(scaled_emmean = emmean - Mu_ran)


# H2 BLUP~BLUE ------------------------------------------------------------
H2reg <- lm(data    = Gpreds,
            formula = BLUP ~ 0 + scaled_emmean) %>% pluck("coefficients")
H2reg #0.8178116


# H2 sumdiv ---------------------------------------------------------------
H2sumdiv <- sum(abs(Gpreds$BLUP)) / sum(abs(Gpreds$scaled_emmean))
H2sumdiv #0.8205183