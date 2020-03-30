library(agridat)
library(dplyr)
library(sommer) 
library(purrr)

### get example data
dat <- john.alpha

### fit model
# random genotype effect
g.ran <- mmer(fixed  = yield ~ rep,
              random =       ~ gen + rep:block,
              data   = dat)

### handle model estimates
vc_g     <- g.ran %>% pluck("sigma") %>% pluck("gen") %>% as.numeric                # genetic variance component
n_g      <- g.ran %>% pluck("U")     %>% pluck("gen") %>% pluck("yield") %>% length # number of genotypes
C22_g    <- g.ran %>% pluck("PevU")  %>% pluck("gen") %>% pluck("yield")            # Prediction error variance matrix for genotypic BLUPs
trC22_g  <- psych::tr(as.matrix(C22_g))                                             # trace
vdBLUP_g <- 2/n_g * (trC22_g - (sum(C22_g)-trC22_g) / (n_g-1))                      # Mean variance of a difference between genotypic BLUPs

### H2 Cullis
H2Cullis <- 1-(vdBLUP_g / 2 / vc_g)
H2Cullis #0.8091336