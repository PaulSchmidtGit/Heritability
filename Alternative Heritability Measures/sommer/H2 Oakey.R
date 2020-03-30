library(agridat)
library(dplyr)
library(sommer)
library(psych)
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
G_g      <- diag(n_g)*vc_g                                                          # subset of G regarding genotypic effects = I * vc.g
C22_g    <- g.ran %>% pluck("PevU")  %>% pluck("gen") %>% pluck("yield")            # Prediction error variance matrix for genotypic BLUPs

ED       <- diag(n_g) - (solve(G_g) %*% C22_g)       # [see p. 813 bottom left in Oakey (2006)]
eM       <- ED %>% eigen                              # obtain eigenvalues

############
# H2 Oakey #
############
# main method [see eq. (7) in Oakey (2006)]
H2Oakey <- sum(eM$values)/(n_g-1) 
H2Oakey # 0.8091336

# approximate method [see p. 813 top right in Oakey (2006)]
H2Oakey_approx <- 1 - tr( as.matrix(solve(G_g) %*% C22_g / n_g ) )
H2Oakey_approx # 0.7754197
