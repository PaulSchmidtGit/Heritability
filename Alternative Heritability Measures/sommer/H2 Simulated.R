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
vc.g    <- as.numeric(g.ran$var.comp[["gen"]])  
A.g     <- g.ran$ZETA[["gen"]]$K 
G.g     <- A.g*vc.g 

vc.rb   <- as.numeric(g.ran$var.comp[["rep:block"]])
A.rb    <- g.ran$ZETA[["rep:block"]]$K
G.rb    <- A.rb*vc.rb 

require(Matrix)
G <- bdiag(G.g,G.rb)
F <- G[diag(G)==vc.g,]
D <- G[diag(G)==vc.g, diag(G)==vc.g]; all(D == G.g)

C22.g   <- g.ran$PEV.u.hat[["gen"]][[1]]
C22.rb  <- g.ran$PEV.u.hat[["rep:block"]][[1]]
C22     <- bdiag(C22.g,C22.rb)

M       <- G - C22
inv.G   <- solve(G)
Q       <- F %*% inv.G %*% M %*% inv.G %*% t(F)
Omega   <- rbind(cbind(D,Q), cbind(Q,Q))
svdout  <- svd(Omega)
Gamma   <- svdout$u %*% diag(sqrt(svdout$d))

##############
# Simulation #
##############
n.g    <- nrow(Gamma)/2                       #;nrow(Gamma)/2==g.ran$dimos[["gen"]][2]
n.sim  <- 10000
h2     <- list()
R      <- list()

  for (i in 1:n.sim){
    z       <- rnorm(n=(2*n.g), mean=0, sd=1)
    w       <- Gamma %*% z
    g.hat   <- w[1      :   n.g ]
    g.true  <- w[(n.g+1):(2*n.g)]
    g.hat.s <- sort(g.hat, decreasing=T)
    #v       <- data.table(g.hat, g.true, ranks);  setorder(v, ranks)
    selmean <- list() 
      for (j in 1:n.g){
        selmean[[j]] <- mean(g.hat.s[1:j])
      }
    R[[i]]  <- unlist(selmean)
    h2[[i]] <- (t(g.hat) %*% g.true)**2 / (t(g.true) %*% g.true %*% t(g.hat) %*% g.hat)
  }

##########
# H2 Sim #
##########
H2Sim <- mean(unlist(h2))
H2Sim # 0.7722301


