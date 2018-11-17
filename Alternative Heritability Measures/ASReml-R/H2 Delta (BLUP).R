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

# BLUPs for genotype main effect
g.pred  <- predict(g.ran, classify="gen", only="gen", sed=T, vcov=T)$pred
BLUPs.g <- data.table(g.pred$pvals[,c(1,2)]); names(BLUPs.g) <- c("gen","BLUP")

##########################
# Handle model estimates #
##########################
library(data.table)
# Genotype information
list.g <- levels(dat$gen) # list of genotype names
n.g    <- length(list.g)  # number of genotypes

# G.g (i.e. estimated G matrix of genotype main effect)
vc.g <- summary(g.ran)$varcomp['gen!gen.var','component'] # VC genotype main effect
G.g.wide <- diag(1, n.g) * vc.g; dimnames(G.g.wide) <- list(list.g, list.g) # G.g matrix
G.g.long <- data.table(reshape::melt(G.g.wide)); names(G.g.long) <- c("gen1","gen2","sigma") # G.g matrix in long format

# Variance of a difference between genotypic BLUPs (based on C22.g/PEV matrix)
vd.g.wide <- g.pred$sed^2; dimnames(vd.g.wide) <- list(list.g, list.g) # C22.g matrix
vd.g.long <- data.table(reshape::melt(vd.g.wide)); names(vd.g.long) <- c("gen1","gen2","vd") # C22.g matrix in long format

# merge BLUPs, G.g and C22.g information into "H2D.blup" table
g.var <- G.g.long[gen1==gen2, .(gen1, sigma)]; names(g.var) <- c("gen","var") # variances G.g
g.cov <- G.g.long[gen1!=gen2]                ; names(g.cov) <- c("gen1","gen2","cov") # covariances of G.g
H2D.blup   <- merge(vd.g.long, g.cov, all=T)
for (i in 1:2){
  temp <- merge(BLUPs.g, g.var, by="gen"); names(temp) <- c(paste0(names(temp),i)) # merge BLUPs and variances for each genotye
  H2D.blup  <- merge(H2D.blup, temp, by=paste0("gen",i)) # merge this for both gen1 and gen2, respectively, to result table
}

# formatting
setcolorder(H2D.blup, c("gen1","BLUP1","gen2","BLUP2","var1","var2","cov","vd"))
H2D.blup <- H2D.blup[order(gen1, gen2)]
H2D.blup[, i.is.j     := gen1==gen2]                          # i=j is not a pair
H2D.blup[, i.larger.j := as.numeric(gen1) > as.numeric(gen2)] # i>j is a duplicate pair  

### Compute H2 Delta based on BLUPs
# H2 Delta ij
H2D.blup[i.is.j==FALSE, Numerator   := var1 + var2 - 2*cov - vd]
H2D.blup[i.is.j==FALSE, Denominator := var1 + var2 - 2*cov     ]
H2D.blup[i.is.j==FALSE,                     H2D.ij := Numerator / Denominator]
# H2 Delta i.
H2D.blup[i.is.j==FALSE,                     H2D.i  := mean(H2D.ij), by="gen1"]
# H2 Delta ..
H2D.blup[i.is.j==FALSE & i.larger.j==FALSE, H2D    := mean(H2D.ij)]

#######################
### H2 Delta (BLUP) ###
#######################

# H2 Delta .. (overall H2)
H2D.. <- as.numeric(unique(na.omit(H2D.blup[,.(H2D)])))
H2D..
# H2 Delta i. (mean H2 per genotype)
H2Di. <- unique(na.omit(H2D.blup[,.(gen1, H2D.i)]))
H2Di.
# H2 Delta ij (H2 per genotype pair)
H2Dij <- unique(na.omit(H2D.blup[i.larger.j==FALSE,.(gen1, gen2, H2D.ij)]))
mergeit <- H2D.blup[,c("gen1","gen2","var1","var2","cov","vd")]
H2Dij <- merge(H2Dij, mergeit, by=c("gen1", "gen2"), all = FALSE)
H2Dij