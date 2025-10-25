library(ape)
library(phytools)

setwd("~/Documents/Research/2025-ContChars/homo-longi")

# dated tree from Ni et al 2021
tree_ni21 <- read.nexus("ni2021.con.tre")
# plot(tree_ni21, main = "Ni et al. 2021")

# continuous and discrete, unlinked WN clock
tree_c_wn2 <- read.nexus("tipdating_wn_2r/data_both.con.tre")

# discretized and discrete, unlinked WN clock
tree_d_wn2 <- read.nexus("tipdating_wn_2r/data_disc.con.tre")

# plot tanglegrams
plot(cophylo(tree_c_wn2, tree_ni21), fsize = 0.6)
plot(cophylo(tree_ni21, tree_d_wn2), fsize = 0.6)  # not included

# compute the Robinson-Foulds distance for reference
library(phangorn)
RF.dist(tree_ni21, tree_c_wn2)
RF.dist(tree_ni21, tree_d_wn2)
