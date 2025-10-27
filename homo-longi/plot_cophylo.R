library(ape)
library(phytools)

setwd("~/Documents/Research/2025-ContChars/homo-longi")

# dated tree from Ni et al 2021
tree_ni21 <- read.nexus("ni2021.con.tre")

# continuous and discrete, linked WN clock
tree_c_wn <- read.nexus("tipdating_wn/data_both.con.tre")

# discretized and discrete, linked WN clock
tree_d_wn <- read.nexus("tipdating_wn/data_disc.con.tre")

# continuous and discrete, unlinked WN clock
tree_c_wn2 <- read.nexus("tipdating_wn_2r/data_both.con.tre")

# discretized and discrete, unlinked WN clock
tree_d_wn2 <- read.nexus("tipdating_wn_2r/data_disc.con.tre")

# plot tanglegrams
plot(cophylo(tree_c_wn2, tree_ni21), fsize = 0.6)


# plot treespaces
library(rwty)

trees_c1 <- load.multi("tipdating_wn/cont", trim = 100)
trees_d1 <- load.multi("tipdating_wn/disc", trim = 100)
trees_c2 <- load.multi("tipdating_wn_2r/cont", trim = 100)
trees_d2 <- load.multi("tipdating_wn_2r/disc", trim = 100)

p <- makeplot.treespace(trees, burnin = 80, fill.color = "lnLike")
p$treespace.points.plot
# p$treespace.heatmap
