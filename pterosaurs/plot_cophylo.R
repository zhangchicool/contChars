library(ape)
library(phytools)

setwd("~/Documents/Research/2025-ContChars/pterosaurs")

# most parsimonious tree from Andres 2021 (TODO)
tree_mp <- read.nexus("andres21.mp.tre")

# continuous and discrete, linked ILN clock
tree_c_iln <- read.nexus("tipdating_iln/data_norm.con.tre")

# discretized and discrete, linked ILN clock
tree_d_iln <- read.nexus("tipdating_iln/data_disc.con.tre")

# continuous and discrete, unlinked ILN clock
tree_c_il2 <- read.nexus("tipdating_iln_2r/data_norm.con.tre")

# discretized and discrete, unlinked ILN clock
tree_d_il2 <- read.nexus("tipdating_iln_2r/data_disc.con.tre")

# plot tanglegrams
plot(cophylo(tree_c_iln, tree_mp), fsize = 0.5)


# plot treespaces
library(rwty)

trees_c1 <- load.multi("tipdating_iln/cont", trim = 100)
trees_d1 <- load.multi("tipdating_iln/disc", trim = 100)
trees_c2 <- load.multi("tipdating_iln_2r/cont", trim = 100)
trees_d2 <- load.multi("tipdating_iln_2r/disc", trim = 100)

p <- makeplot.treespace(trees, burnin = 100, fill.color = "lnLike")
p$treespace.points.plot
# p$treespace.heatmap
