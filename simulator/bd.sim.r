library("TreeSim")

set.seed(123)

# generate birth-death trees
rept  <- 205
tmrca <- 1.0
birth <- 5.0
death <- 4.0
tree  <- sim.bd.age(tmrca, rept, birth, death, mrca=T)

# write trees with <50 internal nodes
for (i in 1:rept) {
	if (tree[[i]]$Nnode < 50) {
		write.tree(tree[[i]], file="bd.trees", append=T)
	}
}

# replace ending ':0;' with ';'
# sed -i 's/:0;/;/g' bd.trees      
