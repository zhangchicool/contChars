library(ape)
library(Quartet)
library(TreeDist)

setwd(system("pwd", intern=T))

# testing purpose only
# true_t <- read.tree(text="(A, B, (C, (D, E)));")
# estm_1 <- read.tree(text="(A, B, C, D, E);")
# estm_2 <- read.tree(text="(A, E, (D, (C, B)));")
# estm_3 <- read.tree(text="(A, B, C, (D, E));")

args = commandArgs(trailingOnly=TRUE)
# read in the two trees for comparison
true_t <- read.tree(file="bd.tre")
estm_t <- read.nexus(file=args[1])

# resolution
if (is.rooted(estm_t)) {
  reso <- (estm_t$Nnode -1) / (Ntip(estm_t) -2)
} else {
  reso <- (estm_t$Nnode -1) / (Ntip(estm_t) -3)
}

# quartet similarity metrics
qstatus <- QuartetStatus(true_t, estm_t)
# quartet divergence: (2d + r1 + r2) / 2Q
d_qt <- QuartetDivergence(qstatus, similarity=F)
# additionally, strict joint assertions (SJA): s / (s + d)
sja <- SimilarityMetrics(qstatus)[,"StrictJointAssertions"]
# SimilarityMetrics(qstatus)
# VisualizeQuartets(true_t, estm_t)

# mutual clustering information distance
d_ci <- TreeDistance(true_t, estm_t)
# ClusteringInfoDistance(true_t, estm_t, reportMatching=T)
# VisualizeMatching(ClusteringInfoDistance, true_t, estm_t)

# Robinson-Foulds distance (using phylo information content)
d_rf <- InfoRobinsonFoulds(true_t, estm_t, normalize=T)
# VisualizeMatching(InfoRobinsonFoulds, true_t, estm_t)
# RobinsonFoulds(true_t, estm_t, reportMatching=T)

# write the metrics to file
write(c(reso, sja, d_qt, d_ci, d_rf),
      file="dist_t.txt", append=T)

# reference
# https://ms609.github.io/Quartet/articles/Using-Quartet.html
# https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html
# https://ms609.github.io/TreeDist/articles/Generalized-RF.html
