# setwd("~/Downloads")

cont <- read.table(file="continuous.txt")

nstate <- 6
nrow   <- 177
ncol   <- 52

disc <- matrix(nrow=nrow, ncol=ncol, byrow=T)

for (j in 1:ncol) {
  c <- as.numeric(cont[,j])
  min <- summary(c)[1] - 0.0001
  max <- summary(c)[6] + 0.0001
  
  for (i in 1:nrow) {
    if(!is.na(c[i])) {
      disc[i,j] <- floor(nstate * (c[i] - min) / (max - min))
    }
  }
}

write.csv(disc, file="cont.csv")
