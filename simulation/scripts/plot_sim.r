library(vioplot)
library(stringr)

setwd("~/Documents/Research/2025-ContChars/simulation")

# calculate root mean square error
calc_rmse <- function(s_true, s_estm) {
  mse <- rep(0, 100)
  for (i in 1:100) {
    tmp <- str_extract_all(s_true[i],
                           "length_mean=\\d+\\.\\d+[eE][+-]?\\d+",
                           simplify=T)
    true <- as.numeric(str_split(tmp, "=", simplify=T)[,2])
    tmp <- str_extract_all(s_estm[i],
                           "length_mean=\\d+\\.\\d+[eE][+-]?\\d+",
                           simplify=T)
    estm <- as.numeric(str_split(tmp, "=", simplify=T)[,2])
    
    if (true[length(true)] == 0)  # for rooted tree
      len <- length(true) - 1
    else
      len <- length(true)
    
    for (j in 1:len)
      mse[i] <- ( (estm[j] - true[j])^2 + mse[i] * (j-1) ) / j
  }
  sqrt(mse)
}

# calculate mean HPD width
calc_mciw <- function(s_true, s_estm) {
  ciw <- rep(0, 100)
  for (i in 1:100) {
    tmp <- str_extract_all(s_true[i],
                           "length_mean=\\d+\\.\\d+[eE][+-]?\\d+",
                           simplify=T)
    true <- as.numeric(str_split(tmp, "=", simplify=T)[,2])
    tmp <- str_extract_all(s_estm[i], 
                         "length_95%HPD=\\{\\d+\\.\\d+[eE][+-]?\\d+,\\d+\\.\\d+[eE][+-]?\\d+",
                         simplify=T)
    hpd <- str_split(str_replace(tmp, "length_95%HPD=\\{", ""), ",", simplify=T)
    
    if (true[length(true)] == 0)  # for rooted tree
      len <- length(true) - 1
    else
      len <- length(true)
    
    for (j in 1:len)
      ciw[i] <- ( (as.numeric(hpd[j,2]) - as.numeric(hpd[j,1])) + ciw[i] * (j-1) ) / j
  }
  return(ciw)
}

## distance metrics ##
dist_nc_d  <- read.table("nonclock/discrete/dist_t.txt",       header=T)
dist_nc_b  <- read.table("nonclock/both/dist_t.txt",           header=T)
dist_nc_c  <- read.table("nonclock/continuous/dist_t.txt",     header=T)
dist_nc_cm <- read.table("nonclock/continuous-mis/dist_t.txt", header=T)

dist_td_d  <- read.table("tipdating/discrete/dist_t.txt",      header=T)
dist_td_b  <- read.table("tipdating/both/dist_t.txt",          header=T)
dist_td_c  <- read.table("tipdating/continuous/dist_t.txt",    header=T)
dist_td_cm <- read.table("tipdating/continuous-mis/dist_t.txt",header=T)

## branch lengths ##
tt_nc_d  <- readLines("nonclock/discrete/true.brl.log")
es_nc_d  <- readLines("nonclock/discrete/estm.brl.log")
tt_nc_b  <- readLines("nonclock/both/true.brl.log")
es_nc_b  <- readLines("nonclock/both/estm.brl.log")
tt_nc_c  <- readLines("nonclock/continuous/true.brl.log")
es_nc_c  <- readLines("nonclock/continuous/estm.brl.log")
tt_nc_cm <- readLines("nonclock/continuous-mis/true.brl.log")
es_nc_cm <- readLines("nonclock/continuous-mis/estm.brl.log")
rmse_nc_d  <- calc_rmse(tt_nc_d,  es_nc_d)
rmse_nc_b  <- calc_rmse(tt_nc_b,  es_nc_b)
rmse_nc_c  <- calc_rmse(tt_nc_c,  es_nc_c)
rmse_nc_cm <- calc_rmse(tt_nc_cm, es_nc_cm)
mciw_nc_d  <- calc_mciw(tt_nc_d,  es_nc_d) 
mciw_nc_b  <- calc_mciw(tt_nc_b,  es_nc_b)
mciw_nc_c  <- calc_mciw(tt_nc_c,  es_nc_c)
mciw_nc_cm <- calc_mciw(tt_nc_cm, es_nc_cm)

tt_td_d  <- readLines("tipdating/discrete/true.brl.log")
es_td_d  <- readLines("tipdating/discrete/estm.brl.log")
tt_td_b  <- readLines("tipdating/both/true.brl.log")
es_td_b  <- readLines("tipdating/both/estm.brl.log")
tt_td_c  <- readLines("tipdating/continuous/true.brl.log")
es_td_c  <- readLines("tipdating/continuous/estm.brl.log")
tt_td_cm <- readLines("tipdating/continuous-mis/true.brl.log")
es_td_cm <- readLines("tipdating/continuous-mis/estm.brl.log")
rmse_td_d  <- calc_rmse(tt_td_d,  es_td_d)
rmse_td_b  <- calc_rmse(tt_td_b,  es_td_b)
rmse_td_c  <- calc_rmse(tt_td_c,  es_td_c)
rmse_td_cm <- calc_rmse(tt_td_cm, es_td_cm)
mciw_td_d  <- calc_mciw(tt_td_d,  es_td_d) 
mciw_td_b  <- calc_mciw(tt_td_b,  es_td_b)
mciw_td_c  <- calc_mciw(tt_td_c,  es_td_c)
mciw_td_cm <- calc_mciw(tt_td_cm, es_td_cm)


par(mfrow=c(3,2), mar=c(2.5,4.5,2.5,0.5))

# 1. plot Quartet distances
vioplot(dist_nc_d$d_qt, dist_nc_b$d_qt, dist_nc_c$d_qt, dist_nc_cm$d_qt,
        xaxt="n", ylim=c(0,0.3))
title(main="nonclock", ylab="Quartet distance")

vioplot(dist_td_d$d_qt, dist_td_b$d_qt, dist_td_c$d_qt, dist_td_cm$d_qt,
        xaxt="n", ylim=c(0,0.3))
title(main="tip-dating")
legend("topright", "topology", bty="n")

# 2. plot RMSE of brls
vioplot(rmse_nc_d, rmse_nc_b, rmse_nc_c, rmse_nc_cm,
        xaxt="n", ylim=c(0,0.1))
title(ylab="RMSE")

vioplot(rmse_td_d, rmse_td_b, rmse_td_c, rmse_td_cm,
        xaxt="n", ylim=c(0,0.1))
legend("topright", "branch lengths", bty="n")

# 3. plot mean CI width of brls
vioplot(mciw_nc_d, mciw_nc_b, mciw_nc_c, mciw_nc_cm,
        xaxt="n", ylim=c(0,0.2))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(ylab="mean HPD width")

vioplot(mciw_td_d, mciw_td_b, mciw_td_c, mciw_td_cm,
        xaxt="n", ylim=c(0,0.2))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
legend("topright", "branch lengths", bty="n")


# plot resolution vs SJA (fig s1)
par(mfrow=c(4,2), mar=c(1,2,1,1))

plot(nc_d$sja ~nc_d$reso,  xlim=c(0.5,1), ylim=c(0.5,1));
legend("bottomleft", "discrete", bty="n")
plot(td_d$sja ~td_d$reso,  xlim=c(0.5,1), ylim=c(0.5,1));
legend("bottomleft", "discrete", bty="n")
plot(nc_b$sja ~nc_b$reso,  xlim=c(0.5,1), ylim=c(0.5,1));
legend("bottomleft", "both", bty="n")
plot(td_b$sja ~td_b$reso,  xlim=c(0.5,1), ylim=c(0.5,1));
legend("bottomleft", "both", bty="n")
plot(nc_c$sja ~nc_c$reso,  xlim=c(0.5,1), ylim=c(0.5,1));
legend("bottomleft", "continuous", bty="n")
plot(td_c$sja ~td_c$reso,  xlim=c(0.5,1), ylim=c(0.5,1));
legend("bottomleft", "continuous", bty="n")
plot(nc_cm$sja~nc_cm$reso, xlim=c(0.5,1), ylim=c(0.5,1));
legend("bottomleft",  "cont_m", bty="n")
legend("bottomright", "nonclock     ", bty="n")
plot(td_cm$sja~td_cm$reso, xlim=c(0.5,1), ylim=c(0.5,1));
legend("bottomleft",  "cont_m", bty="n")
legend("bottomright", "tip-dating   ", bty="n")


# the following are not included
## tree lengths ##
true  <- read.table("simulator/bd.tl.txt", header=T)
tl_d  <- read.table("nonclock/discrete/tl_estm.txt",       header=F)
tl_b  <- read.table("nonclock/both/tl_estm.txt",           header=F)
tl_c  <- read.table("nonclock/continuous/tl_estm.txt",     header=F)
tl_cm <- read.table("nonclock/continuous-mis/tl_estm.txt", header=F)
## tree heights ##
th_d  <- read.table("tipdating/discrete/th_estm.txt",      header=F)
th_b  <- read.table("tipdating/both/th_estm.txt",          header=F)
th_c  <- read.table("tipdating/continuous/th_estm.txt",    header=F)
th_cm <- read.table("tipdating/continuous-mis/th_estm.txt",header=F)
# relative bias
tl_bias_d  <- (tl_d$V2 - true$tLength) / true$tLength
tl_bias_b  <- (tl_b$V2 - true$tLength) / true$tLength
tl_bias_c  <- (tl_c$V2 - true$tLength) / true$tLength
tl_bias_cm <- (tl_cm$V2- true$tLength) / true$tLength
th_bias_d  <- (th_d$V2 - true$tHeight) / true$tHeight
th_bias_b  <- (th_b$V2 - true$tHeight) / true$tHeight
th_bias_c  <- (th_c$V2 - true$tHeight) / true$tHeight
th_bias_cm <- (th_cm$V2- true$tHeight) / true$tHeight
# relative CI width
tl_CIw_d  <- (tl_d$V5 - tl_d$V4) / true$tLength
tl_CIw_b  <- (tl_b$V5 - tl_b$V4) / true$tLength
tl_CIw_c  <- (tl_c$V5 - tl_c$V4) / true$tLength
tl_CIw_cm <- (tl_cm$V5-tl_cm$V4) / true$tLength
th_CIw_d  <- (th_d$V5 - th_d$V4) / true$tHeight
th_CIw_b  <- (th_b$V5 - th_b$V4) / true$tHeight
th_CIw_c  <- (th_c$V5 - th_c$V4) / true$tHeight
th_CIw_cm <- (th_cm$V5-th_cm$V4) / true$tHeight

# plot tree lengths/heights
vioplot(tl_bias_d, tl_bias_b, tl_bias_c, tl_bias_cm,
        xaxt="n", ylim=c(-0.2,0.2))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(ylab="relative bias")
legend("topright", "tree length", bty="n")

vioplot(th_bias_d, th_bias_b, th_bias_c, th_bias_cm,
        xaxt="n", ylim=c(-0.2,0.2))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
legend("topright", "tree height", bty="n")

vioplot(tl_CIw_d, tl_CIw_b, tl_CIw_c, tl_CIw_cm,
        xaxt="n", ylim=c(0,0.8))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(ylab="relative CI width")
legend("topright", "tree length", bty="n")

vioplot(th_CIw_d, th_CIw_b, th_CIw_c, th_CIw_cm,
        xaxt="n", ylim=c(0,0.8))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
legend("topright", "tree height", bty="n")
