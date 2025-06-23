library(vioplot)

setwd("~/Documents/Research/2025-ContChars/simulation")

## distance metrics ##
nc_d  <- read.table("nonclock/discrete/dist_t.txt",        header=T)
nc_b  <- read.table("nonclock/both/dist_t.txt",            header=T)
nc_c  <- read.table("nonclock/continuous/dist_t.txt",      header=T)
nc_cm <- read.table("nonclock/continuous-mis/dist_t.txt",  header=T)

td_d   <- read.table("tipdating/discrete/dist_t.txt",      header=T)
td_b   <- read.table("tipdating/both/dist_t.txt",          header=T)
td_c   <- read.table("tipdating/continuous/dist_t.txt",    header=T)
td_cm  <- read.table("tipdating/continuous-mis/dist_t.txt",header=T)

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

# 0. plot Quartet metrics
vioplot(nc_d$d_qt, nc_b$d_qt, nc_c$d_qt, nc_cm$d_qt,
        xaxt="n", ylim=c(0,0.25))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(main="nonclock", ylab="Quartet")

vioplot(td_d$d_qt, td_b$d_qt, td_c$d_qt, td_cm$d_qt,
        xaxt="n", ylim=c(0,0.25))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(main="tip-dating")


par(mfrow=c(3,2), mar=c(2.5,4.5,2.5,0.5))

# 1. plot MCI distances
vioplot(nc_d$d_ci, nc_b$d_ci, nc_c$d_ci, nc_cm$d_ci,
        xaxt="n", ylim=c(0,0.4))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(main="nonclock", ylab="MCI distance")

vioplot(td_d$d_ci, td_b$d_ci, td_c$d_ci, td_cm$d_ci,
        xaxt="n", ylim=c(0,0.4))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(main="tip-dating")

# 2. plot tree lengths/heights
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
legend("topright", "tree length", bty="n")


# 3. plot resolution vs SJA (new fig)
par(mfrow=c(4,2), mar=c(1,1,1,1))

plot(nc_d$sja ~nc_d$reso,  xlim=c(0.5,1), ylim=c(0.6,1));
legend("bottomleft", "discrete", bty="n")
plot(td_d$sja ~td_d$reso,  xlim=c(0.5,1), ylim=c(0.6,1));
legend("bottomleft", "discrete", bty="n")
plot(nc_b$sja ~nc_b$reso,  xlim=c(0.5,1), ylim=c(0.6,1));
legend("bottomleft", "both", bty="n")
plot(td_b$sja ~td_b$reso,  xlim=c(0.5,1), ylim=c(0.6,1));
legend("bottomleft", "both", bty="n")
plot(nc_c$sja ~nc_c$reso,  xlim=c(0.5,1), ylim=c(0.6,1));
legend("bottomleft", "continuous", bty="n")
plot(td_c$sja ~td_c$reso,  xlim=c(0.5,1), ylim=c(0.6,1));
legend("bottomleft", "continuous", bty="n")
plot(nc_cm$sja~nc_cm$reso, xlim=c(0.5,1), ylim=c(0.6,1));
legend("bottomleft",  "cont_m", bty="n")
legend("bottomright", "nonclock     ", bty="n")
plot(td_cm$sja~td_cm$reso, xlim=c(0.5,1), ylim=c(0.6,1));
legend("bottomleft",  "cont_m", bty="n")
legend("bottomright", "tip-dating   ", bty="n")
