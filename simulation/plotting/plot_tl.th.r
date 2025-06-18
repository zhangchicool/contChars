library(vioplot)

setwd("~/Documents/Research/2025-ContChars/")

true  <- read.table("simulator/bd.tl.txt", header=T)

## tree lengths ##
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

# plot tree lengths and heights
par(mfrow=c(2,2), mar=c(2.5,4.5,2.5,0.5))

vioplot(tl_bias_d, tl_bias_b, tl_bias_c, tl_bias_cm,
        xaxt="n", ylim=c(-0.2,0.2))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(main="tree length (nonclock)", ylab="relative bias")

vioplot(th_bias_d, th_bias_b, th_bias_c, th_bias_cm,
        xaxt="n", ylim=c(-0.2,0.2))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(main="tree height (tip-dating)")

vioplot(tl_CIw_d, tl_CIw_b, tl_CIw_c, tl_CIw_cm,
        xaxt="n", ylim=c(0,0.8))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(ylab="relative CI width")

vioplot(th_CIw_d, th_CIw_b, th_CIw_c, th_CIw_cm,
        xaxt="n", ylim=c(0,0.8))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
