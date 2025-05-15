library(vioplot)

setwd("~/Documents/Research/2025-ContChars/")

## distance metrics ##
nc_d  <- read.table("nonclock/discrete/dist_t.txt",        header=T)
nc_b  <- read.table("nonclock/both/dist_t.txt",            header=T)
nc_c  <- read.table("nonclock/continuous/dist_t.txt",      header=T)
nc_cm <- read.table("nonclock/continuous-mis/dist_t.txt",  header=T)

td_d   <- read.table("tipdating/discrete/dist_t.txt",      header=T)
td_b   <- read.table("tipdating/both/dist_t.txt",          header=T)
td_c  <- read.table("tipdating/continuous/dist_t.txt",     header=T)
td_cm  <- read.table("tipdating/continuous-mis/dist_t.txt",header=T)

# 1. plot Quartet metrics
par(mfrow=c(2,2), mar=c(2.5,4.5,2.5,0.5))

vioplot(nc_d$d_qt, nc_b$d_qt, nc_c$d_qt, nc_cm$d_qt,
        xaxt="n", ylim=c(0,0.25))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(main="distance metrics (nonclock)", ylab="Quartet")

vioplot(td_d$d_qt, td_b$d_qt, td_c$d_qt, td_cm$d_qt,
        xaxt="n", ylim=c(0,0.25))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(main="distance metrics (tip-dating)")

# and MCI distance metrics
vioplot(nc_d$d_ci, nc_b$d_ci, nc_c$d_ci, nc_cm$d_ci,
        xaxt="n", ylim=c(0,0.4))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))
title(ylab="MCI")

vioplot(td_d$d_ci, td_b$d_ci, td_c$d_ci, td_cm$d_ci,
        xaxt="n", ylim=c(0,0.4))
axis(side=1, at=1:4, labels=c("disc","both","cont","cont_m"))


# 2. plot resolution vs SJA
par(mfrow=c(4,2), mar=c(1,1,1,1))

#(new fig) 
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
