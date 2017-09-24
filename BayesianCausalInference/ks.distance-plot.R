rm(list = ls())
setwd("~/Desktop/code-BayesianCausalInference/")
# read 
ks.test.cntl <- read.csv(file = "ks.test.csv")[, 2:6]
threshold <- read.csv(file = "threshold.csv")[, 2:6]

pdf("ks-distance-99.pdf", width = 7, height = 5)
time <- 1:20
dates <- as.Date(time, origin = "2016-03-21")
dates_ten <- dates[seq(1, 20, by = 1)]
par(mfrow=c(2,3), mar=c(3,2.5,2.5,1), mgp=c(1.6,.6,0), oma=c(0,1.5,0,0),
    mai = c(0.3, 0.4, 0.3, 0.1))	
layout(mat = matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), nrow = 2, byrow = TRUE))
# layout.show(n = 5)
# axis(side=2,at=seq(0,8,1),las=2)
plot(dates, threshold[, 1], type = "o", ylim = c(0, 1), col = "lightblue", 
     ylab = "K-S distance", xaxt='n', xlab = " ")
lines(dates, ks.test.cntl[, 1], col = "red", type = "o")
title("Dataset 1")
axis(side=1, at=dates_ten, labels=format(dates_ten, "%b-%d-%y"), las = 1,
     cex.axis=0.8, las = 1, font = 2, tcl = -0.2, padj = -1)
plot(dates, threshold[, 2], type = "o", ylim = c(0, 1), col = "lightblue", 
     ylab = "K-S distance", xaxt='n', xlab = " ")
title("Dataset 2")
lines(dates, ks.test.cntl[, 2], col = "red", type = "o")
axis(side=1, at=dates_ten, labels=format(dates_ten, "%b-%d-%y"), las = 1,
     cex.axis=0.8, las = 1, font = 2, tcl = -0.2, padj = -1)
plot(dates, threshold[, 3], type = "o", ylim = c(0, 1), col = "lightblue", ylab = "K-S distance",
     xaxt='n', xlab = " ")
title("Dataset 3")
lines(dates, ks.test.cntl[, 3], col = "red", type = "o")
axis(side=1, at=dates_ten, labels=format(dates_ten, "%b-%d-%y"), las = 1,
     cex.axis=0.8, las = 1, font = 2, tcl = -0.2, padj = -1)
plot(dates, threshold[, 4], type = "o", ylim = c(0, 1), col = "lightblue", ylab = "K-S distance",
     xaxt='n', xlab = " ")
title("Dataset 4")
lines(dates, ks.test.cntl[, 4], col = "red", type = "o")
axis(side=1, at=dates_ten, labels=format(dates_ten, "%b-%d-%y"), las = 1,
     cex.axis=0.8, las = 1, font = 2, tcl = -0.2, padj = -1)
plot(dates, threshold[, 5], type = "o", ylim = c(0, 1), col = "lightblue", ylab = "K-S distance",
     xaxt='n', xlab = " ")
title("Dataset 5")
lines(dates, ks.test.cntl[, 5], col = "red", type = "o")
axis(side=1, at=dates_ten, labels=format(dates_ten, "%b-%d-%y"), las = 1,
     cex.axis=0.8, las = 1, font = 2, tcl = -0.2, padj = -1)
dev.off()
