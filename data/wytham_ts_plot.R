###################################
## wytham_ts_plot.R by JPJ 10 x 24
###################################

## PURPOSE: to plot the wytham time series
## USAGE: Rscript wytham_ts_plot.R

wytham <- read.csv(here("data/Wytham_tits.csv"), header=TRUE)
pdf("wytham_ts_plot.pdf", height=3, width=9)
par(mar=c(5,5,1,1))
plot(wytham[,1], wytham[,2], xlab="Year", ylab="Number of broods", ylim=c(100,500), pch=21, bg="grey", type="o", cex=1.5, cex.lab=1.5, cex.axis=1.25, las=1)
box(lwd=1.5)
dev.off()



