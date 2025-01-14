#####################################
## au_sable_ts_plot.R by JPJ 6 i 25
#####################################

## PURPOSE: to plot the au sable GPP time series
## USAGE: Rscript au_sable_ts_plot.R

aus <- read.csv("au_sable_river_prepped.csv", header=TRUE)
pdf("au_sable_ts_plot.pdf", height=3, width=9)
par(mar=c(5,5,1,1))
plot(1:dim(aus)[1], aus[,6], type="n", xlab="Year", ylab="GPP", cex.lab=1.75, cex.axis=1.5, xaxt="n", las=1)
axis(1, at=c(mean(c(1,366)), mean(c(367,731)), mean(c(732,1096))), labels=c(2012, 2013, 2014), cex.axis=1.5)
rect(366.5, -6, 731.5, 6, col="grey", lty=0) ## shading the middle year grey
points(1:366, aus[aus[,4]==2012,6], type="l")
points(367:731, aus[aus[,4]==2013,6], type="l")
points(732:1096, aus[aus[,4]==2014,6], type="l")

## marking missing data with red lines
	## each data point has a decimal, so you can find missing data by not matching a "."
for (i in 1:dim(aus)[1]){
	if (grepl(".", aus[i,6])=="FALSE") { segments(i, par("usr")[3], i, par("usr")[3]+0.6, col="red", lwd=1.5) }
}

box(lwd=1)
dev.off()
