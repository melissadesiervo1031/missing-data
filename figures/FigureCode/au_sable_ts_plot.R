#####################################
## au_sable_ts_plot.R by JPJ 6 i 25
#####################################
library(tidyverse)
## PURPOSE: to plot the au sable GPP time series
## USAGE: Rscript au_sable_ts_plot.R

aus <- read.csv("./data/au_sable_river_prepped.csv", header=TRUE)
ausNew <- aus %>% 
  mutate(date = as.POSIXct(date)) 
ggplot() + 
  geom_rect(aes(xmin = as.POSIXct("2014-01-01T00:00:00Z"), xmax = as.POSIXct("2014-12-31T00:00:00Z"), 
                ymin = -6, ymax = 6), fill = "grey80") +
  geom_line(data = ausNew, aes(x = date, y = GPP)) + 
  geom_rug(data = ausNew[is.na(ausNew$GPP),], aes(date), col = "red") +
  theme_classic() + 
  labs(x = "Year", y = "GPP (scaled)")

pdf("./figures/au_sable_ts_plot.pdf", height=3, width=9)
par(mar=c(5,7,1,1))
plot(1:dim(aus)[1], aus[,6], type="n", xlab="Year", ylab="", cex.lab=1.75, cex.axis=1.5, xaxt="n", las=1)
mtext("Scaled GPP", side=2, line=5, cex=1.75)
mtext(expression("(g O"[2]~" m"^-2~"d"^-1~")"), side=2, line=3, cex=1.25)
axis(1, at=c(mean(c(1,366)), mean(c(367,731)), mean(c(732,1096))), labels=c(2012, 2013, 2014), cex.axis=1.5)
rect(731.5, -6, 1096, 6, col="grey", lty=0) ## shading the middle year grey
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



