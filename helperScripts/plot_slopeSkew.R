#install.packages('ggplot2')
library(ggplot2)
setwd('~/DNARepOrigin')
inDat = read.table('DNARepOrigin_output.tab')
inDat = inDat[which(inDat[,2] == 262144),]

inDat = cbind(inDat,Rank=as.numeric(as.factor(inDat[,2])))

pdf("Single_skew.pdf", width=10, height=7)
ggplot(inDat) + geom_rect(aes(xmin=V1-V2/100, xmax=V1+V2/200,ymin=Rank-.5,ymax=Rank+.5, color=V3))
dev.off()

pdf("Single_slope.pdf", width=10, height=7)
ggplot(inDat) + geom_rect(aes(xmin=V1-V2/100, xmax=V1+V2/200,ymin=Rank-.5,ymax=Rank+.5, color=V4))
dev.off()
