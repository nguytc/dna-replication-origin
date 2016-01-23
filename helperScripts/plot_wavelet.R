setwd('/home/support/tnguyen/DNARepOrigin/waveletOutput')
library(ggplot2)

test = scan('0_details.txt', what='character', sep="\n")
test2 = lapply(test, function(x){return(unlist(strsplit(x,"\t")))})

plotArray = data.frame()
for(i in 1:length(test2)){
  SiteCount = length(test2[[i]])
  SiteIntervals = seq(0,1,by=1/SiteCount)
  plotArray= rbind(plotArray,data.frame(Rank=i, PosMin=SiteIntervals[-1*length(SiteIntervals)], PosMax=SiteIntervals[-1], Coeffs=as.numeric(test2[[i]])/sd(as.numeric(test2[[i]]))))
}

pdf("0_details.pdf",height=12, width=12)
ggplot(plotArray) + geom_rect(aes(xmin=PosMin,xmax=PosMax,ymin=Rank-1, ymax=Rank, fill=Coeffs)) + theme_bw()
dev.off()

#test = read.table('swt_details.txt',sep="\t")
#test2 = lapply(test, function(x){return(unlist(strsplit(x,"\t")))})
