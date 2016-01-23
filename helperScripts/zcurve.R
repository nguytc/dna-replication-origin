library(ggplot2)

# MULTIPLE DATA SETS
data <- read.table('zcurve_xplusy.tab', header=FALSE)

print(data[0])

xaxis=(1:1000)/20
yaxis=(1:1000)/11

sm_obj = smooth.spline(x=xaxis, y=yaxis, spar=.6)

print(sm_obj)

x_val[which(diff(sign(predict(sm_obj, xaxis, deriv=1)$y)) != 0)]