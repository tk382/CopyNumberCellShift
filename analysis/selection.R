#####load packages and functions#####
library(lme4)
library(dplyr)
library(lmvar)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(CopyNumberCellShift)
setwd('/Users/tae/Dropbox/TaeProject/CopyNumber/CopyNumberCellShift')
Rcpp::sourceCpp('src/cnvJoint.cpp')

#####set Y ######
j =  22 #chromosome ID
load("../summary.100kb.normalized.RData")
X = data.100kb[data.100kb[, 1]%in%paste0("chr", j), -c(1:3)]
Y = as.matrix(apply(X, 2, as.numeric))

wts = as.numeric(defaultWeights_c(nrow(Y)))
steps = min(nrow(Y),ncol(Y))-1
res = cnvJoint(Y = as.matrix(Y), wts = wts, steps = steps, maxloop = 30, verbose=TRUE)

for (i in 1:32){
  plot(Y[,i], type = 'l', ylim = c(-5,5), main=i)
  lines(res$thetalist[,i,15], col = 'red', lwd=2)
  Sys.sleep(0.5)
}

for (i in 1:30){
  plot(Y[,31], type = 'l', ylim = c(-5,5), main=i)
  lines(res$thetalist[,31,i], col = 'red', lwd=2)
  #points(res$philist[,i], pch = 'o', col = 'blue')
  Sys.sleep(0.5)
}


##chr22 result##
#starts showing peak around 60 starting at k=13
#starts showing peak around 110 starting at 16

##chr21 result##
#starts showing changepoint around 10 at 8
#starts showing the second peak around 10 at 11
#starts showing changepoint around 360 at 23
#starts showing changepoint around 40 at 26

##chr20 result##
#aic picks 13
#bic picks 13
#starts showing changepoint around 300 at 5
#starts showing double-changepoint at around 260 at 8
