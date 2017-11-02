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

#load data
load("../summary.100kb.normalized.RData")

get_result = function(chr){
  X = data.100kb[data.100kb[, 1]%in%paste0("chr", chr), -c(1:3)]
  Y = as.matrix(apply(X, 2, as.numeric))
  wts = as.numeric(defaultWeights_c(nrow(Y)))
  steps = min(nrow(Y),ncol(Y))-1
  res=1
  res = cnvJoint(Y = as.matrix(Y), wts = wts, steps = steps, maxloop = 30, verbose=TRUE)
  return(res)
}

res22 = get_result(22)
res21 = get_result(21)
res20 = get_result(20)


thetalist = res21$thetalist
philist   = res21$philist
xilist    = res21$xilist
bkplist   = res21$bkp
#first see the change in theta for each k
sumsq1 = sumsq2 = rep(0,29)
for (i in 2:30){
  theta1 = thetalist[,,i-1]
  theta2 = thetalist[,,i]
  sumsq2[i-1] = sum((colMeans(theta1)-colMeans(theta2))^2)
  sumsq1[i-1] = sum((theta1-theta2)^2)
}
plot(sumsq2)
plot(sumsq1)

#then see the likelihood
errorsq = rep(1,30); lik = rep(1,30)
for (k in 1:30){
  resid = Y-thetalist[,,k]*philist[,k] - philist[,k]%*%t(xilist[k,])
  #boxplot(as.numeric(resid), ylim = c(-10,10), main = i)
  #hist(as.numeric(resid), xlim = c(-9,9),main=i)
  errorsq[k] = sum(as.numeric(resid)^2)/(nrow(resid)*ncol(resid))
  lik[k] = sum(dnorm(as.numeric(resid), 0, sd(as.numeric(resid)), log = TRUE))
  #Sys.sleep(0.5)
}
plot(errorsq, xlab = 'number of changepoints k', ylab = 'RSS',
     main='ENCODE data chromosome 20')

#get and plot deltasum
deltasum = numeric(29)
for (k in 2:30){
  cp1 = bkplist[k-1,1:(k-1)]
  cp2 = bkplist[k, 1:k]
  comp = cp2[is.na(pmatch(cp2, cp1))]
  delta = thetalist[comp+1,,k] - thetalist[comp,,k]
  deltasum[k] = sum(as.numeric(delta)^2)
}
par(mar=c(6,6,4,1))
plot(deltasum, main = expression("change in"~Theta~"by k"),
     xlab = "k'th change point added",
     ylab = expression(Sigma[ij]~(Theta[ijk]-Theta[ij(k-1)])^2))
deltasumdf = data.frame(deltasum = deltasum, ind = 1:30)
ggplot(deltasumdf, aes(x=ind, y=deltasum))+geom_point()+
  xlab("k'th change point added") +
  ylab(expression(Sigma[ij]~(Theta[ijk]-Theta[ij(k-1)])^2~"for new change point i"))+
  ggtitle(expression("(a) change in "~Theta~"by k"))


k = 29 #select k
#which change point I should focus on?
cp1 = bkplist[k-1,1:(k-1)]
cp2 = bkplist[k, 1:k]
comp = cp2[is.na(pmatch(cp2, cp1))]
for (l in 1:length(comp)){
  plot(thetalist[comp[l],,k])
  points(thetalist[comp[l]+1,,k], col = 'red')
}
l=1
probeind = max((comp[l]-20),0):min(nrow(Y),(comp[l]+20))
plot(Y[probeind,10])
lines(thetalist[probeind, 10, k], col='red', lwd=2)
lines(thetalist[probeind, 10, k-1], col='blue', lwd=2)


#which sample to plot?
l=1
for (j in 1:30){
  plot(Y[probeind,j], main=j, ylim = c(-8,5))
  lines(thetalist[probeind, j, k], col='red', lwd=2)
  lines(thetalist[probeind, j, k-1], col='blue', lwd=2)
}

# I'll plot 17 and 21
make_sample_plot = function(l,j,k,plotlab){
  probeind = max((comp[l]-20),1):min(nrow(Y),(comp[l]+20))
  df = data.frame(x =  probeind,
                        Y = Y[probeind,j],
                        cp2 = thetalist[probeind,j,k],
                        cp1 = thetalist[probeind,j,k-1])
  df2 = data.frame(cp2 = thetalist[probeind,j,k],
                         cp1 = thetalist[probeind, j, k-1])
  df2 = melt(df2)
  lab2 = paste0('cp',k); lab1 = paste0('cp',k-1)
  df2$variable = rep(c(lab1, lab2), each=length(probeind))
  df2$x = rep(probeind, 2)
  gg = ggplot()+
    geom_point(data=df,aes(x=x, y=Y))+
    geom_line(data=df2, aes(x=x, y=value, col=variable))+
    ggtitle(paste0(plotlab,'sample ',j))+theme(legend.title=element_blank())+
    xlab('probe index')+ylab(expression(Theta))
  return(gg)
}
gg10 = make_sample_plot(1,10,29, '(c)')
gg4 = make_sample_plot(1,4,29, '(b)')
grid.arrange(gg4, gg10, nrow=1)


k=25
j=17

j=21
samp21df = data.frame(x =  (comp[l]-20):(comp[l]+20),
                      Y = Y[(comp[l]-20):(comp[l]+20),j],
                      cp25 = thetalist[(comp[l]-20):(comp[l]+20),j,k],
                      cp24 = thetalist[(comp[l]-20):(comp[l]+20), j, k-1])
samp21df2 = data.frame(cp25 = thetalist[(comp[l]-20):(comp[l]+20),j,k],
                       cp24 = thetalist[(comp[l]-20):(comp[l]+20), j, k-1])
samp21df2 = melt(samp21df2)
samp21df2$x = rep((comp[l]-20):(comp[l]+20), 2)
#samp17df2$x = rep(x = (comp[l]-20):(comp[l]+20),3)
gg21 = ggplot()+
  geom_point(data=samp21df,aes(x=x, y=Y))+
  geom_line(data=samp21df2, aes(x=x, y=value, col=variable))+
  ggtitle('(c) sample 21')+theme(legend.title=element_blank())+
  xlab('probe index')+ylab(expression(Theta))


grid.arrange(gg17, gg21, nrow=1)


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
