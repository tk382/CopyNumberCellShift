library(lme4)
library(dplyr)
library(lmvar)
library(cnvJoint)
setwd('/Users/tae/Dropbox/TaeProject/CopyNumber/cnvJoint')
Rcpp::sourceCpp('src/all_functions.cpp')

j =  22# chromosome ID
load("../summary.100kb.normalized.RData")
X = data.100kb[data.100kb[, 1]%in%paste0("chr", j), -c(1:3)]
Y = as.matrix(apply(X, 2, as.numeric))
newY = Y[-c(87,25,56,382,24,6,58), ]
wts = defaultWeights_c(nrow(Y))
steps = ncol(Y)-1

newY = Y[-c(87,25,56,382,24,6,58), ]

n=10
p=100
cell_var = TRUE
probe_var=FALSE
li = simulationdata(p,n,probe_var,cell_var)
Y = li$Y
truecp = li$cp
truedelta = li$delta
truetheta = li$theta
truephi = li$phi
wts = defaultWeights_c(nrow(Y))
steps = ncol(Y)-1
for (i in 1:n){
  plot(Y[,i])
  Sys.sleep(0.2)
}
wts = defaultWeights_c(nrow(Y))
steps = ncol(Y)-1


cnvJoint2 = function(Y, wts, steps, maxloop, maxiter){
  n = ncol(Y); p = nrow(Y)
  delta = rep(1, n)
  oldphi = rep(1, p)
  it=1
  diff=10
  Jvec = 1:n
  Ivec = 1:p
  while(diff>1e-5 & it < maxiter){
    res = cnv_c(Y, delta, wts, steps)
    phi = res$bic$phi
    theta = res$bic$theta
    yhat = (as.numeric(phi)*theta)
    residual = Y-yhat
    d = data.frame(res = as.numeric(residual), j = rep(Jvec, each=p), i = rep(Ivec, n))
    dd = d %>% group_by(j) %>% summarize(sd = sd(res)*sqrt((p-1)/p))
    delta = dd$sd
    diff = sum((oldphi-phi)^2)
    oldphi = phi
    it = it+1
    plot(res$aic$theta[,1], ylim=c(-15, 15))
    points(truetheta[,1], col = 'red')
    plot(delta, ylim = c(0,10))
    points(truedelta, col = 'red')
    gc()
  }
}
cnvJoint2(Y,wts,steps,100,100)


#
# cnvJoint = function(Y, wts, steps, maxloop, maxiter){
#   oldgam = rep(0, ncol(Y))
#   oldphi = rep(0,nrow(Y))
#   diff=1
#   it=1
#   n = ncol(Y); p = nrow(Y)
#   while(diff>1e-04 & it < maxiter){
#     res = cnv_c(sweep(Y,2,oldgam), wts = wts, steps=steps, maxloop = 30)
#     phi = res$aic$phi
#     theta = res$aic$theta
#     yhat = (as.numeric(phi)*theta)
#     residual = Y-yhat
#     #Ivec = 1:p;
#     Jvec = 1:n
#     d = data.frame(res = as.numeric(residual), j = rep(Jvec, each=p))
#     #d = data.frame(res = as.numeric(residual), i = rep(Ivec, n), j = rep(Jvec, each=p))
#     mod = lmer(res ~ (1|j), d)
#     gam = (ranef(mod)$j)[,1]
#     diff = sum((oldphi-phi)^2)
#     oldphi = phi
#     oldgam = gam
#     it = it+1
#   }
#   return(list(res = res, ranef = gam, mod = mod))
# }

# res = cnv_c(Y, wts, steps, 30)
# out = cnvJoint(Y, wts, steps, 30, 50)
#
# #show that out is better than res
# plot(out$res$bicerror, type = 'l', lwd = 2, main = 'bic error')
# lines(res$bicerror, col = 'red', lwd=2)
# legend('topright', legend = c('new','old'), col = c('black', 'red'), pch = 1)
#
# plot(out$res$aicerror, type = 'l', lwd = 2, main = 'aic error')
# lines(res$aicerror, col = 'red', lwd = 2)
# legend('topright', legend = c('new','old'), col = c('black', 'red'), pch = 1)
#
#
# #do we really need non-constant variance for each probe?
# phi = as.numeric(out$res$aic$phi)
# theta = out$res$aic$theta
# yhat = phi*theta
# eps = rep(0, nrow(Y))
# for (i in 1:nrow(Y)){
#   eps[i] = var(Y[i,]-gam-yhat[i,])
# }
# plot(eps)
# #plot(eps, ylim = c(0,1))
#
#
# #difference in estimation of theta and phi
# plot(rowMeans(Y), ylim = c(-5, 5))
# points(rowMeans(res$aic$theta), col = 'blue')
# lines(res$aic$phi, col = 'red')
#
# plot(rowMeans(Y), ylim = c(-5, 5))
# points(rowMeans(theta), col = 'blue')
# lines(phi, col = 'red')
#
#
#
# plot(Y[,1]-rep(gam[1], nrow(Y)), ylim = c(-5, 5), type = 'l')
# for (i in 2:32){
#   lines(Y[,i]-rep(gam[i],nrow(Y)))
# }
# lines(rowMeans(theta), col = 'blue', lwd = 2)
# lines(rowMeans(res$aic$theta), col = 'red', lwd = 2)
#
