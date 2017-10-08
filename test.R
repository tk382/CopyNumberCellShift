library(lme4)
library(dplyr)
library(lmvar)
library(ggplot2)
library(reshape2)
library(gridExtra)
setwd('/Users/tae/Dropbox/TaeProject/CopyNumber/cnvJoint')
Rcpp::sourceCpp('src/all_functions.cpp')


#Run R functions separately - simdata2
li = simdata2(382, 32)
Y = as.matrix(li$Y)
for (j in 1:32){
  plot((li$theta[,j] + li$xi[j])*li$phi, ylim = c(-6,6))
  points((res$aic$theta[,j] + res$aic$xi[j]) * res$aic$phi, col='red')
  Sys.sleep(0.5)
}
plot(li$theta[,1])
points(res$aic$theta[,1], col = 'red')

plot(li$phi, type = 'l')
lines(-res$aic$phi, col = 'red')

plot(li$xi, type = 'l')
lines(res$aic$xi[1,], col = 'red')

j =  22 #chromosome ID
load("../summary.100kb.normalized.RData")
X = data.100kb[data.100kb[, 1]%in%paste0("chr", j), -c(1:3)]
Y = as.matrix(apply(X, 2, as.numeric))
wts = defaultWeights_c(nrow(Y))
steps = ncol(Y)-1

res = cnv_c(Y, wts, steps, 30)
res_old = cnv_c_old(Y,wts,steps,30)

plot(res$aic$xi[1,], ylim=c(-1,1),
     ylab='xi',
     xlab='samples')
points(c(3,5,12,27), res$aic$xi[1,c(3,5,12,27)],
       col='red')

###########comparing error
ind = 1:32
aicdat = data.frame(x = ind, aic_after = res$aicerror[ind], aic_before=res_old$aicerror[ind])
aicdat2 = melt(aicdat, id.vars='x')
aicplot = ggplot(aicdat2, aes(x=x, y=value, col=variable))+
  geom_line() + ggtitle('aic')+theme(legend.position="none") + labs(x='changepoints',y='aic')

bicdat = data.frame(x = ind, bic_after = res$bicerror[ind], bic_before=res_old$bicerror[ind])
bicdat2 = melt(bicdat, id.vars='x')
bicplot = ggplot(bicdat2, aes(x=x, y=value, col=variable, guide=FALSE))+
  geom_line() + ggtitle('bic') + theme(legend.position="none")+ labs(x='changepoints',y='bic')

rssdat = data.frame(x=ind, after = res$rss[ind], before = res_old$rss[ind])
rssdat2 = melt(rssdat, id.vars='x')
rssplot = ggplot(rssdat2, aes(x=x, y=value, col=variable))+
  geom_line() + ggtitle('rss')+ theme(legend.position="none")+ labs(x='changepoints',y='rss')

grid.arrange(aicplot, bicplot, rssplot, ncol=3)


#thetaphi
residuals = data.frame(x=1:nrow(Y),
                       withxi=(res$aic$theta[,1]) * res$aic$phi,
                     withoutxi=res_old$aic$theta[,1]*res_old$aic$phi

)
for (i in 2:32){
  residuals = rbind(residuals, data.frame(x=1:nrow(Y),
                                          withxi=(res$aic$theta[,i]) * res$aic$phi,
                      withoutxi=res_old$aic$theta[,i]*res_old$aic$phi

  ))
}
resid = melt(residuals, id.vars='x')
resid$sample = rep(1:32, each=nrow(Y))
picksample = c(1:4, 29:32)
pick_resid = resid[resid$sample%in% (picksample), ]
ggplot(pick_resid, aes(x=x, y=value, col=variable))+
  geom_line(alpha = 0.7)+
  facet_wrap(~sample, nrow=4)+
  ylim(-5,5)+
  ylab(expression(Theta*phi))+
  xlab("probes")

#zoom in to one sample
#normal
picksample=3
pick_resid = resid[resid$sample%in% (picksample), ]
plot_1=ggplot(pick_resid, aes(x=x, y=value, col=variable))+
  geom_point(alpha = 0.3)+ylim(-4,4)+ggtitle("sample3")

#abnormal
picksample=30
pick_resid = resid[resid$sample%in% (picksample), ]
plot_25=ggplot(pick_resid, aes(x=x, y=value, col=variable))+
  geom_point(alpha = 0.3)+ylim(-4,4)+ggtitle("sample30")

grid.arrange(plot_1, plot_25)


####Y-xiphi
newmat = Y-res$aic$phi %*% res$aic$xi
newmat2 = melt(newmat)
newmat2$Var2 = as.factor(rep(1:32, each=nrow(Y)))
picksample = c(3, 30)
picked_newmat = newmat2[newmat2$Var2 %in% (picksample), ]
yminusxiphi=ggplot(picked_newmat, aes(x=Var1, y=value, col=Var2))+
  geom_point(alpha = 0.9, size=0.5)+
  ylim(-5,5)+
  ggtitle(expression(Y[ij]~-~xi[j]*phi[i]~"for j=3 and 30"))+
  xlab("probes")

newmat3 = melt(Y-newmat)
newmat3$Var2 = as.factor(rep(1:32, each=nrow(Y)))
picksample=c(3, 30)
picked_newmat = newmat3[newmat3$Var2 %in% (picksample), ]
xiphi=ggplot(picked_newmat, aes(x=Var1, y=value, col=Var2))+
  geom_point(alpha = 0.9, size=0.5)+
  ylim(-5,5)+
  ggtitle(expression(xi[j]*phi[i]~"for j= 3 and 30"))+
  xlab("probes")

grid.arrange(yminusxiphi, xiphi)

####show thetaphi and xiphi side by side######
bigper = data.frame(x=1:nrow(Y), y=res$aic$phi * res$aic$xi[1])
for (i in 2:32){
  perturbed = data.frame(x=1:nrow(Y), y=res$aic$phi * res$aic$xi[i])
  bigper = rbind(bigper, perturbed)
}
bigper2 = melt(bigper, id.vars='x')
bigper2$sample = as.factor(rep(1:32, each=nrow(Y)))
temp = bigper2[bigper2$sample%in%(c(1,3,30,31)),]
xiphi =ggplot(temp, aes(x=x, y=value, col=variable))+
  geom_point()+
  facet_wrap(~sample, nrow=4)+
  theme(legend.position="none")+
  ggtitle(expression(theta[ij]*phi[i])) +
  ylim(-6,6)


bigper = data.frame(x=1:nrow(Y), y=res$aic$theta[,1]*res$aic$phi)
for (i in 2:32){
  perturbed = data.frame(x=1:nrow(Y), y=res$aic$theta[,i]*res$aic$phi)
  bigper = rbind(bigper, perturbed)
}
bigper2 = melt(bigper, id.vars='x')
bigper2$sample = as.factor(rep(1:32, each=nrow(Y)))
temp2 = bigper2[bigper2$sample%in%(c(1,3,30,31)),]
thetaphi =ggplot(temp2, aes(x=x,y=value))+
  geom_point()+
  facet_wrap(~sample, nrow=4)+
  ggtitle(expression(xi[j]*phi[i]))+
  theme(legend.position="none")+
  ylim(-6,6)

grid.arrange(thetaphi, xiphi, ncol=2)




#compare theta alone
i=1
onlytheta = data.frame(x = 1:nrow(Y), thetaplusxi = res$aic$theta[,i]+res$aic$xi[i],
                       theta=res$aic$theta[,i],
                       oldtheta=res_old$aic$theta[,i])
onlytheta2=melt(onlytheta,id.vars='x')
onlytheta2$sample = i
final=onlytheta2
for (i in 2:32){
  onlytheta = data.frame(x = 1:nrow(Y), thetaplusxi = res$aic$theta[,i]+res$aic$xi[i],
                         theta=res$aic$theta[,i],
                         oldtheta=res_old$aic$theta[,i])
  onlytheta2=melt(onlytheta,id.vars='x')
  onlytheta2$sample = i
  final=rbind(final,onlytheta2)
}
samplepick=c(1:4, 29:32)
ggplot(final[final$sample%in%(samplepick),], aes(x=x, y=value,col=variable))+
  geom_line(size=1, alpha=0.5)+
  facet_wrap(~sample,nrow=4)


# #original data
# for (i in 1:32){
#   plot(Y[,i], ylim = c(-5,5))
#   abline(h=0, col = 'red', lwd=2)
#   Sys.sleep(0.5)
# }
#
# #create dataframe
# for (i in 1:32){
#   d = data.frame(x = 1:382,
#                  y = (res$aic$theta[,i] + res$aic$xi[i]) * res$aic$phi,
#                  y_old = res_old$aic$theta[,i]*res_old$aic$phi,
#                  truey = Y[,i])
#   d2 = melt(d, id.vars = 'x')
#   ggplot(d2, aes(x=x,y=value, col=variable))+geom_point()
# }
#
# rss1 = ic_c(k=length(res$aic$bkp),
#      Y, phi=res$aic$phi,
#      xi = res$aic$xi,
#      theta = res$aic$theta,
#      p = nrow(Y),
#      n = ncol(Y))
#
# rss2 = ic_c_old(k=length(res_old$aic$bkp),
#                 Y, phi=res_old$aic$phi,
#                 theta=res_old$aic$theta,
#                 p=nrow(Y), n=ncol(Y))
#
# #(theta+xi)*phi
# for (i in 1:32){
#   plot((res$aic$theta[,i] + res$aic$xi[i]) * res$aic$phi, ylim= c(-5,5))
#   abline(h=0, col = 'red', lwd=2)
#   points(Y[,i], col = 'red', pch='*')
#   Sys.sleep(0.5)
# }
#
# #theta only
# for (i in 28:32){
#   plot(res$aic$theta[,i]+res$aic$xi[i], ylim=c(-5,5))
#   abline(h=0, col = 'red', lwd=2)
#   Sys.sleep(0.5)
# }
#
#
# #Plot RSS
# plot(res$aicerror, ylim = c(-1.7,-1))
# plot(res$bicerror)




set.seed(1991)
n=10
p=100
cell_var = TRUE
probe_var=FALSE
#li = simulationdata(p,n,probe_var,cell_var)
li = simulationdata_ranef(p,n)
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
    res = cnv_c(Y, wts, steps, 10000)
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


cnvJoint_shift = function(Y, wts, steps, maxloop, maxiter){
  oldgam = rep(0, ncol(Y))
  oldphi = rep(0,nrow(Y))
  diff=1
  it=1
  n = ncol(Y); p = nrow(Y)
  while(diff>-1e-04 & it < maxiter){
    res = cnv_c(sweep(Y,2,oldgam), wts = wts, steps=steps, maxloop = 30)
    phi = res$aic$phi
    theta = res$aic$theta
    yhat = (as.numeric(phi)*theta)
    residual = Y-yhat
    #Ivec = 1:p;
    Jvec = 1:n
    d = data.frame(res = as.numeric(residual), j = rep(Jvec, each=p))
    mod = lmer(res ~ (1|j)-1, d)
    gam = (ranef(mod)$j)[,1]
    diff = sum((oldphi-phi)^2)
    plot(truetheta[,1], ylim = c(-3, 3))
    points(theta[,1], col = 'red')
    oldphi = phi
    oldgam = gam
    it = it+1
    print(diff)
  }
  return(list(res = res, ranef = gam, mod = mod))
}
s = sort(truephi, index.return=TRUE)
plot(s$x)
points(phi[s$ix], col = 'red')
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


plot(res$aic$theta[,1], type = 'l', ylim = c(-5,5))
for (i in 2:21){
 lines(res$aic$theta[,i])
}

plot(res_old$aic$theta[,1], type = 'l', ylim = c(-5,5))
for (i in 2:21){
  lines(res_old$aic$theta[,i])
}

plot(res$aic$theta[,22], type = 'l', ylim = c(-5,5))
for (i in 22:32){
  lines(res$aic$theta[,i])
}

plot(res_old$aic$theta[,22], type = 'l', ylim = c(-5,5))
for (i in 22:32){
  lines(res_old$aic$theta[,i])
}



get_intra_cluster_var = function(cluster){
  diffmat_new = matrix(0, length(cluster), length(cluster))
  diffmat_old = matrix(0, length(cluster), length(cluster))
  for (i in cluster[1:(length(cluster)-1)]){
    for (j in (i+1):cluster[length(cluster)]){
      vec1 = res$aic$theta[,i]
      vec2 = res$aic$theta[,j]
      diffmat_new[i,j] = diffmat_new[j,i] = sum((vec1-vec2)^2)

      vec1 = res_old$aic$theta[,i]
      vec2 = res_old$aic$theta[,j]
      diffmat_old[i,j] = diffmat_old[j,i] = sum((vec1-vec2)^2)
    }
  }
  return(list(
    new = sum(diffmat_new)/(length(cluster)*(length(cluster)-1)),
         old = sum(diffmat_old)/(length(cluster)*(length(cluster)-1))))
}
get_intra_cluster_var(cluster2)
