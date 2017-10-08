library(lme4)
library(dplyr)
library(lmvar)
library(ggplot2)
library(reshape2)
library(gridExtra)
setwd('/Users/tae/Dropbox/TaeProject/CopyNumber/CopyNumberCellShift')
Rcpp::sourceCpp('src/all_functions.cpp')


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
