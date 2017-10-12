#####load packages and functions#####
library(lme4)
library(dplyr)
library(lmvar)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(MBAmethyl)
setwd('/Users/tae/Dropbox/TaeProject/CopyNumber/CopyNumberCellShift')
Rcpp::sourceCpp('src/all_functions.cpp')

#####set Y ######
j =  22 #chromosome ID
load("../summary.100kb.normalized.RData")
X = data.100kb[data.100kb[, 1]%in%paste0("chr", j), -c(1:3)]
Y = as.matrix(apply(X, 2, as.numeric))


######get clusters and organize Y#####
#clust <- hclust(dist(t(Y)))
#plot(clust)
#clusterCut <- cutree(clust, 3)
k = kmeans(t(Y), 3)
clusterCut = k$cluster
ind = sort(clusterCut, index.return = TRUE)$ix
Y = Y[,ind]
#Y = cbind(Y[,clusterCut==1], Y[, clusterCut==2], Y[,clusterCut==3])

######get wts, steps, and run gfLars#####
wts = as.numeric(defaultWeights_c(nrow(Y)))
steps = min(nrow(Y),ncol(Y))-1
#res_mba = MBAmethyl(as.matrix(Y), wts, steps)
res = cnv_c(as.matrix(Y), wts, steps, 30)
res_old = cnv_c_old(as.matrix(Y),wts,steps,30)

#####comparing error#####
ind = 1:steps  #change this to zoom in - maybe 15 to 30
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


######thetaphi#####
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
  geom_point(alpha = 0.3)+ylim(-4,4)+ggtitle(paste0("sample", picksample))

#abnormal
picksample=5
pick_resid = resid[resid$sample%in% (picksample), ]
plot_25=ggplot(pick_resid, aes(x=x, y=value, col=variable))+
  geom_point(alpha = 0.3)+ylim(-4,4)+ggtitle(paste("sample",picksample))

grid.arrange(plot_1, plot_25)


#####Y-xiphi######
newmat = Y-res$aic$phi %*% res$aic$xi
newmat2 = melt(newmat)
newmat2$Var2 = as.factor(rep(1:ncol(Y), each=nrow(Y)))
picksample = c(1, 90)
picked_newmat = newmat2[newmat2$Var2 %in% (picksample), ]
yminusxiphi=ggplot(picked_newmat, aes(x=variable, y=value, col=Var2))+
  geom_point(alpha = 0.9, size=0.5)+
  ylim(-5,5)+
  ggtitle(expression(Y[ij]~-~xi[j]*phi[i]~"for j=3 and 30"))+
  xlab("probes")

newmat3 = melt(Y-newmat)
newmat3$Var2 = as.factor(rep(1:ncol(Y), each=nrow(Y)))
picksample=c(1, 90)
picked_newmat = newmat3[newmat3$Var2 %in% (picksample), ]
xiphi=ggplot(picked_newmat, aes(x=variable, y=value, col=Var2))+
  geom_point(alpha = 0.9, size=0.5)+
  ylim(-5,5)+
  ggtitle(expression(xi[j]*phi[i]~"for j= 3 and 30"))+
  xlab("probes")

grid.arrange(yminusxiphi, xiphi)

#####show thetaphi and xiphi side by side#######
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




########only theta##########
i=1
onlytheta = data.frame(x = 1:nrow(Y), thetaplusxi = res$aic$theta[,i]+res$aic$xi[i],
                       theta=res$aic$theta[,i],
                       oldtheta=res_old$aic$theta[,i])
onlytheta2=melt(onlytheta,id.vars='x')
onlytheta2$sample = i
final=onlytheta2
for (i in 2:ncol(Y)){
  onlytheta = data.frame(x = 1:nrow(Y), thetaplusxi = res$aic$theta[,i]+res$aic$xi[i],
                         theta=res$aic$theta[,i],
                         oldtheta=res_old$aic$theta[,i])
  onlytheta2=melt(onlytheta,id.vars='x')
  onlytheta2$sample = i
  final=rbind(final,onlytheta2)
}
samplepick=c(17:24)
ggplot(final[final$sample%in%(samplepick),], aes(x=x, y=value,col=variable))+
  geom_line(size=1, alpha=0.5)+
  facet_wrap(~sample,nrow=4)


#######distance matrix#####
get_distance_matrix = function(cluster){
  diffmat_new = matrix(0, length(cluster), length(cluster))
  diffmat_old = matrix(0, length(cluster), length(cluster))
  for (i in 1:(length(cluster)-1)){
    for (j in (i+1):length(cluster)){
      vec1 = res$aic$theta[,cluster[i]]
      vec2 = res$aic$theta[,cluster[j]]
      diffmat_new[i,j] = diffmat_new[j,i] = sum((vec1-vec2)^2)

      vec1 = res_old$aic$theta[,i]
      vec2 = res_old$aic$theta[,j]
      diffmat_old[i,j] = diffmat_old[j,i] = sum((vec1-vec2)^2)
    }
  }
  return(list(
    new = sum(diffmat_new)/(length(cluster)*(length(cluster)-1)),
         old = sum(diffmat_old)/(length(cluster)*(length(cluster)-1)),
          old_diffmat = diffmat_old,
    new_diffmat = diffmat_new))
}

#####compute distance matrix######
intra_cluster = which(clusterCut==1)
outtemp = get_distance_matrix(intra_cluster)

across_cluster = 1:ncol(Y)
out = get_intra_cluster_var(across_cluster)


#####visualize distance matrix#####
new_diffmat = out$new_diffmat
old_diffmat = out$old_diffmat
nd = melt(new_diffmat)
od = melt(old_diffmat)
odplot = ggplot(od, aes(Var1, Var2)) +
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient(low = "white", high = "indianred",
                      guide= guide_colorbar(title.hjust = 0.2, title=''),
                      limits = c(0,max(nd$value))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(title='MALBAC circulating lung tumor cells')
ndplot = ggplot(nd, aes(Var1, Var2)) +
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient(low = "white", high = "indianred",
                      guide= guide_colorbar(title.hjust = 0.2, title=''),
                      limits = c(0,max(nd$value)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(title='')
grid.arrange(odplot, ndplot, nrow=1)


