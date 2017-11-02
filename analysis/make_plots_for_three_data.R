#####load packages and functions#####
library(lme4)
library(dplyr)
library(lmvar)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(CopyNumberCellShift)
setwd('/Users/tae/Dropbox/TaeProject/CopyNumber/CopyNumberCellShift')
Rcpp::sourceCpp('src/all_functions.cpp')

#####set Y ######
j =  22 #chromosome ID
load("../summary.100kb.normalized.RData")
X = data.100kb[data.100kb[, 1]%in%paste0("chr", j), -c(1:3)]
orig_Y = as.matrix(apply(X, 2, as.numeric))

######polygenomic breast tumor#########
orig = read.table('../polygenomic_breast_tumor_copynumber.txt', header = TRUE, stringsAsFactors = FALSE)
orig = orig[orig$CHR=='chr22',-(1:3)]
Y_center = orig - 2
diffvec = apply(Y_center, 2, function(x) sum(x^2)/nrow(Y_center))
ind = which(diffvec < 0.01)

orig2 = read.table('../polygenomic_breast_tumor_counts.txt', header=TRUE, stringsAsFactors = FALSE)
orig2 = orig2[orig2$CHR=='chr22',-(1:3)]

vec = rep(0,length(ind))
for (j in 1:length(ind)){
  i = ind[j]
  temp = orig2[,i]
  cutoff = quantile(temp, c(0.25, 0.75))
  temp2 = temp[temp > cutoff[1] & temp < cutoff[2]]
  vec[j] = (median(temp2))
}
nor = mean(vec)
poly_Y = log(orig2 / nor)

######circulating lung tumor cells#########
orig = read.table('../circulating_lung_tumor_copynumber.txt', header = TRUE, stringsAsFactors = FALSE)
orig = orig[orig$CHR=='chr22',-(1:3)]
Y_center = orig - 2
diffvec = apply(Y_center, 2, function(x) sum(x^2)/nrow(Y_center))
ind = which(diffvec < 0.01)

orig2 = read.table('../circulating_lung_tumor_counts.txt', header=TRUE, stringsAsFactors = FALSE)
orig2 = orig2[orig2$CHR=='chr22',-(1:3)]

vec = rep(0,length(ind))
for (j in 1:length(ind)){
  i = ind[j]
  temp = orig2[,i]
  cutoff = quantile(temp, c(0.25, 0.75))
  temp2 = temp[temp >= cutoff[1] & temp <= cutoff[2]]
  vec[j] = (median(temp2))
}
nor = mean(vec)
lung_Y = log(orig2 / nor)




organize_by_cluster = function(Y){
  ######get clusters and organize Y#####
  clust <- hclust(dist(t(Y)))
  plot(clust)
  clusterCut <- cutree(clust, 2)
  ind = sort(clusterCut, index.return = TRUE)$ix
  Y = Y[,ind]
  Y = cbind(Y[,clusterCut==1], Y[, clusterCut==2], Y[,clusterCut==3], Y[,clusterCut==4], Y[,clusterCut==5])
  return(Y)
}
orig_Y = organize_by_cluster(orig_Y)
lung_Y = organize_by_cluster(lung_Y)
poly_Y = organize_by_cluster(poly_Y)

get_result = function(Y){
  wts = as.numeric(defaultWeights_c(nrow(Y)))
  steps = min(nrow(Y),ncol(Y))-8
  #res_mba = MBAmethyl(as.matrix(Y), wts, steps)
  res = cnv_c(as.matrix(Y), wts, steps, 30)
  res_old = cnv_c_old(as.matrix(Y),wts,steps,30)
  return(list(res=res, res_old = res_old))
}
orig_li = get_result(orig_Y)
lung_li = get_result(lung_Y)
poly_li = get_result(poly_Y)
orig_res = orig_li$res; orig_res_old = orig_li$res_old
lung_res = lung_li$res; lung_res_old = lung_li$res_old
poly_res = poly_li$res; poly_res_old = poly_li$res_old



######only theta#####
only_theta_plot = function(res, res_old, Y, label){
  i=1
  onlytheta = data.frame(x = 1:nrow(Y),
                         CNVjoint=res$aic$theta[,i],
                         MBAmethyl=res_old$aic$theta[,i])
  onlytheta2=melt(onlytheta,id.vars='x')
  onlytheta2$sample = i
  final=onlytheta2
  for (i in 2:ncol(Y)){
    onlytheta = data.frame(x = 1:nrow(Y),
                           CNVjoint=res$aic$theta[,i],
                           MBAmethyl=res_old$aic$theta[,i])
    onlytheta2=melt(onlytheta,id.vars='x')
    onlytheta2$sample = i
    final=rbind(final,onlytheta2)
  }
  samplepick=c(1:4, (ncol(Y)-3) : ncol(Y))
  g = ggplot(final[final$sample%in%(samplepick),], aes(x=x, y=value,col=variable))+
    geom_line(size=0.5, alpha=1)+
    facet_wrap(~sample,nrow=4)+
    ggtitle(label)+
    guides(colour = guide_legend(title=''))
  return(g)
}
orig_g = only_theta_plot(orig_res, orig_res_old, orig_Y, 'original')
lung_g = only_theta_plot(lung_res, lung_res_old, lung_Y, 'MALBAC')
poly_g = only_theta_plot(poly_res, poly_res_old, poly_Y, 'DOP-PCR')

pdf('onlytheta.pdf', width=12, height=8)
grid.arrange(orig_g, lung_g, poly_g, ncol=3)
dev.off()

#####compute distance matrix######
get_distance_matrix = function(cluster, res, res_old){
  diffmat_new = matrix(0, length(cluster), length(cluster))
  diffmat_old = matrix(0, length(cluster), length(cluster))
  for (i in 1:(length(cluster)-1)){
    for (j in (i+1):length(cluster)){
      vec1 = res$aic$theta[,cluster[i]]
      vec2 = res$aic$theta[,cluster[j]]
      diffmat_new[i,j] = diffmat_new[j,i] = sum((vec1-vec2)^2)/length(vec1)
      vec1 = res_old$aic$theta[,i]
      vec2 = res_old$aic$theta[,j]
      diffmat_old[i,j] = diffmat_old[j,i] = sum((vec1-vec2)^2)/length(vec1)
    }
  }
  return(list(
    new = sum(diffmat_new)/(length(cluster)*(length(cluster)-1)),
    old = sum(diffmat_old)/(length(cluster)*(length(cluster)-1)),
    old_diffmat = diffmat_old,
    new_diffmat = diffmat_new))
}


visualize_distance = function(Y, res, res_old, label){
  across_cluster = 1:ncol(Y)
  out = get_distance_matrix(across_cluster, res, res_old)

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
    labs(title=label)

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
  g = grid.arrange(odplot, ndplot, nrow=1)
  return(g)
}

orig_dist = get_distance_matrix(1:ncol(orig_Y), orig_res, orig_res_old)
orig_dist_normal = get_distance_matrix(1:21, orig_res, orig_res_old)
orig_dist_normal$old
orig_dist_normal$new
orig_dist_tumor = get_distance_matrix(22:32, orig_res, orig_res_old)
orig_dist_tumor$old
orig_dist_tumor$new

g1 = visualize_distance(orig_Y, orig_res, orig_res_old, 'original')
g2 = visualize_distance(lung_Y, lung_res, lung_res_old, 'MALBAC')
g3 = visualize_distance(poly_Y, poly_res, poly_res_old, 'DOP-PCR')
gfinal = grid.arrange(g1,g2,g3, nrow=3)
pdf("diffmat.pdf", width = 8, height = 12) # Open a new pdf file
grid.arrange(g1,g2,g3, nrow=3) # Write the grid.arrange in the file
dev.off()



