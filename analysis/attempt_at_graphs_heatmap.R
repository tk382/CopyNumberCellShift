library(lme4)
library(dplyr)
library(lmvar)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gplots)

setwd('/Users/tae/Dropbox/TaeProject/CopyNumber/CopyNumberCellShift')
Rcpp::sourceCpp('src/all_functions.cpp')

#####set Y ######
j =  22 #chromosome ID
load("../summary.100kb.normalized.RData")
X = data.100kb[data.100kb[, 1]%in%paste0("chr", j), -c(1:3)]
encode_Y = as.matrix(X)
encode_g = heatmap.2(encode_Y, tracecol = NA, Rowv = FALSE, Colv = TRUE, dendrogram = "column")
encode_Y = encode_Y[, encode_g$colInd]

lung_encode = read.table('../circulating_lung_tumor_counts.txt', header=TRUE, stringsAsFactors = FALSE)
lung = as.matrix(lung_encode[lung_encode$CHR=='chr22',-(1:3)])
lung_g = heatmap.2(lung, Rowv=FALSE, Colv = TRUE, dendrogram = "column", tracecol = NA)
lung_Y = lung[, lung_g$colInd]
lung_Y = log(lung_Y)

poly_encode = read.table('../polygenomic_breast_tumor_counts.txt', header = TRUE, stringsAsFactors = FALSE)
poly = as.matrix(poly_encode[poly_encode$CHR=='chr22', -(1:3)])
poly_g = heatmap.2(poly, tracecol = NA, Rowv=FALSE, Colv=TRUE, dendrogram = "column")
poly_Y = poly[, poly_g$colInd]
poly_Y = log(poly_Y)

get_result = function(Y){
  wts = as.numeric(defaultWeights_c(nrow(Y)))
  steps = min(nrow(Y),ncol(Y))-8
  res = cnv_c(as.matrix(Y), wts, steps, 30)
  res$aic$theta = res$aic$theta * as.numeric(sign(res$aic$phi))
  res$bic$theta = res$bic$theta * as.numeric(sign(res$bic$phi))
  res_MBAmethyl = cnv_c_old(as.matrix(Y),wts,steps,30)
  res_MBAmethyl$aic$theta = res_MBAmethyl$aic$theta * as.numeric(sign(res_MBAmethyl$aic$phi))
  res_MBAmethyl$bic$theta = res_MBAmethyl$bic$theta * as.numeric(sign(res_MBAmethyl$bic$phi))
  return(list(res=res, res_MBAmethyl = res_MBAmethyl))
}

encode_li = get_result(encode_Y)
lung_li = get_result(lung_Y)
poly_li = get_result(poly_Y)
encode_res = encode_li$res; encode_res_MBAmethyl = encode_li$res_MBAmethyl
lung_res = lung_li$res; lung_res_MBAmethyl = lung_li$res_MBAmethyl
poly_res = poly_li$res; poly_res_MBAmethyl = poly_li$res_MBAmethyl

encode_col = colnames(encode_Y)
colnames(encode_res$aic$theta) =
  colnames(encode_res$bic$theta) =
  colnames(encode_res_MBAmethyl$aic$theta) =
  colnames(encode_res_MBAmethyl$bic$theta) = encode_col

lung_col = colnames(lung_Y)
colnames(lung_res$aic$theta) =
  colnames(lung_res$bic$theta) =
  colnames(lung_res_MBAmethyl$aic$theta) =
  colnames(lung_res_MBAmethyl$bic$theta) = lung_col

poly_col = colnames(poly_Y)
colnames(poly_res$aic$theta) =
  colnames(poly_res$bic$theta) =
  colnames(poly_res_MBAmethyl$aic$theta) =
  colnames(poly_res_MBAmethyl$bic$theta) = poly_col

#######heatmap with clustering#########
pdf('analysis/writeup/paper/encode_cnvJoint_theta_heat.pdf', width=12, height=8)
encode_cnvJoint = heatmap.2(encode_res$aic$theta,
                     Colv=TRUE,
                     Rowv=FALSE,
                     dendrogram = "column",
                     main = 'ENCODE_cnvJoint',
                     tracecol = NA,
                     labRow=FALSE,
                     denscol="black",
                     key.title = NA,
                     keysize = 1
)
dev.off()
encode_ordered_cnvJoint = encode_res$aic$theta[,encode_cnvJoint$colInd ]
pdf('analysis/writeup/paper/encode_MBAmethyl_theta_heat.pdf', width=12, height=8)
encode_MBAmethyl = heatmap.2(encode_res_MBAmethyl$aic$theta,
          Colv=TRUE,
          Rowv=FALSE,
          dendrogram = "column",
          main = 'ENCODE_MBAmethyl',
          tracecol=NA,
          labRow=FALSE,
          denscol="black",
          key.title = NA,
          keysize = 1)
dev.off()
encode_ordered_MBAmethyl = encode_res_MBAmethyl$aic$theta[, encode_MBAmethyl$colInd]
pdf('analysis/writeup/paper/poly_cnvJoint_theta_heat.pdf', width=12, height=8)
poly_cnvJoint = heatmap.2(poly_res$aic$theta,
          Colv = TRUE,
          Rowv = FALSE,
          dendrogram = "column",
          main = 'DOP-PCR_cnvJoint',
          tracecol = NA,
          labRow=FALSE,
          denscol = "black",
          key.title = NA,
          keysize = 1)
dev.off()
poly_ordered_cnvJoint = poly_res$aic$theta[,poly_cnvJoint$colInd]
pdf('analysis/writeup/paper/poly_MBAmethyl_theta_heat.pdf', width=12, height=8)
poly_MBAmethyl = heatmap.2(poly_res_MBAmethyl$aic$theta,
          Colv = TRUE, Rowv = FALSE,
          dendrogram = "column",
          main = "DOP-PCR_MBAmethyl",
          tracecol = NA,
          labRow=FALSE,
          denscol = "black",
          key.title = NA,
          keysize = 1)
dev.off()
poly_ordered_MBAmethyl = poly_res_MBAmethyl$aic$theta[,poly_MBAmethyl$colInd]
pdf('analysis/writeup/paper/lung_cnvJoint_theta_heat.pdf', width=12, height=8)
lung_cnvJoint = heatmap.2(lung_res$aic$theta,
          Colv = TRUE, Rowv = FALSE,
          dendrogram = "column",
          main = "MALBAC_cnvJoint",
          tracecol = NA,
          labRow=FALSE,
          denscol="black",
          key.title = NA,
          keysize = 1)
dev.off()
lung_ordered_cnvJoint = lung_res$aic$theta[,lung_cnvJoint$colInd]
pdf('analysis/writeup/paper/lung_MBAmethyl_theta_heat.pdf', width=12, height=8)
lung_MBAmethyl = heatmap.2(lung_res_MBAmethyl$aic$theta,
          Colv = TRUE, Rowv = FALSE,
          dendrogram = "column",
          main = "MALBAC_MBAmethyl",
          labRow=FALSE,
          tracecol = NA,
          xlab = NA,
          denscol="black",
          key.title = NA,
          keysize = 1)
dev.off()
lung_ordered_MBAmethyl = lung_res_MBAmethyl$aic$theta[,lung_MBAmethyl$colInd]

#######heatmap without clustering##########
encode_cnvJoint = heatmap.2(encode_res$aic$theta,
                     Colv = FALSE, Rowv = FALSE,
                     dendrogram = "none",
                     main = 'ENCODE_cnvJoint',
                     tracecol = NA,
                     labRow=FALSE,
                     denscol="black",
                     key.title = NA,
                     keysize = 1
)
encode_MBAmethyl = heatmap.2(encode_res_MBAmethyl$aic$theta,
                     Colv = FALSE, Rowv = FALSE,
                     dendrogram = "none",
                     main = 'ENCODE_MBAmethyl',
                     tracecol=NA,
                     labRow=FALSE,
                     denscol="black",
                     key.title = NA,
                     keysize = 1)
poly_cnvJoint = heatmap.2(poly_res$aic$theta,
                     Colv = FALSE, Rowv = FALSE,
                     dendrogram = "none",
                     main = 'DOP-PCR_cnvJoint',
                     tracecol = NA,
                     labRow=FALSE,
                     denscol = "black",
                     key.title = NA,
                     keysize = 1)
poly_MBAmethyl = heatmap.2(poly_res_MBAmethyl$aic$theta,
                     Colv = FALSE, Rowv = FALSE,
                     dendrogram = "none",
                     main = "DOP-PCR_MBAmethyl",
                     tracecol = NA,
                     labRow=FALSE,
                     denscol = "black",
                     key.title = NA,
                     keysize = 1)

lung_cnvJoint = heatmap.2(lung_res$aic$theta,
                     Colv = FALSE, Rowv = FALSE,
                     dendrogram = "none",
                     main = "MALBAC_cnvJoint",
                     tracecol = NA,
                     labRow=FALSE,
                     denscol="black",
                     key.title = NA,
                     keysize = 1)

lung_MBAmethyl = heatmap.2(lung_res_MBAmethyl$aic$theta,
                     Colv = FALSE, Rowv = FALSE,
                     dendrogram = "none",
                     main = "MALBAC_MBAmethyl",
                     labRow=FALSE,
                     tracecol = NA,
                     xlab = NA,
                     denscol="black",
                     key.title = NA,
                     keysize = 1)


####visualize distance####
get_distance_matrix = function(cluster, mat){
  diffmat_cnvJoint = matrix(0, length(cluster), length(cluster))
  for (i in 1:(length(cluster)-1)){
    for (j in (i+1):length(cluster)){
      vec1 = mat[,cluster[i]]
      vec2 = mat[,cluster[j]]
      diffmat_cnvJoint[i,j] = diffmat_cnvJoint[j,i] = sum((vec1-vec2)^2)/length(vec1)
    }
  }
  return(list(
    diff = sum(diffmat_cnvJoint)/(length(cluster)*(length(cluster)-1)),
    diffmat = diffmat_cnvJoint))
}

get_cross_distance_sum = function(cluster1, cluster2, mat){
  diff=0
  for (i in 1:length(cluster1)){
    for (j in 1:length(cluster2))
      vec1 = mat[,cluster1[i]]
      vec2 = mat[,cluster2[j]]
      diff = diff + sum((vec1-vec2)^2)/length(vec1)
  }
  return=diff
}

visualize_distance = function(mat1, mat2, titlecnvJoint, titleMBAmethyl, clust){
  t=get_distance_matrix(clust, mat1)
  t_MBAmethyl=get_distance_matrix(clust, mat2)
  nd = melt(t$diffmat)
  od = melt(t_MBAmethyl$diffmat)
  g1 = ggplot(nd, aes(Var1, Var2)) +
    geom_tile(aes(fill = value),colour = "white") +
    scale_fill_gradient(low = "white", high = "indianred",
                        guide= guide_colorbar(title.hjust = 0.2, title=''),
                        limits = c(0,max(max(nd$value), max(od$value))))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(title=titlecnvJoint)
  g2 = ggplot(od, aes(Var1, Var2)) +
    geom_tile(aes(fill = value),colour = "white") +
    scale_fill_gradient(low = "white", high = "indianred",
                        guide= guide_colorbar(title.hjust = 0.2, title=''),
                        limits = c(0,max(max(nd$value), max(od$value))))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(title=titleMBAmethyl)
  gg = grid.arrange(g2, g1, nrow=1)
  return(graph = gg)
}
pdf('analysis/writeup/paper/encode_distmat_heat.pdf', width=11, height=5)
visualize_distance(encode_ordered_cnvJoint, encode_ordered_MBAmethyl, 'ENCODE_cnvJoint', 'ENCODE_MBAmethyl', 1:32)
dev.off()
pdf('analysis/writeup/paper/poly_distmat_heat.pdf', width=11, height=5)
visualize_distance(poly_ordered_cnvJoint, poly_ordered_MBAmethyl, 'DOPPCR_cnvJoint', 'DOPPCR_MBAmethyl', 1:100)
dev.off()
pdf('analysis/writeup/paper/lung_distmat_heat.pdf', width=11, height=5)
visualize_distance(lung_ordered_cnvJoint, lung_ordered_MBAmethyl, 'MALBAC_cnvJoint', 'MALBAC_MBAmethyl', 1:29)
dev.off()
####get clusters####

get_splits = function(res){
  g = heatmap.2(res$aic$theta, Colv=TRUE, Rowv=FALSE, dendrogram = "column")
  print(g$colDendrogram[[1]])
  print(g$colDendrogram[[2]])
  return(g$colInd)
}
encode_colInd = get_splits(encode_res)
encode_cl1 = encode_colInd[1:21]
encode_cl2 = encode_colInd[22:32]
encode_MBAmethyl_colInd = get_splits(encode_res_MBAmethyl)
encode_MBAmethyl_cl1 = encode_MBAmethyl_colInd[1:21]
encode_MBAmethyl_cl2 = encode_MBAmethyl_colInd[22:32]

poly_colInd = get_splits(poly_res)
poly_cl1 = poly_colInd[1:63]
poly_cl2 = poly_colInd[64:100]
poly_MBAmethyl_colInd = get_splits(poly_res_MBAmethyl)
poly_MBAmethyl_cl1 = poly_MBAmethyl_colInd[1:64]
poly_MBAmethyl_cl2 = poly_MBAmethyl_colInd[65:100]

lung_colInd = get_splits(lung_res)
lung_cl1 = lung_colInd[1:8]
lung_cl2 = lung_colInd[9:29]
lung_MBAmethyl_colInd = get_splits(lung_res_MBAmethyl)
lung_MBAmethyl_cl1 = lung_colInd[1:11]
lung_MBAmethyl_cl2 = lung_colInd[12:29]

encode_intradist =
  poly_intradist =
  lung_intradist = data.frame(cnvJoint = c(0,0), MBAmethyl = c(0,0))

encode_intradist$cnvJoint[1] = get_distance_matrix(encode_cl1, encode_res$aic$theta)$diff
encode_intradist$cnvJoint[2] = get_distance_matrix(encode_cl2, encode_res$aic$theta)$diff
encode_intradist$MBAmethyl[1] = get_distance_matrix(encode_MBAmethyl_cl1, encode_res_MBAmethyl$aic$theta)$diff
encode_intradist$MBAmethyl[2] = get_distance_matrix(encode_MBAmethyl_cl2, encode_res_MBAmethyl$aic$theta)$diff
encode_intradist$clustindex = c('ENCODEclust1', 'ENCODEclust2')
encode_intradist = melt(encode_intradist)
poly_intradist$cnvJoint[1] = get_distance_matrix(poly_cl1, poly_res$aic$theta)$diff
poly_intradist$cnvJoint[2] = get_distance_matrix(poly_cl2, poly_res$aic$theta)$diff
poly_intradist$MBAmethyl[1] = get_distance_matrix(poly_MBAmethyl_cl1, poly_res_MBAmethyl$aic$theta)$diff
poly_intradist$MBAmethyl[2] = get_distance_matrix(poly_MBAmethyl_cl2, poly_res_MBAmethyl$aic$theta)$diff
poly_intradist$clustindex = c('DOPPCRclust1', 'DOPPCRclust2')
poly_intradist = melt(poly_intradist)
lung_intradist$cnvJoint[1] = get_distance_matrix(lung_cl1, lung_res$aic$theta)$diff
lung_intradist$cnvJoint[2] = get_distance_matrix(lung_cl2, lung_res$aic$theta)$diff
lung_intradist$MBAmethyl[1] = get_distance_matrix(lung_MBAmethyl_cl1, lung_res_MBAmethyl$aic$theta)$diff
lung_intradist$MBAmethyl[2] = get_distance_matrix(lung_MBAmethyl_cl2, lung_res_MBAmethyl$aic$theta)$diff
lung_intradist$clustindex = c('MALBACclust1', 'MALBACclust2')
lung_intradist = melt(lung_intradist)

intradist = rbind(encode_intradist, poly_intradist, lung_intradist)
pdf('analysis/writeup/paper/intraclusters_bar.pdf', width=12, height=6)
ggplot(intradist,aes(x=clustindex,y=value,fill=factor(variable)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="",
                      labels=c("cnvJoint", "MBAmethyl"))+
  xlab("clusters")+ylab("average euclidean distance")+
  ggtitle("(a)")
dev.off()

crossdist = data.frame(clustindex=character(6), variable=character(6), value=numeric(6))
crossdist$clustindex = rep(c('ENCODE', 'DOPPCR', 'MALBAC'), each = 2)
crossdist$variable = rep(c('cnvJoint','MBAmethyl'), 3)
crossdist$value = c(
  get_cross_distance_sum(encode_cl1, encode_cl2, encode_res$aic$theta),
  get_cross_distance_sum(encode_MBAmethyl_cl1, encode_MBAmethyl_cl2, encode_res_MBAmethyl$aic$theta),
  get_cross_distance_sum(poly_cl1, poly_cl2, poly_res$aic$theta),
  get_cross_distance_sum(poly_MBAmethyl_cl1, poly_MBAmethyl_cl2, poly_res_MBAmethyl$aic$theta),
  get_cross_distance_sum(lung_cl1, lung_cl2, lung_res$aic$theta),
  get_cross_distance_sum(lung_MBAmethyl_cl1, lung_MBAmethyl_cl2, lung_res_MBAmethyl$aic$theta)
)
pdf('analysis/writeup/paper/acrossclusters_bar.pdf', width=12, height=6)
ggplot(crossdist,aes(x=clustindex,y=value,fill=factor(variable)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="",
                      labels=c("cnvJoint", "MBAmethyl"))+
  xlab("clusters")+ylab("average euclidean distance")+
  ggtitle("(b)")
dev.off()

####plot xi####

encode_xi = encode_res$aic$xi[1,]
poly_xi = poly_res$aic$xi[1,]
lung_xi = lung_res$aic$xi[1,]

xi = data.frame(x=c(encode_xi,poly_xi,lung_xi),
                lab = c(rep('ENCODE', 32), rep('DOP-PCR', 100), rep('MALBAC', 29)))
xi$lab = factor(xi$lab, levels=c('ENCODE','DOP-PCR','MALBAC'),ordered = TRUE)
ggplot(xi, aes(x=lab, y=x)) + geom_boxplot()+
  ggtitle(expression("distribution of"~xi))


####observe spikes####
for (i in 1:32){
  plot(encode_res$aic$theta[,i], ylim = c(-5,5), type = 'l', main = i)
  #lines(encode_res_MBAmethyl$aic$theta[,i], col = 'red')
  lines(encode_Y[, i], col = 'red')
  Sys.sleep(0.5)
}

