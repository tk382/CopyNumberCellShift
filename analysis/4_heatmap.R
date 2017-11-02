#####load packages and functions#####
library(lme4)
library(dplyr)
library(lmvar)
library(ggplot2)
library(reshape2)
library(gplots)
library(gridExtra)
#library(CopyNumberCellShift)
setwd('/Users/tae/Dropbox/TaeProject/CopyNumber/CopyNumberCellShift')
Rcpp::sourceCpp('src/cnvJoint.cpp')
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)


#####set Y ######
j =  22 #chromosome ID
load("../summary.100kb.normalized.RData")
X = data.100kb[data.100kb[, 1]%in%paste0("chr", j), -c(1:3)]
Y = as.matrix(apply(X, 2, as.numeric))

wts = as.numeric(defaultWeights_c(nrow(Y)))
steps = min(nrow(Y),ncol(Y))-1
res = cnvJoint(Y = as.matrix(Y), wts = wts, steps = steps, maxloop = 30, verbose=TRUE)
res_mba = cnv_c_old(as.matrix(Y), wts, steps, 30)




heatmap.2(Y,
          Colv=FALSE,
          Rowv=FALSE,
          dendrogram = "none",
          main = 'Raw Data',
          tracecol = NA,
          labRow=FALSE,
          col=my_palette,
          breaks=seq(-8,8,length.out=300),
          denscol="black",
          key.title = NA,
          keysize = 2
)
colnames(res$bic$theta) = colnames(Y)
heatmap.2(res$bic$theta,
          Colv=FALSE,
          Rowv=FALSE,
          dendrogram = "none",
          main = 'New Result',
          col = my_palette,
          breaks=seq(-8,8,length.out=300),
          tracecol = NA,
          labRow=FALSE,
          denscol="black",
          key.title = NA,
          keysize = 2
          )

colnames(res_mba$bic$theta) = colnames(Y)
heatmap.2(res_mba$bic$theta,
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram = "none",
          main= 'Old Result',
          tracecol = NA,
          col = my_palette,
          breaks=seq(-8,8,length.out=300),
          labRow = FALSE,
          denscol="black",
          key.title = NA,
          keysize=2)

gink = read.table('../Seb_processed_ginkgo/SegFixed.txt', header = TRUE)
gink = gink[gink$CHR=='chr22', ]
colnames(gink)[-(1:3)] = colnames(Y)
gink2 = as.matrix(gink[,-(1:3)])
gink2 = gink2 * 2
gink2 = log(gink2/2+0.01)
gink2 = gink2 + 0.01
heatmap.2(as.matrix(gink2),
          Colv = FALSE,
          Rowv = FALSE,
          col= my_palette,
          breaks=seq(-8,8,length.out=300),
          dendrogram = "none",
          main= 'Ginkgo Result',
          tracecol = NA,
          labRow = FALSE,
          denscol="black",
          key.title = NA,
          keysize=2)
