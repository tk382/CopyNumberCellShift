######polygenomic breast tumor#########
orig = read.table('/Users/tae/Dropbox/TaeProject/CopyNumber/polygenomic_breast_tumor_copynumber.txt', header = TRUE, stringsAsFactors = FALSE)
orig = orig[orig$CHR=='chr22',-(1:3)]
Y_center = orig - 2
diffvec = apply(Y_center, 2, function(x) sum(x^2)/nrow(orig))
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
Y = log(orig2 / nor)

######circulating lung tumor cells#########
orig = read.table('../circulating_lung_tumor_copynumber.txt', header = TRUE, stringsAsFactors = FALSE)
orig = orig[orig$CHR=='chr22',-(1:3)]
Y_center = orig - 2
diffvec = apply(Y_center, 2, function(x) sum(x^2)/nrow(Y))
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
Y = log(orig2 / nor)

