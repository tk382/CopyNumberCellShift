#post processing
postprocess = function(Y, theta, phi, xi){
  out = theta
  resid = Y - theta*phi - phi %*% t(xi)
  for (j in 1:ncol(Y)){
    mod = lm(resid[,j]~1)
    cooksd = cooks.distance(mod)
    ind = which(cooksd>10*mean(cooksd, na.rm=TRUE))
    out[ind, j] = out[ind,j]+resid[ind,j]/phi[ind]
  }
  return(out)
}
