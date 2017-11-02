#post processing
postprocess = function(Y, theta, phi, xi){
  out = theta
  resid = Y - theta*phi - phi %*% t(xi)
  for (j in 1:ncol(Y)){
    mod = lm(Y[,j]~1)
    cooksd = cooks.distance(mod)
    ind = which(cooksd>5*mean(cooksd, na.rm=TRUE))
    out[ind, j] = out[ind,j]+resid[ind,j]/phi[ind]
  }
  return(out)
}
newtheta = postprocess(Y, theta, phi, xi)
plot(Y[,1], ylim = c(-10,5))
lines(newtheta[,1], col = 'red', lwd=2)
lines(theta[,1], col='blue')
