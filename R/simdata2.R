simdata2 = function(p,n){
  k = round(n/2)
  sigma = runif(p, .1, 1)
  phi = rnorm(p)
  theta = matrix(0, p, n)
  cp = sample(1:(p-1), k)
  cp = c(0, cp, p)
  cp = sort(cp)
  for (c in 2:length(cp)){
    change = rbinom(1, 1, 0.2)
    theta[(cp[c-1]+1):cp[c], ] = change*sample(seq(-3, 3), 1)
  }
  #theta[cnv_probes, cancer_cells] = theta[cnv_probes, non_cancer_cells[1]] + rep(0.5, length(cnv_probes))
  Y = matrix(0, p, n)
  phi = rnorm(p,0,1)
  xi = rnorm(n,0,1)
  for (i in 1:p){
    for (j in 1:n){
      Y[i,j] = (theta[i,j] + xi[j]) * phi[i] + rnorm(1)
    }
  }
  return(list(Y=Y, theta=theta, phi=phi, xi=xi, cp = cp))
}
