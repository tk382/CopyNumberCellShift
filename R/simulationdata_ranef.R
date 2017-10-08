simulationdata_ranef = function(p, n){
  k = round(n/2)
  delta = rnorm(n, 0, 1)
  phi = rnorm(p)
  theta = matrix(0, p, n)
  cp = sample(1:(p-1), k)
  cp = c(0, cp, p)
  cp = sort(cp)
  for (c in 2:length(cp)){
    theta[(cp[c-1]+1):cp[c], ] = sample(seq(-3, 3), 1)
  }
  Y = matrix(0, p, n)
  error = matrix(rnorm(p*n), p, n)
  for (j in 1:n){
    error[,j] = error[,j] + delta[j]
  }
  Y = phi*theta+error
  return(list(Y=Y, 
              cp = cp, 
              theta = theta, 
              phi = phi, 
              delta = delta, 
              sigma = sigma)
         )
}
