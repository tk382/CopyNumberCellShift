simulationdata = function(p, n, probe_var, cell_var){
  k = round(n/2)
  delta = runif(n, .1, 1)
  sigma = runif(p, .1, 1)
  phi = rnorm(p)
  theta = matrix(0, p, n)
  cp = sample(1:(p-1), k)
  cp = c(0, cp, p)
  cp = sort(cp)
  for (c in 2:length(cp)){
      theta[(cp[c-1]+1):cp[c], ] = sample(seq(-3, 3), 1)
  }
  cancer_cells = sample(1:n, round(n/3))
  non_cancer_cells = (1:n)[-cancer_cells]
  cnv_probes = round(p/3):round(p/2)
  #theta[cnv_probes, cancer_cells] = theta[cnv_probes, non_cancer_cells[1]] + rep(0.5, length(cnv_probes))
  Y = matrix(0, p, n)
  error = matrix(0, p, n)
  if(cell_var==TRUE & probe_var==TRUE){
    for (i in 1:p){
      for (j in 1:n){
        error[i,j] = rnorm(1, 0, delta[j] + sigma[i])
      }
    }
  }
  if(cell_var==TRUE & probe_var==FALSE){
    for (j in 1:n){
        error[,j] = rnorm(p, 0,delta[j])
    }
  }
  if(cell_var==FALSE && probe_var == TRUE){
    for (i in 1:p){
      error[i, ] = rnorm(n, 0, sigma[i])
    }
  }
  if(cell_var==FALSE & probe_var == FALSE){
    error = matrix(rnorm(n*p, 0, 1), nrow=p)
  }
  Y = phi*theta+error
  return(list(Y=Y, cp = cp, theta = theta, phi = phi, delta = delta, sigma = sigma, cancer_cells = cancer_cells))
}

