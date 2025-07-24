
data {
  
  int<lower=2> nCause;
  int<lower=1> nAlgo;
  int<lower=0> aj[nAlgo,nCause];
  
  simplex[nCause] p_uncalib;
  simplex[nCause] Mmat_bysimplex[nAlgo,nCause];
  real<lower=0> pss;
  
}

parameters {
  
  simplex[nCause] p_calib;
  
}

transformed parameters {
  
  matrix[nCause,nCause] Mmat[nAlgo];
  simplex[nCause] q[nAlgo];
  vector[nAlgo] loglik;
  
  for(k in 1:nAlgo){
    
    for(i in 1:nCause){
      
      Mmat[k][i,] = Mmat_bysimplex[k,i]';
      
    }
  
    q[k] = (Mmat[k]')*p_calib;
  
    loglik[k] = multinomial_lpmf(aj[k,] | q[k]);
  
  }
  
}

model {
  
  target += dirichlet_lpdf(p_calib | 1 + nCause*pss*p_uncalib);
  
  target += loglik;
  
}
