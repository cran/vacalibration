
data {
  
  int<lower=2> nCause;
  int<lower=1> nAlgo;
  int<lower=0> aj[nAlgo,nCause];
  int<lower=2> nTocalib[nAlgo];
  int<lower=1,upper=nCause> idtocalib[sum(nTocalib)];
  int<lower=0,upper=(nAlgo*nCause)> nCumtocalib[nAlgo+1];
  
  simplex[nCause] p_uncalib;
  vector<lower=0>[nCause] Mmatprior_asDirich[nAlgo,nCause];
  real<lower=0> pss;
  
}

parameters {
  
  simplex[nCause] Mmat_bysimplex[nAlgo,nCause];
  simplex[nCause] p_calib;
  
}

transformed parameters {
  
  matrix[nCause,nCause] Mmat[nAlgo];
  simplex[nCause] q[nAlgo];
  vector[nAlgo] loglik;
  
  for(k in 1:nAlgo){
    
    Mmat[k] = diag_matrix(rep_vector(1.0, nCause));
    
    for(i in 1:nTocalib[k]){
      
      Mmat[k][idtocalib[nCumtocalib[k]+i],idtocalib[(nCumtocalib[k]+1):nCumtocalib[k+1]]] = ((Mmat_bysimplex[k,idtocalib[nCumtocalib[k]+i]])[idtocalib[(nCumtocalib[k]+1):nCumtocalib[k+1]]]/sum((Mmat_bysimplex[k, idtocalib[nCumtocalib[k]+i]])[idtocalib[(nCumtocalib[k]+1):nCumtocalib[k+1]]]))';
      
    }
  
    q[k] = (Mmat[k]')*p_calib;
  
    loglik[k] = multinomial_lpmf(aj[k,] | q[k]);
  
  }
  
}

model {
  
  target += dirichlet_lpdf(p_calib | 1 + nCause*pss*p_uncalib);
  
  for(k in 1:nAlgo){
  
    for(i in 1:nCause){
      
      target += dirichlet_lpdf(Mmat_bysimplex[k,i] | Mmatprior_asDirich[k,i]);
      
    }
    
  }
  
  target += loglik;
  
}
