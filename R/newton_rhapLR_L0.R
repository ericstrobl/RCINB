
newton_rhapLR_L0 <- function(X,Y,Z,dp){
  
  LL = list()
  LL$X = X; LL$Y = Y; 
  LL$Z = Z; 
  d = ncol(X)-dp
  LL$W = rep(1,d+dp)
  
  n = nrow(X)
  LL$L = log(n)/n; LL$d = d; LL$dp = dp
  
  sol_alpha = multiroot(derivsp_alpha_L2_LR, rep(0,d+dp), parms = LL)
  alpha = as.matrix(sol_alpha$root)
  # LL$W = c(abs(alpha[seq_len(d)]),rep(1,dp))
  
  # delta = 1
  # imax = 0
  # while ( (delta > 1E-8) & (imax<1) ) {
  #   imax = imax + 1
  # 
  #   sol_alpha = multiroot(derivsp_alpha_L2_LR, alpha, parms = LL)
  #   alphan = sparsify(as.matrix(sol_alpha$root))
  # 
  #   delta = sum(abs(alphan - alpha))
  #   alpha = alphan
  #   LL$W = c(abs(alpha[seq_len(d)]),rep(1,dp)) # L2 regularization for ppl
  # }

  return( list(alpha = alpha, Z = Z) )
  
}

derivsp_alpha_L2_LR <- function(alpha, LL){
  prod = c(exp(-LL$Z%*%alpha)) ##
  MM = colMeans(c(LL$Y)*LL$X)*LL$W - colMeans((1/(1+prod))*LL$Z)*LL$W
  alpha[seq_len(LL$d)] = MM[seq_len(LL$d)] - LL$L*alpha[seq_len(LL$d)]
  
  alphap = alpha[(LL$d+1):(LL$d+LL$dp)]; alphap = (alphap-mean(alphap))*(1-(1/LL$dp))
  alpha[(LL$d+1):(LL$d+LL$dp)] = MM[(LL$d+1):(LL$d+LL$dp)] - LL$L*alphap
  
  return( alpha )
}
