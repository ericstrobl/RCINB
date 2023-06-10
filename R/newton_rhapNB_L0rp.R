newton_rhapNB_L0rp <- function(X,Y,Z,Cb,dp,warmStart){
  
  LL = list()
  LL$X = X; LL$Y = Y; 
  LL$Z = Z; 
  d = ncol(X)-dp
  if (d>0){LL$max = apply(X[,seq_len(d),drop=FALSE],2,max)}
  LL$W = rep(1,d+dp)
  LL$Cb = Cb
  
  n = nrow(X)
  LL$L = 1*log(n)/n; LL$d = d; LL$dp = dp
  
  sol_alpha = multiroot(derivsp_alpha_L2rp, warmStart, parms = LL, positive=TRUE, rtol = 1e-5, atol = 1e-6, ctol = 1e-6)
  # print(sol_alpha)
  alpha = -as.matrix(sol_alpha$root); 
  if (LL$d>0){alpha[1:LL$d] = alpha[1:LL$d] + 0/LL$max}
  alpha[(LL$d+1):(LL$d+LL$dp)] = alpha[(LL$d+1):(LL$d+LL$dp)] + 5
  LL$W = c(abs(alpha[seq_len(d)]),rep(1,dp))
  
  # if (LL$d>0){print(alpha[1:LL$d])}
  
  delta = 1
  imax = 0
  while ((delta > 1E-8)&(imax<1)){
    imax = imax+1
    sol_alpha = multiroot(derivsp_alpha_L2rp, alpha, parms = LL, positive=TRUE, rtol = 1e-5, atol = 1e-6, ctol = 1e-6)
    alphan = sparsify(as.matrix(-sol_alpha$root));
    if (LL$d>0){alphan[1:LL$d] = alphan[1:LL$d] + 0/LL$max}
    alphan[(LL$d+1):(LL$d+LL$dp)] = alphan[(LL$d+1):(LL$d+LL$dp)] + 5
    
    delta = sum(abs(alphan - alpha))
    alpha = alphan
    LL$W = c(abs(alpha[seq_len(d)]),rep(1,dp)) # L2 regularization for ppl
  }
  
  r = exp(alpha[(d+1):(d+dp)])
  
  return( list(alpha = alpha, Z = Z, r =r) )
  
}

derivsp_alpha_L2rp <- function(alpha, LL){
  alpha = -alpha
  if (LL$d>0){alpha[1:LL$d] = alpha[1:LL$d] + 0/LL$max}
  alpha[(LL$d+1):(LL$d+LL$dp)] = alpha[(LL$d+1):(LL$d+LL$dp)] + 5
  prod = exp(LL$Z%*%alpha)*LL$Cb ##
  MM = colMeans(c(LL$Y)*LL$X)*LL$W - colMeans(c(prod)*LL$Z)*LL$W
  alpha[seq_len(LL$d)] = MM[seq_len(LL$d)] - LL$L*alpha[seq_len(LL$d)]
  
  alphap = alpha[(LL$d+1):(LL$d+LL$dp)]; alphap = (alphap-mean(alphap))*(1-(1/LL$dp))
  alpha[(LL$d+1):(LL$d+LL$dp)] = MM[(LL$d+1):(LL$d+LL$dp)] - LL$L*alphap
  
  return( alpha )
}
 
