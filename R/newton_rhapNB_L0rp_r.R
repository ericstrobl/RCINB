
newton_rhapNB_L0rp_r <- function(X,Y,Z,Cb,dp,warmStart){
  
  LL = list()
  LL$X = X; LL$Y = Y; 
  LL$Z = Z; 
  d = ncol(X)-dp
  if (d>0){LL$max = apply(X[,seq_len(d),drop=FALSE],2,max)}
  LL$W = rep(1,d+dp)
  LL$Cb = Cb
  
  n = nrow(X)
  LL$L = 1*log(n)/n; LL$d = d; LL$dp = dp
  
  sol_alpha = multiroot(derivsp_alpha_L2rp_r, warmStart, parms = LL, positive=TRUE, rtol = 1e-5, atol = 1e-6, ctol = 1e-6)
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
    sol_alpha = multiroot(derivsp_alpha_L2rp_r, alpha, parms = LL, positive=TRUE, rtol = 1e-5, atol = 1e-6, ctol = 1e-6)
    alphan = sparsify(as.matrix(-sol_alpha$root));
    if (LL$d>0){alphan[1:LL$d] = alphan[1:LL$d] + 0/LL$max}
    alphan[(LL$d+1):(LL$d+LL$dp)] = alphan[(LL$d+1):(LL$d+LL$dp)] + 5
    
    delta = sum(abs(alphan - alpha))
    alpha = alphan
    LL$W = c(abs(alpha[seq_len(d)]),rep(1,dp)) # L2 regularization for ppl
  }
  
  LL$prod = exp(Z%*%alpha)*Cb
  sol_r = multiroot(derivsp_r_var_r, 0.1, parms = LL,positive=TRUE, rtol = 1e-5, atol = 1e-6, ctol = 1e-6)
  r = sol_r$root
  
  return( list(alpha = alpha, Z = Z, r = r) )
  
}

sparsify <- function(alpha){
  ix = which(abs(alpha)<1E-2)
  alpha[ix] = 0
  return(alpha)
}

derivsp_alpha_L2rp_r <- function(alpha, LL){
  alpha = -alpha
  if (LL$d>0){alpha[1:LL$d] = alpha[1:LL$d] + 0/LL$max}
  alpha[(LL$d+1):(LL$d+LL$dp)] = alpha[(LL$d+1):(LL$d+LL$dp)] + 5
  prod = exp(LL$Z%*%alpha)*LL$Cb ##
  MM = colMeans(c(LL$Y)*LL$X)*LL$W - colMeans(c(prod)*LL$Z)*LL$W
  alpha[seq_len(LL$d)] = MM[seq_len(LL$d)] - LL$L*alpha[seq_len(LL$d)]
  
  alphap = alpha[(LL$d+1):(LL$d+LL$dp)]; alphap = (alphap-mean(alphap))*(1-(1/LL$dp))
  alpha[(LL$d+1):(LL$d+LL$dp)] = MM[(LL$d+1):(LL$d+LL$dp)] - LL$L*alphap
  
  # prod = exp(LL$Z%*%alpha)
  # alpha = colMeans(c(LL$Y)*LL$X)*LL$W-colMeans(c(prod)*LL$Z)*LL$W - LL$L*alpha
  
  return( alpha )
}

derivsp_r_var_r <- function(r, LL){
  # r = r + 1E-3
  # r = 1/(1+exp(-r))*10+1E-10 # squash, strictly positive (0,100)
  # print(ncol((digamma(Xp %*% r+LL$Y+1E-10)-digamma(Xp %*% r+1E-10)+log(Xp %*% r+1E-10))*Xp))
  
  r1 = mean( digamma(r+LL$Y+1E-10)-digamma(r+1E-10)+log(r+1E-10) )- mean( log(r+LL$prod+1E-10) )
  # print(r1)
  return( r1 )
}