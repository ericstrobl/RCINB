test_NB_L0_mom00rp2sC3_fast_r <- function(X,Y,Z,Cb,dp,warmStart){
  
  # test_NB_L0_mom00rp2sC2 <- function(X,Y,C,Ip,shapes=NULL,rates=NULL,DAG=NULL){
  # Ip is vector of person indices, e.g., c(1,1,1,2,2,2,2,2,2,3,3,3)
  
  require(Rfast)
  require(rootSolve)
  require(MASS)
  
  if (length(X)>0) X = as.matrix(X)
  
  # remove outliers
  ix = which(Y>200)
  if (length(ix)>0){
    X = X[-ix,,drop=FALSE]
    Y = Y[-ix]
    d = ncol(Z)-dp
    if (d == 0){
      Cb = Cb[-ix]
      Z = Z[-ix,,drop=FALSE]
    }
  }
  
  # learn model
  # nb = trainNB_L0rpC3(X,Y,C,Xp,shapes,rates,DAG)
  nb = newton_rhapNB_L0rp_r(X,Y,Z,Cb,dp,warmStart)
  alpha = nb$alpha
  r = nb$r
  
  # print(r)
  
  d = ncol(Z); n = length(Y); m = nrow(Z)
  
  mu = pmax(c(exp(Z%*%alpha)*Cb),1E-10)
  scores1 = cbind( c(Y)*X , digamma(r+Y+1E-10) )
  scores2 = cbind( -c(mu)*Z , -log(r+mu) )
  
  
  a = 1; b = 1
  stats1 = cbind(exp(-a*Y),sin(b*Y))
  stats2 = cbind(((exp(a)*r)/(exp(a)*(mu + r) - mu))^r, Re(-0.5i*((r/(-exp(1i*b)*mu + mu + r))^r - ((exp(1i*b)*r)/(-mu + exp(1i*b)*(mu + r)))^r)))
  # page 47 in Rayner, Thas, Best. Smooth Goodness Tests of Fit 2nd Edition
  e = ncol(stats1)
  A = matrix(0,e+d+1,e+d+1)
  
  A[(e+1):(e+d),(e+1):(e+d)] = (-t(mu*Z)%*%Z)/m
  
  A[1,(e+1):(e+d)] = -colMeans( -((exp(a) - 1)*r*Z*mu*(r/((exp(a) - 1)*mu*exp(-a) + r))^r)/(exp(a)*(mu + r) - mu) )
  A[1,e+d+1] = -mean( (exp(-a)*((exp(a)*r)/(exp(a)*(mu + r) - mu))^(r + 1)*((exp(a)*(mu + r) - mu)*log((exp(a)*r)/(exp(a)*(mu + r) - mu)) + (exp(a) - 1)*mu))/r )
  A[2,(e+1):(e+d)] = -colMeans(Re( -0.5i*((exp(1i*b)*r^2*(-Z*mu + Z*mu*exp(1i*b))*((exp(1i*b)*r)/(-mu + exp(1i*b)*(mu + r)))^(r - 1))/(-mu + exp(1i*b)*(mu + r))^2 - (r^2*(Z*mu - Z*mu*exp(1i*b))*(r/(-mu*exp(1i*b) + mu + r))^(r - 1))/(-mu*exp(1i*b) + mu + r)^2)   ))
  A[2,e+d+1] = -mean(Re( -0.5i*((r/(-exp(1i*b)*mu + mu + r))^r*((-exp(1i*b)*mu + mu + r)*(1/(-exp(1i*b)*mu + mu + r) - r/(-exp(1i*b)*mu + mu + r)^2) + log(r/(-exp(1i*b)*mu + mu + r))) - ((exp(1i*b)*r)/(-mu + exp(1i*b)*(mu + r)))^r*(exp(-1i*b)*(-mu + exp(1i*b)*(mu + r))*(exp(1i*b)/(-mu + exp(1i*b)*(mu + r)) - (exp(2i*b)*r)/(-mu + exp(1i*b)*(mu + r))^2) + log((exp(1i*b)*r)/(-mu + exp(1i*b)*(mu + r))))) ))
  
  A[(e+1):(e+d),(e+1):(e+d)] = (-t(mu*Z)%*%Z)/m
  A[(e+d+1),(e+1):(e+d)] = colMeans(-mu/(r+mu)*Z)
  A[e+d+1,e+d+1] = mean(trigamma(r+Y+1e-10))-trigamma(r+1e-10)+1/r-mean(1/(r+mu))
  A = -A
  
  if ((d-dp)>0){
    B1 = cov(cbind(stats1,scores1))
    B2 = cov(cbind(-stats2,scores2))
    B = B1/n+B2/m
  } else{
    B = cov(cbind(stats1-stats2,scores1+scores2))/n
  }
  
  g = 1:e; b = (e+1):(e+d+1)
  Sigma = B[g,g,drop=FALSE]-
    A[g,b,drop=FALSE] %*% ginv(A[b,b,drop=FALSE]) %*% B[b,g,drop=FALSE] -
    B[g,b,drop=FALSE] %*% t(ginv(A[b,b,drop=FALSE])) %*% t(A[g,b,drop=FALSE]) +
    A[g,b,drop=FALSE] %*% ginv(A[b,b,drop=FALSE]) %*% B[b,b,drop=FALSE] %*% t(ginv(A[b,b,drop=FALSE])) %*% t(A[g,b,drop=FALSE])
  
  h = colMeans(stats1)-colMeans(stats2)
  S = t(h) %*% ginv(Sigma) %*% h
  
  p = 1-pchisq(S,length(h))
  
  # return( list(alpha = alpha, r = c(exp(Xp%*%alphap)), p = p, S = S) )
  return( list(alpha = alpha, r = nb$r, p = p, S = S) )
  # return( list(alpha = alpha, r = nb$r, p = 0, S = S) )
}
