RCI_NB_G <- function(data,Yi,C,Xp,alpha=0.02,Gtr){
  
  dp = ncol(Xp)
  X = data[,-Yi,drop=FALSE]
  d = ncol(X)
  
  Y = data[,Yi,drop=FALSE]
  
  iS = c() # predictors
  iV = seq_len(d) # targets
  
  n = nrow(X)
  shapes = rep(0,d); rates = matrix(1,dp,d) # for gamma distribution for simulation
  G = matrix(0,d,d) # start with null graph
  warmStart = matrix(0,d,d+dp)
  
  Z = X; Zp = Xp; Cb = C
  iP = (d+1):(d+dp);iO=1
  while ( length(iV)>0 ){ # until no targets remain
    iR = c() # indices to remove from iV
    s_max = Inf; cnt = 0
    if (iO==1){
      iit = which(colSums(Gtr[,-Yi])==0) # first iter, get roots; only works if Y is in the last index
      iO = iO+1
    } else{
      iit = which(colSums(Gtr[,-Yi])!=0) # other iter, 
    }
    for (v in iit){
      
      iSs = which(Gtr[,v]!=0) # parents of v
      if ( sum(colSums(Gtr[,iSs,drop=FALSE])!=0)==0 ){ # parents of v without parents themselves
        
        oo = test_NB_L0_mom00rp2sC3_fast_r( cbind(X[,iSs,drop=FALSE]/C, Xp), X[,v,drop=FALSE],
                                          cbind(Z[,iSs,drop=FALSE], Zp), Cb, dp, warmStart[v,c(iSs,iP)]) # regress predictors on one target
        
        iR = c(iR,v) # remove index v
        
        shapes[v] = oo$r
        rates[,v] = oo$r/exp(oo$alpha[(length(iSs)+1):(length(iSs)+dp)])
        
        G[iSs,v] = oo$alpha[seq_len(length(iSs))]
      }
      
      
    }
    
    # print(oo$alpha)
    iV = setdiff(iV, iR) # remove iR from iV
    
    # sample Z
    iS = c(iS,iR) # include iR into predictors
    
    oo = sample_gamma_var(n, (length(iS)-length(iR)+1):length(iS), Z[,iS,drop=FALSE], Zp, Cb, 
                          shapes=shapes[iS],rates=rates[,iS,drop=FALSE],DAG=G[iS,iS,drop=FALSE])
    if (nrow(Z)==n){
      Z = matrix(0,nrow(oo$Z),ncol(Z))
      E = matrix(0,nrow(oo$Z),ncol(Z))
    }
    Z[,iR] = oo$Z; Zp = oo$Zp; Cb = oo$Cb; E[,iR] = oo$E
    
    Gtr[,iR] = 0 # remove parents
    
  }
  
  iS = which(Gtr[,Yi]!=0)
  
  nb = newton_rhapLR_L0( cbind(X[,iS,drop=FALSE]/C, Xp), Y, cbind(Z[,iS,drop=FALSE], Zp), dp) # train logistic regressor
  
  Gn = matrix(0,d+1,d+1)
  ix = 1:(d+1); ix = ix[-Yi]
  Gn[ix,ix] = G; Gn[ix[iS],Yi] = nb$alpha[seq_len(length(iS))]
  
  #   
  ## LOGISTIC REGRESSION
  # print('predictive model')
  shaps = matrix(0,dp,d); mod = list()
  if (sum(abs(Gn[ix[iS],Yi]))>0){
    require(xgboost)
    Yn = Z %*% Gn[-Yi,Yi,drop=FALSE] + Zp %*% as.matrix(nb$alpha[(length(iS)+1):(length(iS)+dp)])
    
    mod = train_xgboost_class(E[,iS,drop=FALSE],Yn)
    
    Yn1 = logistic(Z %*% Gn[-Yi,Yi,drop=FALSE] + Zp %*% as.matrix(nb$alpha[(length(iS)+1):(length(iS)+dp)]))
    Yn1 = rbinom(length(Yn1),size=1,prob=Yn1)
    
    require(treeshap)
    unified_mod <- xgboost.unify(mod, E[,iS,drop=FALSE])
    tdata = as.data.frame(E[,iS,drop=FALSE])
    # tdata = as.data.frame(data$E[,-DAG$Y])
    colnames(tdata) = 0:(ncol(tdata)-1)
    shap <- as.matrix(treeshap(unified_mod, tdata, verbose = FALSE)$shaps)
    
    iY = which(Yn1==1)
    for (j in seq_len(dp)){
      ix = which(Zp[,j]==1); ix = intersect(ix,iY)
      shaps[j,iS] = colMeans(shap[ix,,drop=FALSE])
    }
  }
  #
  return(list(G=Gn,shapes=shapes,rates=rates, mod = mod, Z=Z, Zp=Zp, E=E, shaps = shaps))
}
