RCI_NB <- function(data,Yi,C,Xp,alpha=0.02,Gt=NULL){
  
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
  iP = (d+1):(d+dp); 
  while ( length(iV)>0 ){ # until no targets remain
    iR = c() # indices to remove from iV
    s_max = Inf; cnt = 0
    for (v in iV){
      # print(v)
      
      if (length(iS)>0){
        iSs = sparsify_it(X,iS,X[,v,drop=FALSE],C,Xp,alpha)
      } else{
        iSs = iS
      }
      
      oo = test_NB_L0_mom00rp2sC3_fast_r( cbind(X[,iSs,drop=FALSE]/C, Xp), X[,v,drop=FALSE],
                                        cbind(Z[,iSs,drop=FALSE], Zp), Cb, dp, warmStart[v,c(iSs,iP)]) # regress predictors on one target
      
      if (oo$S < s_max){ # in case no p-values greater than Type I error rate
        vf = v
        shapef = oo$r 
        ratef = oo$r/exp(oo$alpha[(length(iSs)+1):(length(iSs)+dp)])
        alphaf = oo$alpha[seq_len(length(iSs))]
        s_max = oo$S
        iSsf = iSs
      }
      
    }
    
    if (cnt==0){ # if none removed, then remove one with largest p-value
      iR = vf # remove index vf
      shapes[vf] = shapef
      rates[,vf] = ratef
      G[iSsf,vf] = alphaf[seq_len(length(iSsf))]
    }
    
    iV = setdiff(iV, iR) # remove iR from iV
    # print(iR)
    
    # sample Z
    iS = c(iS,iR) # include iR into predictors
    oo = sample_gamma_var(n, (length(iS)-length(iR)+1):length(iS), Z[,iS,drop=FALSE], Zp, Cb, 
                          shapes=shapes[iS],rates=rates[,iS,drop=FALSE],DAG=G[iS,iS,drop=FALSE])
    if (nrow(Z)==n){
      Z = matrix(0,nrow(oo$Z),ncol(Z))
      E = matrix(0,nrow(oo$Z),ncol(Z))
    }
    Z[,iR] = oo$Z; Zp = oo$Zp; Cb = oo$Cb; E[,iR] = oo$E

  }
  
  #sparsify - remove vertices from iS that are NOT correlated with Y
  
  iS = sparsify_it(X,iS,data[,Yi],C,Xp,alpha)

  
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
    
    Yn1 = logistic(Z %*% Gn[-Yi,Yi,drop=FALSE] + Zp %*% as.matrix(nb$alpha[(length(iS)+1):(length(iS)+dp)]))
    Yn1 = rbinom(length(Yn1),size=1,prob=Yn1)
    
    for (j in seq_len(dp)){ # perform TreeShap given patient
      ix = which(Zp[,j]==1)
      mod = train_xgboost_class(E[ix,iS,drop=FALSE],Yn[ix])
    
      require(treeshap)
      unified_mod <- xgboost.unify(mod, E[ix,iS,drop=FALSE])
      tdata = as.data.frame(E[ix,iS,drop=FALSE])
      # tdata = as.data.frame(data$E[,-DAG$Y])
      colnames(tdata) = 0:(ncol(tdata)-1)
      shap <- as.matrix(treeshap(unified_mod, tdata, verbose = FALSE)$shaps)
    
      iY = which(Yn1[ix]==1)
      ix1 = intersect(1:length(ix),iY)
      shaps[j,iS] = colMeans(shap[ix1,,drop=FALSE])
      
    }
  }
  #
  return(list(G=Gn,shapes=shapes,rates=rates, mod = mod, Z=Z, Zp=Zp, E=E, shaps = shaps))
}

sample_gamma_var <- function(n,ix,Z,Zp,Cb,shapes=NULL,rates=NULL,DAG=NULL){
  
  dp = ncol(Zp)
  d = ncol(Z)
  if (n==nrow(Z)){
    ns = colSums(Zp)
    ms = floor(ns^(1.1))
    Zp = matrix(0,0,dp)
    for (p in seq_len(dp)){
      Zpn = matrix(0,ms[p],dp); Zpn[,p]=1
      Zp = rbind(Zp,Zpn)
    }
    Cb = sample(Cb,nrow(Zp),replace=TRUE)
    Z = matrix(0,nrow(Zp),d)
  }
  m = nrow(Zp)
  
  # simulate V
  for (i in ix){
    Z[,i] = rgamma(m,shape= shapes[i] ,rate= Zp %*% rates[,i,drop=FALSE])
  }
  E = Z
  if (!is.null(DAG)){
    Z[,1:d] = sample_NB_DAG_err(m,Z[,1:d,drop=FALSE],DAG)
  }
  
  
  return(list(Z=Z[,ix,drop=FALSE],Zp=Zp,Cb=Cb,E=E[,ix]))
  
}


sparsify_it <- function(X,iS,Y,C,Xp,alpha=0.02){
  
  iT = c()
  for (v in iS){
    p = earth_test(X[,v], Y, cbind(C,Xp))$p.value
    iT = c(iT, p)
  }
  eT = iS[which(iT<(alpha*5))] # < 0.1
  if (length(eT)>0){
    iS = eT
  } else{
    iS = iS[which(iT == min(iT))]
  }
  
  return(iS)
  # return(iT)
  
}

