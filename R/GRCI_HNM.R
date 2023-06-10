GRCI_HNM <- function(X,Y,C,Xp,outL=NULL){
  require(RANN)
  require(gtools)
  
  p = ncol(X)
  if (is.null(outL)){
    ptm <- proc.time()
    # X = normalizeData(X)
    
    ## SKELETON DISCOVERY
    # print('skeleton')
    # require(pcalg)
    # suffStat = list()
    # suffStat$data = X
    # G = pc(suffStat, earth_wrap, alpha=0.10, p=ncol(X))
    # G = as(G@graph, "matrix")
    # G = ((G + t(G))>0)
    G = matrix(TRUE,p,p)
    diag(G) = FALSE
    
    ## EXTRACT ERRORS
    # print('errors')
    
    outL = DirectHNM_fast_Y_HNM(X,Y,C,Xp,G)
    timeG = (proc.time() - ptm)[3]
  }
  
  
  ## LOGISTIC REGRESSION
  # print('predictive model')
  beta = glm.fit(cbind(X,Xp),Y,family=binomial())$coefficients[1:ncol(X)]
  Yn = X %*% as.matrix(beta)
  E = outL$X
  colnames(E) = NULL
  # Yn = logistic(Z %*% Gn[-Yi,Yi,drop=FALSE] + Zp %*% as.matrix(nb$alpha[(length(iS)+1):(length(iS)+dp)]))
  # Yn = rbinom(length(Yn),size=1,prob=Yn)
  mod = train_xgboost_class(E,Yn)
  
  require(treeshap)
  unified_mod <- xgboost.unify(mod, E)
  tdata = as.data.frame(E)
  # tdata = as.data.frame(data$E[,-DAG$Y])
  colnames(tdata) = 0:(ncol(tdata)-1)
  shap <- as.matrix(treeshap(unified_mod, tdata, verbose = FALSE)$shaps)
  
  dp = ncol(Xp)
  shaps = matrix(0,dp,ncol(X))
  for (j in seq_len(dp)){
    ix = which(Xp[,j]==1); ix = intersect(ix,which(Y==1))
    shaps[j,] = colMeans(shap[ix,])
  }
  
  
  return(list(E=outL$X, order=1:ncol(X), shaps=shaps, G=G, outL = outL, timeG = timeG))
}

DirectHNM_fast_Y_HNM <- function(X,Y,C,Xp,G){
  X = as.matrix(X)
  # E = X
  
  K = c()
  S = rep(-Inf,ncol(X))
  U = 1:ncol(X)
  update = U
  
  # penalty = 1
  # w = rdirichlet(200,rep(1,ncol(X))*penalty) # random projections
  # # w = normalize_rep(w)
  # w = t(w)
  
  repeat{ 
    s_out = FindSink_HNM(X,U,S,C,Xp,G,update)
    sink = s_out$sink
    S = s_out$S
    
    K = c(K,sink)
    U = U[-which(U==sink)]
    
    if (length(U)==0){ ###
      break ###
    } ###
    
    if (sum(G[sink,])){
      # X[,sink] = PartialOut(X[,G[sink,]],X[,sink],w[G[sink,],,drop=FALSE]) # partial out neighbors from sink
      X[,sink] = nb_Pearson_HNM(X[,G[sink,],drop=FALSE],X[,sink],C,Xp)
    }
    update = intersect(U,which(G[sink,]))
    G[sink,] = FALSE; G[,sink]=FALSE # remove node
    
  }
  
  return(list(K=K,X=X)) #output
}


FindSink_HNM <- function(X,U,S,C,Xp,G,update,w=NULL){
  r = length(update)
  
  ## BASELINE SCORES
  for (i in seq_len(r)){
    V = which(G[update[i],])
    if (length(V)==0){ # if no neighbors, then break because score is -Inf
      next
    }
    
    S[update[i]] = CompareG_HNM(X,C,Xp,update[i],V,w) # baseline
    
  }
  
  sink = U[S[U]==min(S[U])][1] 
  
  return(list(sink=sink,S=S)) #output
  
}

CompareG_HNM <- function(X,C,Xp,i,js,w=NULL){
  
  if (identical(i,js)){  ####
    return(0) ####
  }
  
  # rXY = PartialOut(X[,js],X[,i],w[js,,drop=FALSE])
  # rXY = lm.fit(X[,js,drop=FALSE],X[,i])$residuals
  rXY = nb_Pearson_HNM(X[,js,drop=FALSE],X[,i],C,Xp)
  
  score = c()
  for (k in 1:length(js)){
    score = c(score, hoeffding.D.test(X[,js[k]]+rnorm(length(X[,js[k]]))*1E-6,
                                      c(rXY)+rnorm(length(rXY))*1E-6,precision=1)$scale)
  }
  score = max(score)
  
  return(score) #13  ####
}


nb_Pearson_HNM <- function(X,Y,C,Xp){
  
  X = as.matrix(X)
  dp = ncol(Xp); d = ncol(X)
  warmStart = rep(0,d+dp)
  nb = newton_rhapNB_L0rp_r(cbind(X,Xp),Y,cbind(X,Xp),C,dp,warmStart)
  # nb = newton_rhapNB_L0rp(cbind(X,Xp),Y,cbind(X,Xp),C,dp,warmStart)
  
  mu = exp(cbind(X,Xp) %*% nb$alpha)*C
  # r = exp(Xp %*% as.matrix(nb$alpha[(d+1):(d+dp)]));
  r = nb$r
  Var = mu+(1/r)*mu^2
  resid = (Y-mu)/sqrt(Var) # pearson residual
  
  return(resid)
  
}