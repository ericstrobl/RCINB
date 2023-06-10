sample_NB_DAG_r <- function(nsamps, DAG){
  
  G = DAG$graph
  r = nrow(G)
  
  # construct Xp
  nsamps_p = floor(nsamps/DAG$dp)
  Xp = matrix(0,nsamps,DAG$dp)
  for (u in seq_len(DAG$dp)){
    if (u != DAG$dp){
      ip = (((u-1)*nsamps_p)+1):(((u-1)*nsamps_p)+nsamps_p) 
    } else{
      ip = (((u-1)*nsamps_p)+1):nsamps
    }
    Xp[ip,u] = 1
  }
  
  Y = DAG$Y
  
  err=matrix(0,nsamps,r)
  for (i in 1:r){
    err[,i]=matrix(rgamma(nsamps,shape=DAG$rs[i],rate=DAG$rs[i]/exp(Xp %*% DAG$offs[,i,drop=FALSE])),nsamps,1)
    # err[,i]=rgamma(nsamps,shape=DAG$rs[i],rate=DAG$rs[i])
  }
  
  err[,Y]=1 #error for diagnosis is zero (exponential of zero is one)
  # data=normalizeData(err) # dont need to normalize data here
  data = err
  
  done=which(colSums(G)==0) # variables without parents
  stop=0;
  while (stop==0){
    for (s in done){
      ch=which(G[s,]==1) 
      for (c in ch){
        if (c %in% done){
          next
        }
        pa=which(G[,c]==1) 
        
        h=intersect(pa,done)
        
        if (setequal(h,pa)){ # if all parents already done
          
          A = (data[,h,drop=FALSE]%*%DAG$weights[h,c,drop=FALSE]) # do not need to include offset because the gamma error terms introduce the offset
          if (c!=Y){
            data[,c]=exp(A)*err[,c]
          } else{
            data[,c]=A
          }
          
          done=unique(c(done, c))
        }
      }
    }
    
    if (length(done) == r){
      stop=1;
    }
  }
  
  C = rgamma(nsamps,1,1)
  # C = rep(1,nsamps)
  for (i in setdiff(1:r,Y)){ # poisson transformation for everything but Y
    data[,i]=rpois(nsamps,data[,i]*C)
  }
  
  Y0 = data[,Y]
  pY = logistic(Y0)
  data[,Y] = rbinom(nsamps,size=1,prob=pY) 
  
  return(list(data=data, E=err, Y0=Y0, C=C, Xp=Xp))
}


logistic <- function(X){
  
  return(1/(1+exp(-X)))
}