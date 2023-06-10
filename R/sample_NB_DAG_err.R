sample_NB_DAG_err <- function(nsamps, err, G){
  
  G = as.matrix(G)
  r = nrow(G)
  data = err
  
  done=which(colSums(G)==0) # variables without parents
  stop=0;
  while (stop==0){
    for (s in done){
      ch=which(G[s,]!=0) 
      for (c in ch){
        if (c %in% done){
          next
        }
        pa=which(G[,c]!=0) 
        
        h=intersect(pa,done)
        
        if (setequal(h,pa)){ # if all parents already done
          
          A = (data[,h,drop=FALSE]%*%G[h,c,drop=FALSE])
          # A = pmin(A,1) ###
          data[,c]=exp(A)*err[,c] # also multiply by offsets here
          
          done=unique(c(done, c))
        }
      }
    }
    
    if (length(done) == r){
      stop=1;
    }
  }
  
  
  return(data)
}
