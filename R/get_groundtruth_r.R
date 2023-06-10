get_groundtruth_r <- function(data,DAG){
  
  id = c()
  for (d in 1:ncol(data$data)){
    if (d != DAG$Y){
      if (isAnc_fast_LE(DAG$graph,d,DAG$Y)){
        id = c(id, d)
      }
    }
  }
  
  shapsT = matrix(0,ncol(data$Xp),ncol(data$data))
  for (j in seq_len(ncol(data$Xp))){
    ix = which(data$Xp[,j]==1);
    
    mod = train_xgboost_class(data$E[ix,id,drop=FALSE],data$Y0[ix])
    require(treeshap)
    unified_mod <- xgboost.unify(mod, data$E[ix,id,drop=FALSE])
    tdata = as.data.frame(data$E[ix,id,drop=FALSE])
    colnames(tdata) = 0:(ncol(tdata)-1)
    shapT <- as.matrix(treeshap(unified_mod, tdata, verbose = FALSE)$shaps)
    
    iY = which(data$data[ix,DAG$Y]==1)
    
    ix1 = intersect(1:length(ix),iY)
    
    shapsT[j,id] = colMeans(shapT[ix1,,drop=FALSE])
  }
  
  return(shapsT[,-DAG$Y,drop=FALSE])
}