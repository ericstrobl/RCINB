organize_shaps <- function(shap,Xp,Y){
  
  d = ncol(shap); dp = ncol(Xp)
  shaps = matrix(0,dp,d)
  iY = which(Y==1)
  for (j in seq_len(dp)){
    ix = which(Xp[,j]==1); ix = intersect(ix,iY)
    shaps[j,] = colMeans(shap[ix,])
  }

return(shaps)
}