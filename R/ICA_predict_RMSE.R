ICA_predict_RMSE <- function(X,Y){
  require(RcppHungarian)
  require(ica)
  require(randomForest)
  
  c = ncol(X)
  time = proc.time()
  mm = icafast(X, c)
  time = (proc.time() - time)[3]
  
  id = HungarianSolver(t(1/abs(mm$M)))$pairs[,2]
  E = mm$S[,order(id)]
  
  # mm$W = t(mm$W)
  # id = c()
  # for (i in 1:c){
  #   Mi = abs(mm$W[,i]) # what is the variable in X that contributes most to source S[,i]
  #   id  = c(id,which(Mi == max(Mi))[1])
  # }
  # E = mm$S[,order(id)]
  
  ## patient-specific root causes
  beta = glm.fit(cbind(E,1),Y,family=binomial())$coefficients[1:ncol(E)]
  # print(beta)
  if (length(beta)==1){
    scores = E * beta
  } else{
    scores = E %*% diag(beta)
  }
  
  return(list(scores = scores, order = 1:c, E = E, time = time))
  
}