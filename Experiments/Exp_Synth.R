require(pcalg)
require(xgboost)
require(rootSolve)
require(independence)
require(treeshap)

ns = c(10000,100000)
ps = c(8,13)

Gs = lapply(1:reps, function (.) lapply(1:length(ps), function(.) vector("list",length(ns))))
RCINB_res = Gs
RCI_res = Gs
GRCI_ANM_res = Gs
GRCI_HNM_res = Gs
ICA_res = Gs

for (i in 1:250){
  for (p in 1:length(ps)){
    
    DAG = generate_NB_DAG_r(ps[p],2)
    print(DAG$Y)
    plot(as(DAG$graph,"graphNEL"))
    dataAll = sample_NB_DAG_r(ns[length(ns)], DAG)
    shapsT = get_groundtruth_r(dataAll,DAG)
    
    for (n in 1:length(ns)){

      if (n!=length(ns)){
        data = subset_data(dataAll,DAG,ns[n])
      } else{
        data = dataAll
      }
      # ground truth
      
      Gs[[i]][[p]][[n]]$shapsT = shapsT
      Gs[[i]][[p]][[n]]$DAG = DAG
      
      attempt = 0
      out = NULL
      while( is.null(out) && attempt <= 3 ) {
        attempt = attempt + 1
        ptm <- proc.time()
        try(out <- RCI_NB(data$data,DAG$Y,data$C,data$Xp),silent=TRUE)
      }
      out$order = seq_len(ps[p]-1)
      RCINB_res[[i]][[p]][[n]]$time = (proc.time() - ptm)[3]
      RCINB_res[[i]][[p]][[n]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
      RCINB_res[[i]][[p]][[n]]$shaps = out$shaps
      
      ptm <- proc.time()
      datan = lm.fit(cbind(log(data$C),data$Xp), data$data[,-DAG$Y])$residuals
      out = RCI(datan,data$data[,DAG$Y]); scores = matrix(0,ns[n],ps[p]-1); scores[,out$order] = out$scores; out$order = seq_len(ps[p]-1)
      out$shaps = organize_shaps(scores,data$Xp,data$data[,DAG$Y])
      RCI_res[[i]][[p]][[n]]$time = (proc.time() - ptm)[3]
      RCI_res[[i]][[p]][[n]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
      RCI_res[[i]][[p]][[n]]$shaps = out$shaps
      

      ptm <- proc.time()
      out = GRCI_HNM(data$data[,-DAG$Y],data$data[,DAG$Y],data$C,data$Xp); outG = out; timeG = out$timeG
      GRCI_HNM_res[[i]][[p]][[n]]$time = (proc.time() - ptm)[3]
      GRCI_HNM_res[[i]][[p]][[n]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
      GRCI_HNM_res[[i]][[p]][[n]]$shaps = out$shaps
      
      ptm <- proc.time()
      out = GRCI_ANM(data$data[,-DAG$Y],data$data[,DAG$Y],data$C,data$Xp); outG = out; timeG = out$timeG
      GRCI_ANM_res[[i]][[p]][[n]]$time = (proc.time() - ptm)[3]
      GRCI_ANM_res[[i]][[p]][[n]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
      GRCI_ANM_res[[i]][[p]][[n]]$shaps = out$shaps

      ptm <- proc.time()
      out = ICA_predict_RMSE(datan,data$data[,DAG$Y])
      out$shaps = organize_shaps(out$scores,data$Xp,data$data[,DAG$Y])
      ICA_res[[i]][[p]][[n]]$time = (proc.time() - ptm)[3]
      ICA_res[[i]][[p]][[n]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
      ICA_res[[i]][[p]][[n]]$shaps = out$shaps

      
      save(file="Results_synth.RData", Gs, RCINB_res, RCI_res, 
           ICA_res, GRCI_ANM_res, GRCI_HNM_res)
    }
  }
}


## RBO
imax = 250
RBO_RCINBm = array(0,c(imax,length(ps),length(ns)))
RBO_RCIm = RBO_RCINBm
RBO_HNMm = RBO_RCINBm
RBO_ICAm = RBO_RCINBm
RBO_ANMm = RBO_RCINBm

for (i in 1:imax){
  for (p in 1:length(ps)){
    for (n in 1:length(ns)){
      
      RBO_RCINBm[i,p,n] = RCINB_res[[i]][[p]][[n]]$RMSE
      RBO_RCIm[i,p,n] = RCI_res[[i]][[p]][[n]]$RMSE
      RBO_ICAm[i,p,n] = ICA_res[[i]][[p]][[n]]$RMSE
      RBO_HNMm[i,p,n] = GRCI_HNM_res[[i]][[p]][[n]]$RMSE
      RBO_ANMm[i,p,n] = GRCI_ANM_res[[i]][[p]][[n]]$RMSE

      # RBO_RCINBm[i,p,n] = RCINB_res[[i]][[p]][[n]]$time
      # RBO_RCIm[i,p,n] = RCI_res[[i]][[p]][[n]]$time
      # RBO_ICAm[i,p,n] = ICA_res[[i]][[p]][[n]]$time
      # RBO_HNMm[i,p,n] = GRCI_HNM_res[[i]][[p]][[n]]$time
      # RBO_ANMm[i,p,n] = GRCI_ANM_res[[i]][[p]][[n]]$time
    }
  }
}

for (p in 1:length(ps)){
  for (n in 1:length(ns)){
    print(mean(RBO_RCINBm[1:imax,p,n],na.rm=TRUE))
  }
}
