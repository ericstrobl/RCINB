load("Hua.RData")
require("independence")

RCINB_res = vector("list",500)
RCI_res = RCINB_res
GRCI_ANM_res = RCINB_res
GRCI_HNM_res = RCINB_res
ICA_res = RCINB_res

Gtr = matrix(0,6,7)
Gtr[1,2]=1;Gtr[2,3]=1;Gtr[3,6]=1; Gtr[2,5]=1; Gtr[2,4]=1
Gtr[,7]=1
n = nrow(data)
p = ncol(data)

for (i in 1:50){
  print(i)
  
  shapsT = RCI_NBC2_r_G(cbind(data,Y),ncol(data)+1,C/length(C),as.matrix(Xp),Gtr=Gtr)$shaps

  ib = sample(1:n,n,replace=TRUE)
  dataft = cbind(data[ib,],Y[ib])
  Ct = C[ib]; Xpt = Xp[ib,]
  
  ptm <- proc.time()
  out <- RCI_NB(dataft,ncol(dataft),Ct/length(Ct),as.matrix(Xpt))
  out$order = seq_len(p)
  RCINB_res[[i]]$time = (proc.time() - ptm)[3]
  RCINB_res[[i]]$RBO = eval_RBO(shapsT,out)
  RCINB_res[[i]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
  RCINB_res[[i]]$shaps = out$shaps
  
  ptm <- proc.time()
  datan = lm.fit(cbind(log(Ct/length(Ct)),as.matrix(Xpt)), dataft[,-(p+1)])$residuals
  out = RCI(datan,dataft[,p+1]); scores = matrix(0,n,p); scores[,out$order] = out$scores; out$order = seq_len(p)
  out$shaps = organize_shaps(scores,Xpt,dataft[,p+1])
  RCI_res[[i]]$time = (proc.time() - ptm)[3]
  RCI_res[[i]]$RBO = eval_RBO(shapsT,out)
  RCI_res[[i]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
  RCI_res[[i]]$shaps = out$shaps
  
  ptm <- proc.time()
  out = GRCI_HNM(dataft[,-(p+1)],dataft[,p+1],Ct/length(Ct),as.matrix(Xpt)); outG = out; timeG = out$timeG
  GRCI_HNM_res[[i]]$time = (proc.time() - ptm)[3]
  GRCI_HNM_res[[i]]$RBO = eval_RBO(shapsT,out)
  GRCI_HNM_res[[i]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
  GRCI_HNM_res[[i]]$shaps = out$shaps
  
  ptm <- proc.time()
  out = GRCI_ANM(dataft[,-(p+1)],dataft[,p+1],Ct/length(Ct),as.matrix(Xpt)); outG = out; timeG = out$timeG
  GRCI_ANM_res[[i]]$time = (proc.time() - ptm)[3]
  GRCI_ANM_res[[i]]$RBO = eval_RBO(shapsT,out)
  GRCI_ANM_res[[i]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
  GRCI_ANM_res[[i]]$shaps = out$shaps
  
  ptm <- proc.time()
  out = ICA_predict(datan,dataft[,p+1])
  out$shaps = organize_shaps(out$scores,Xpt,dataft[,p+1])
  ICA_res[[i]]$time = (proc.time() - ptm)[3]
  ICA_res[[i]]$RBO = eval_RBO(shapsT,out)
  ICA_res[[i]]$RMSE = sqrt(mean((out$shaps-shapsT)^2))
  ICA_res[[i]]$shaps = out$shaps
  
}

save(file="Results_Hua.RData", RCINB_res, RCI_res, GRCI_HNM_res,
     GRCI_ANM_res, ICA_res)
}

nrep = 50
res_mat = matrix(0,nrep,5)
for (t in 1:nrep){

  res_mat[t,1] = RCINB_res[[t]]$RMSE
  res_mat[t,2] = RCI_res[[t]]$RMSE
  res_mat[t,3] = ICA_res[[t]]$RMSE
  res_mat[t,4] = GRCI_HNM_res[[t]]$RMSE
  res_mat[t,5] = GRCI_ANM_res[[t]]$RMSE
  
  # res_mat[t,1] = RCINB_res[[t]]$time
  # res_mat[t,2] = RCI_res[[t]]$time
  # res_mat[t,3] = ICA_res[[t]]$time
  # res_mat[t,4] = GRCI_HNM_res[[t]]$time
  # res_mat[t,5] = GRCI_ANM_res[[t]]$time
}

print(colMeans(res_mat))

