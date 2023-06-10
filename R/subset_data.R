subset_data <- function(data,DAG,new_samps){
  
  datan = data
  nsamps = nrow(data$data)
  
  # construct Xp
  nsamps_p = floor(nsamps/DAG$dp)
  new_samps_p = floor(new_samps/DAG$dp)
  resid = new_samps - new_samps_p*DAG$dp
  Xp = matrix(0,nsamps,DAG$dp)
  ips = c()
  for (u in seq_len(DAG$dp)){
    if (u != DAG$dp){
      ip = (((u-1)*nsamps_p)+1):(((u-1)*nsamps_p)+new_samps_p)
    } else{
      ip = (((u-1)*nsamps_p)+1):(((u-1)*nsamps_p)+new_samps_p+resid)
    }
    ips = c(ips, ip)
  }
  
  datan = data
  datan$Xp = data$Xp[ips,,drop=FALSE]
  datan$E = data$E[ips,,drop=FALSE]
  datan$data = data$data[ips,,drop=FALSE]
  datan$Y0 = data$Y0[ips]
  datan$C = data$C[ips]
  
  return(datan)
}