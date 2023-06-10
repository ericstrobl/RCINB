generate_NB_DAG_r <- function(p,en){
  
  N = p*p - p;
  
  DAG = list()
  
  samplesB = rbinom(N/2,1, en/(p-1) ); # sample edges
  graph = matrix(0,p,p)
  graph[upper.tri(graph, diag=FALSE)] <- samplesB; # put in edges in upper triangular
  
  ord = sample(1:p,p,replace=FALSE) # permute order
  DAG$graph = graph[ord,ord]
  
  Ys = which( (rowSums(DAG$graph)==0) & (colSums(DAG$graph)>0) ) # no children, some observed parents
  DAG$Y = Ys[sample(length(Ys),1)]
  
  weights = matrix(-0.75*runif(p^2)-0.25,p,p)
  DAG$weights = weights*DAG$graph
  DAG$weights[,DAG$Y] = DAG$weights[,DAG$Y]*sample(c(-1,1),p,replace=TRUE) # target can have negative and positive coefficients

  DAG$dp = sample(2:10,1)
  DAG$rs = runif(p,0.1,1) ###
  DAG$offs = matrix(-0.75*runif(DAG$dp*p)-0.25,DAG$dp,p)
  
  return(DAG)
}
