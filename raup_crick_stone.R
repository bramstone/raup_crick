#weighted random sampling implementation from Efraimidis & Spirakis created by user Dinre on Stackoverflow 
#   https://stackoverflow.com/questions/15113650/faster-weighted-sampling-without-replacement
#PS Efraimidis, PG Spirakis. Weighted random sampling with a reservoir. Information Processing Letters 97: 181-185.
#This provides a much faster way to implement weighted random sampling without replacement
#There is also an Rcpp implementation on the above link from user Julius
wr.sample <- function(data, n, prob){
  
  key <- runif(length(data)) ^ (1 / prob)
  return(data[order(key, decreasing=TRUE)][1:n])
}
#v.wr.sample <- Vectorize(wr.sample,vectorize.args="n") #put above null-creation for-loop
#nulls <- v.wr.sample(1:gamma,alpha_levels,occur) #put in null-creation for-loop in place of v_sample_int_crank

#My implementation of Chase et al.'s Raup Crick code, much is copied but the null distribution method is changed
#JM Chase, NJB Kraft, KG Smith, M Vellend, BD Inouye. 2011. Using null models to disentangle variation in community 
#dissimilarity from variation in alpha-diversity. Ecosphere 2: 1-11.
#
#classic_metric = FALSE standardizes the metric to range from -1 to 1 instead of the original 0 to 1 as per Chase et al.
#similarity = FALSE represents an index of dissimilarity, higher values indicate more dissimilarity; TRUE 
#   multiplies by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper

raup.crick=function(spXsite, reps=9999, parallel=TRUE, classic_metric=FALSE, similarity=FALSE,
  split_ties=TRUE, as.distance.matrix=TRUE, set_all_species_equal=FALSE, row_names) {
  require(wrswoR)
  require(vegan)
  if(is.null(row_names)){
    row.names(spXsite) <- spXsite[,row_names]
    spXsite <- spXsite[,-row_names]
  }
  #make abundances into ocurrences
  spXsite <- ceiling(spXsite/max(spXsite))
  storage.mode(spXsite) <- "integer"
  #create an occurrence vector - used to give frequency weights for null sampling
  gamma <- ncol(spXsite)
  occur <- colSums(spXsite)
  #NOT recommended- this sets all species to occur with equal frequency in the null model
  if(set_all_species_equal){
    occur <- rep(1,gamma)
    warning('This sets all species frequencies as equal and is the same as vegdist(x,method="raup") but much slower.')
  }
  pair_count <- nrow(spXsite)*(nrow(spXsite)-1)/2
  alpha_levels <- rowSums(spXsite)
  v_sample_int_crank <- Vectorize(sample_int_crank,vectorize.args="size")
  #build a null distribution of the expected (null) number of shared species for a pair of alpha values:
  if(!parallel) {
    n_shared_null <- matrix(integer(reps*pair_count),ncol=reps,nrow=pair_count)
    for(i in 1:reps){
      cat(paste("\rNull bootstrap iteration",i,"of",reps,sep=" "))
      nulls <- v_sample_int_crank(gamma,alpha_levels,occur)
      nulls <- t(sapply(nulls,function(x) {y <- integer(gamma); y[x] <- 1; as.integer(y)}))
      n_shared_null[,i] <- designdist(nulls,method="J",terms="binary")
      flush.console()
    }
    rm(nulls)
  }
  if(parallel) {
    options(warn=-1)
    require(foreach)
    require(doParallel) #for Windows
    nCores <- detectCores()
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    message("Generating bootstrap null values")
    #foreach() from the foreach package, %dopar% from the doParallel package
    n_shared_null <- foreach(icount(reps),.combine=cbind) %dopar% {
      nulls <- v_sample_int_crank(gamma,alpha_levels,occur)
      nulls <- t(sapply(nulls,function(x) {y <- integer(gamma); y[x] <- 1; as.integer(y)}))
      as.integer(vegan::designdist(nulls,method="J",terms="binary"))
    }
    stopCluster(cl)
    rm(nulls)
    message("complete")
  }
  n_shared_obs <- designdist(spXsite,method="J",terms="binary")
  n_shared_null <- split(n_shared_null,factor(1:nrow(n_shared_null)))
  #How many time were more species shared in the null distributions than were observed?
  n_greater_than_null <- mapply(function(x,y){sum(x>y)}, n_shared_null, n_shared_obs,USE.NAMES=F)
  n_equal_to_null <- mapply(function(x,y){sum(x==y)}, n_shared_null, n_shared_obs,USE.NAMES=F)
  if(split_ties) {n_equal_to_null <- n_equal_to_null/2}
  #combine
  rc <- (n_greater_than_null+n_equal_to_null)/reps 
  #final formatting
  rc_mat <- matrix(NA, nrow=nrow(spXsite), ncol=nrow(spXsite), dimnames=list(rownames(spXsite),rownames(spXsite)))
  rc_mat[lower.tri(rc_mat)] <- rc
  rc_mat[upper.tri(rc_mat)] <- rc
  diag(rc_mat) <- 0
  #adjustments
  if(!classic_metric) {rc_mat <- (rc_mat-.5)*2}
  if(similarity & !classic_metric) {rc_mat <- rc_mat*-1}
  if(similarity & classic_metric) {rc_mat <- 1-rc_mat}
  if(as.distance.matrix) {return(as.dist(rc_mat))} else {return(rc_mat)}
}