raup.crick <- function(spXsite, reps=999, parallel=TRUE, classic_metric=FALSE, similarity=FALSE,
  split_ties=TRUE, as.distance.matrix=TRUE, set_all_species_equal=FALSE, rownames_in_col1=FALSE) {
  require(wrswoR, quietly=T)
  if(rownames_in_col1){
    row.names(spXsite) <- spXsite[,1]
    spXsite <- spXsite[,-1]
  }
  #make abundances into ocurrences
  spXsite <- ceiling(spXsite / max(spXsite))
  storage.mode(spXsite) <- "integer"
  #create an occurrence vector - used to give frequency weights for null sampling
  gamma <- ncol(spXsite)
  occur <- colSums(spXsite)
  #NOT recommended- this sets all species to occur with equal frequency in the null model
  if(set_all_species_equal){
    occur <- rep(1, gamma)
    warning('This sets all species frequencies as equal and is the same as vegdist(x,method="raup") but much slower.')
  }
  pair_count <- nrow(spXsite) * (nrow(spXsite) - 1) / 2
  alpha_levels <- rowSums(spXsite)
  v_sample_int_crank <- Vectorize(sample_int_crank, vectorize.args="size")
  #build a null distribution of the expected (null) number of shared species for a pair of alpha values:
  if(!parallel) {
    n_shared_null <- matrix(integer(reps * pair_count), ncol=reps, nrow=pair_count)
    for(i in 1:reps){
      cat(paste("\rGenerating bootstrap null values; iteration", i, "of",reps, sep=" "))
      nulls <- v_sample_int_crank(gamma, alpha_levels, occur)
      nulls <- t(sapply(nulls, function(x) {y <- integer(gamma); y[x] <- 1; as.integer(y)}))
      n_shared_null[,i] <- as.integer(as.dist(tcrossprod(nulls)))
      flush.console()
    }
    rm(nulls)
  }
  if(parallel) {
    options(warn=-1)
    require(foreach, quietly=T)
    require(doParallel, quietly=T) #for Windows
    nCores <- detectCores()
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    message("Generating bootstrap null values")
    #foreach() from the foreach package, %dopar% from the doParallel package
    n_shared_null <- foreach(icount(reps), .combine=cbind) %dopar% {
      nulls <- v_sample_int_crank(gamma, alpha_levels, occur)
      nulls <- t(sapply(nulls, function(x) {y <- integer(gamma); y[x] <- 1; as.integer(y)}))
      as.integer(as.dist(tcrossprod(nulls)))
    }
    stopCluster(cl)
    rm(nulls)
    message("complete")
  }
  n_shared_obs <- as.dist(tcrossprod(spXsite))
  n_shared_null <- split(n_shared_null, factor(1:nrow(n_shared_null)))
  #How many time were more species shared in the null distributions than were observed?
  n_greater_than_null <- mapply(function(x,y) {sum(x>y)}, n_shared_null, n_shared_obs, USE.NAMES=F)
  #How many time were the same number of species shared in the null distributions as were observed?
  n_equal_to_null <- mapply(function(x,y) {sum(x==y)}, n_shared_null, n_shared_obs, USE.NAMES=F)
  if(split_ties) {n_equal_to_null <- n_equal_to_null / 2}
  #combine
  rc <- (n_greater_than_null + n_equal_to_null) / reps 
  #final formatting
  rc_mat <- matrix(NA, nrow=nrow(spXsite), ncol=nrow(spXsite), dimnames=list(rownames(spXsite), rownames(spXsite)))
  rc_mat[lower.tri(rc_mat)] <- rc
  rc_mat[upper.tri(rc_mat)] <- rc
  diag(rc_mat) <- 0
  #Adjustments; default is [-1,1] series of dissimilarites 
  #(i.e., -1 indicates sharing more species than any null distribution, 1 indicates sharing less than any null distribution)
  if(!classic_metric) {rc_mat <- (rc_mat - .5) * 2}
  if(similarity & !classic_metric) {rc_mat <- rc_mat * -1}
  if(similarity & classic_metric) {rc_mat <- 1 - rc_mat}
  if(as.distance.matrix) {return(as.dist(rc_mat))} else {return(rc_mat)}
}
