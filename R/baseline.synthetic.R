#' estimate basline copy numbers by synthetic normal cell using GMM to estimate sigmas
#'
#' @param norm.mat smoothed data matrix; genes in rows; cell names in columns.
#' @param min.cells minimal number of cells per cluster.
#' @param n.cores number of cores for parallel computing.
#'
#' @return 1) relative gene expression; 2) synthetic baseline profiles; 3) clustering results.
#'
#' @examples

#' test.relt <- baseline.synthetic(norm.mat=orm.mat.smooth, min.cells=10, n.cores)
#'
#' norm.mat.relat <- test.relt$expr.relat
#' @export


baseline.synthetic <- function(norm.mat=norm.mat, min.cells=10, n.cores){ 

 d <- parallelDist::parDist(t(norm.mat), threads = n.cores) ##use smooth and segmented data to detect intra-normal cells
  km <- 6
  fit <- hclust(d, method="ward.D2")
  ct <- cutree(fit, k=km)

  while(!all(table(ct)>min.cells)){
    km <- km -1
    ct <- cutree(fit, k=km)
    if(km==2){
      break
    }
  }

  
  expr.relat <- NULL
  syn <- NULL
  for(i in min(ct):max(ct)){
    data.c1 <- norm.mat[, which(ct==i)]
    sd1 <- apply(data.c1,1,sd)
    set.seed(123)
    syn.norm <- sapply(sd1,function(x)(x<- rnorm(1,mean = 0,sd=x)))
    relat1 <- data.c1 -syn.norm
    expr.relat <- rbind(expr.relat, t(relat1))
    syn <- cbind(syn,syn.norm)
    i <- i+1
  }

  reslt <- list(data.frame(t(expr.relat)), data.frame(syn), ct)
  names(reslt) <- c("expr.relat","syn.normal", "cl")
  
  return(reslt)
}
