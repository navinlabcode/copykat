#' find a cluster of diploid cells with integrative clustering method
#'
#' @param norm.mat.smooth smoothed data matrix; genes in rows; cell names in columns.
#' @param min.cells minimal number of cells per cluster.
#' @param n.cores number of cores for parallel computing.
#'
#' @return 1) predefined diploid cell names; 2) clustering results; 3) inferred baseline.
#'
#' @examples

#' test.bnc <- baseline.norm.cl(norm.mat.smooth=norm.mat.smooth, min.cells=5, n.cores=10)
#'
#' test.bnc.cells <- test.bnc$preN
#' @export


baseline.norm.cl <- function(norm.mat.smooth, min.cells=5, n.cores=n.cores){

  d <- parallelDist::parDist(t(norm.mat.smooth), threads = n.cores) ##use smooth and segmented data to detect intra-normal cells
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

  SDM <-NULL
  SSD <-NULL
  for(i in min(ct):max(ct)){

    data.c <- apply(norm.mat.smooth[, which(ct==i)],1, median)
    sx <- max(c(0.05, 0.5*sd(data.c)))
    GM3 <- mixtools::normalmixEM(data.c, lambda = rep(1,3)/3, mu = c(-0.2, 0, 0.2), sigma = sx,arbvar=FALSE,ECM=FALSE,maxit=5000)
    SDM <- c(SDM, GM3$sigma[1])
    SSD <- c(SSD, sd(data.c))
       i <- i+1
      }

  wn <- mean(cluster::silhouette(cutree(fit, k=2), d)[, "sil_width"])

  ####
 PDt <- pf(max(SDM)^2/min(SDM)^2, nrow(norm.mat.smooth), nrow(norm.mat.smooth), lower.tail = FALSE)
  #PDt <- dt((min(SDM)-max(SDM))/mad(SDM),df=km-1)

 # print(c("low sigma pvalue:", PDt))
  #print(c("low sd pvalue:", dt((min(SSD)-max(SSD))/mad(SSD),df=km-1)))

  if(wn <= 0.15|(!all(table(ct)>min.cells))| PDt > 0.05){
    WNS <- "unclassified.prediction"
    print("low confidence in classification")
  }else {
    WNS <- ""
  }
    basel <- apply(norm.mat.smooth[, which(ct %in% which(SDM==min(SDM)))], 1, median)
    preN <- colnames(norm.mat.smooth)[which(ct %in% which(SDM==min(SDM)))]

  ### return both baseline and warning message
  RE <- list(basel, WNS, preN, ct)
  names(RE) <- c("basel", "WNS", "preN", "cl")
  return(RE)
    }


