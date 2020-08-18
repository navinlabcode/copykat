#' pre-define a group of normal cells with GMM.
#'
#' @param CNA.mat smoothed data matrix; genes in rows; cell names in columns.
#' @param max.normal find the first number diploid cells to save efforts.
#' @param mu.cut diploid baseline cutoff.
#' @param Nfraq.cut minimal fractoins of genomes with CNAs.
#'
#' @return 1) predefined diploid cell names; 2) clustering results; 3) inferred baseline.
#'
#' @examples
#' test.gmm <- baseline.GMM(CNA.mat=smooth.com, max.normal=30, mu.cut=0.05, Nfraq.cut=0.99)
#'
#' test.gmm.cells <- test.bnc$preN
#' @export
baseline.GMM <- function(CNA.mat, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99, RE.before=basa, n.cores=1){

     N.normal <-NULL
     for(m in 1:ncol(CNA.mat)){

      print(paste("cell: ", m, sep=""))
      sam <- CNA.mat[, m]
      sg <- max(c(0.05, 0.5*sd(sam)));
      GM3 <- mixtools::normalmixEM(sam, lambda = rep(1,3)/3, mu = c(-0.2, 0, 0.2),sigma = sg,arbvar=FALSE,ECM=FALSE,maxit=500);#maxrestarts=10; arbmean=TRUE; arbvar=TRUE;epsilon=0.01

      ###decide normal or tumor
      s <- sum(abs(GM3$mu)<=mu.cut)

      if(s>=1){
        frq <- sum(GM3$lambda[which(abs(GM3$mu)<=mu.cut)])
     #   print(paste("N.fraq ", frq, sep=""))
     #   print(paste("sigma: ", GM3$sigma[1], sep=""))
        if(frq> Nfraq.cut){
          pred <- "diploid"
        }else{pred<-"aneuploid"}

      }else {pred <- "aneuploid"}
    #  print(paste("pred: ", pred, sep=""))
       N.normal<- c(N.normal,pred)

      if(sum(N.normal=="diploid")>=max.normal){break}
      m<- m+1
    }

    names(N.normal) <- colnames(CNA.mat)[1:length(N.normal)]
    preN <- names(N.normal)[which(N.normal=="diploid")]

    d <- parallelDist::parDist(t(CNA.mat), threads = n.cores) ##use smooth and segmented data to detect intra-normal cells
    km <- 6
    fit <- hclust(d, method="ward.D2")
    ct <- cutree(fit, k=km)


    if(length(preN) >2){
      WNS <- ""
      basel <- apply(CNA.mat[, which(colnames(CNA.mat) %in% preN)], 1, mean)

      RE <- list(basel, WNS, preN, ct)
      names(RE) <- c("basel", "WNS", "preN", "cl")
      return(RE)
    }else{
      return(RE.before) ##found this bug
    }

}

