#' MCMC segmentation.
#'
#' @param fttmat transformed normalized data matrix; genes in rows; cell names in columns.
#' @param clu predefined clusters of cells.
#' @param sam.name given sample name.
#' @param bins minimal number of bins per segment.
#' @param cut.cor minimal KS difference to call breakpoints.
#' @param n.cores number of cores for parallel computing.
#'
#' @return 1) segmented CNA data matrix; 2) vector of breakpoints.
#'
#' @examples
#' test.mcmc <- CNA.MCMC(clu,fttmat, bins, cut.cor, n.cores=10)
#'
#' test.mcmc.seg.mat <- test.mamc$logCNA
#' @export
CNA.MCMC <- function(clu,fttmat, bins, cut.cor, n.cores){
  CON<- NULL
  for(i in min(clu):max(clu)){
    data.c <- apply(fttmat[, which(clu==i)],1, median)
    CON <- cbind(CON, data.c)
    i <- i+1
  }

  norm.mat.sm <- exp(CON)
  n <- nrow(norm.mat.sm)

  BR <- NULL

  for(c in 1:ncol(norm.mat.sm)){

    breks <- c(seq(1, as.integer(n/bins-1)*bins, bins),n)
    bre <- NULL

    for (i in 1:(length(breks)-2)){
      #i<-42
      a1<-  max(mean(norm.mat.sm[breks[i]:breks[i+1],c]), 0.001)
      posterior1 <-MCMCpack::MCpoissongamma(norm.mat.sm[breks[i]:breks[i+1],c], a1, 1, mc=1000)


      a2 <- max(mean(norm.mat.sm[(breks[i+1]+1):breks[i+2],c]), 0.001)
      posterior2 <-MCMCpack::MCpoissongamma(norm.mat.sm[(breks[i+1]+1):breks[i+2],c], a2, 1, mc=1000)

      if (ks.test(posterior1,posterior2)$statistic[[1]] > cut.cor){
        bre <- c(bre, breks[i+1])
       }

        i<- i+1
        }

    breks <- sort(unique(c(1, bre, n)))
    BR <- sort(unique(c(BR, breks)))
    c<-c+1
  }

  #print(paste(length(BR), " breakpoints", sep=""))

  ###CNA
  norm.mat.sm <- exp(fttmat)

  seg <- function(z){
      x<-numeric(n)
      for (i in 1:(length(BR)-1)){
        a<- max(mean(norm.mat.sm[BR[i]:BR[i+1],z]), 0.001)
        posterior1 <-MCMCpack::MCpoissongamma(norm.mat.sm[BR[i]:BR[i+1],z], a, 1, mc=1000)
        x[BR[i]:BR[i+1]]<-mean(posterior1)
         i<- i+1
      }
      x<-log(x)

  }

  seg.test <- parallel::mclapply(1:ncol(norm.mat.sm), seg, mc.cores = n.cores)
  logCNA <- matrix(unlist(seg.test), ncol = ncol(norm.mat.sm), byrow = FALSE)

  res <- list(logCNA, BR)
  names(res) <- c("logCNA","breaks")
  return(res)
}

