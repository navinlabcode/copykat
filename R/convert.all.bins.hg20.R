#' conver gene by cell matrix to genomic bins by cells matrix.
#'
#' @param DNA.mat input target bins, provided by copycat, with 220KB windows.
#' @param RNA.mat RNA data matrix with genes in rows.
#' @param full.anno annotation file for all known genes, automatically loaded in copycat.
#' @param n.cores number of cores for parallel computing.
#'
#' @return adjusted datamatrix with bins in rows, cells in columns.
#'
#' @examples
#' test.cab <- convert.all.bins(DNA.mat, RNA.mat, n.cores=10)
#'
#' test.cab.uber <- test.cab$RNAadj
#' @export
convert.all.bins.hg20 <- function(DNA.mat, RNA.mat, n.cores){
##make list obj for each window
         DNA <- DNA.mat[-which(DNA.mat$chrom==24),]; dim(DNA)
         end <- DNA$chrompos
         start <- c(0, end[-length(end)])
          ls.all <- list()
          for(i in 1:nrow(DNA)){
          sub.anno <- full.anno[which(full.anno$chromosome_name==DNA$chrom[i]),]
          cent.gene <- 0.5*(sub.anno$start_position+sub.anno$end_position)
          x <- sub.anno$hgnc_symbol[which(cent.gene<=end[i] & cent.gene>= start[i])]
          if(length(x)==0){x <- "NA"}
          ls.all[[i]] <- x
          i<- i+1
          }

          ##convert gene to bin
          RNA <- RNA.mat[, 8:ncol(RNA.mat)]

         ###adj
          R.ADJ <- function(i){
            shr <- intersect(ls.all[[i]], RNA.mat$hgnc_symbol)
            if(length(shr)>0){
              Aj <- apply(RNA[which(RNA.mat$hgnc_symbol %in% shr), ], 2, median)
            }

          }

          test.adj <- parallel::mclapply(1:nrow(DNA), R.ADJ, mc.cores = n.cores)
          RNA.aj <- matrix(unlist(test.adj), ncol = ncol(RNA), byrow = TRUE)
          colnames(RNA.aj) <- colnames(RNA.mat)[8:ncol(RNA.mat)]

          mind <- which(test.adj=="NULL")

           if(length(mind)>1){
           ind <- 1:nrow(DNA)
           Rw <- ind[-which(ind %in% mind)]

           FK.again <- function(i){
            fkI <- abs(Rw-mind[i])
            fk <-  RNA.aj[which(fkI==min(fkI))[1], ]
          }

          tt.FK <-  parallel::mclapply(1:length(mind), FK.again, mc.cores = n.cores)
          FK <- matrix(unlist(tt.FK), ncol = ncol(RNA), byrow = TRUE)
          colnames(FK) <- colnames(RNA.mat)[8:ncol(RNA.mat)]

           RNA.aj <- cbind(Rw, RNA.aj)
           FK <- cbind(mind, FK)
           RNA.co <- data.frame(rbind(RNA.aj, FK))

          RNA.com <- RNA.co[order(RNA.co$Rw, decreasing = FALSE), ]
          RNA.adj <- cbind(DNA[, 1:3], RNA.com[, 2:ncol(RNA.com)])
         } else{
           RNA.adj <- cbind(DNA[, 1:3], RNA.aj)
         }

        reslt <- list(DNA, RNA.adj)
        names(reslt) <- c("DNA.adj", "RNA.adj")
        return(reslt)

        }


