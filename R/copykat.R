#' copycat main_func.
#'
#' @param rawmat raw data matrix; genes in rows; cell names in columns.
#' @param  id.type gene id type: Symbol or Ensemble.
#' @param  cell.line if the data are from pure cell line,put "yes"; if cellline data are a mixture of tumor and normal cells, still put "no".
#' @param LOW.DR minimal population fractoins of genes for smoothing.
#' @param UP.DR minimal population fractoins of genes for segmentation.
#' @param win.size minimal window sizes for segmentation.
#' @param norm.cell.names a vector of normal cell names.
#' @param KS.cut segmentation parameters, input 0 to 1; larger looser criteria.
#' @param sam.name sample name.
#' @param n.cores number of cores for parallel computing.
#' @param ngene.chr minimal number of genes per chromosome for cell filtering.
#' @param distance  distance methods include euclidean, and correlation coverted distance include pearson and spearman.
#' @return 1) aneuploid/diploid prediction results; 2) CNA results in 220KB windows; 3) heatmap; 4) hclustering object.
#'
#' @examples
#' test.ck <- copykat(rawmat=rawdata, sam.name="test", n.cores=10)
#'
#' test.pred <- test.ck$prediction
#' @export


copykat <- function(rawmat=rawdata, id.type="S", cell.line="no", ngene.chr=5,LOW.DR=0.05, UP.DR=0.1, win.size=25, norm.cell.names="", KS.cut=0.1, sam.name="", distance="euclidean", n.cores=1){
  start_time <- Sys.time()
  set.seed(1)
  sample.name <- paste(sam.name,"_copykat_", sep="")

  print("running copykat v1.0.4")
  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep=""))

  genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))

  if(sum(genes.raw> 200)==0) stop("none cells have more than 200 genes")
  if(sum(genes.raw<100)>1){
    rawmat <- rawmat[, -which(genes.raw< 200)]
    print(paste("filtered out ", sum(genes.raw<=200), " cells with less than 200 genes; remaining ", ncol(rawmat), " cells", sep=""))
  }
  ##
  der<- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)

  if(sum(der>LOW.DR)>=1){
    rawmat <- rawmat[which(der > LOW.DR), ]; print(paste(nrow(rawmat)," genes past LOW.DR filtering", sep=""))
  }

  WNS1 <- "data quality is ok"
  if(nrow(rawmat) < 7000){
    WNS1 <- "low data quality"
    UP.DR<- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }

  print("step 2: annotations gene coordinates ...")
  anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  anno.mat <- anno.mat[order(anno.mat$abspos, decreasing = FALSE),]

# print(paste(nrow(anno.mat)," genes annotated", sep=""))

  ### module 3 removing genes that are involved in cell cycling
  HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
  toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
  if(length(toRev)>0){
    anno.mat <- anno.mat[-toRev, ]
  }

#  print(paste(nrow(anno.mat)," genes after rm cell cycle genes", sep=""))
  ### secondary filtering
  ToRemov2 <- NULL
  for(i in 8:ncol(anno.mat)){
    cell <- cbind(anno.mat$chromosome_name, anno.mat[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    } else if(length(rle(cell[,1])$length)<23|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i<- i+1
  }

  if(length(ToRemov2)==(ncol(anno.mat)-7)) stop("all cells are filtered")

  if(length(ToRemov2)>0){
    anno.mat <-anno.mat[, -which(colnames(anno.mat) %in% ToRemov2)]
  }

 # print(paste("filtered out ", length(ToRemov2), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat<- log(sqrt(rawmat3)+sqrt(rawmat3+1))
  norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))
  colnames(norm.mat) <-  colnames(rawmat3)

  #print(paste("A total of ", ncol(norm.mat), " cells, ", nrow(norm.mat), " genes after preprocessing", sep=""))

  ##smooth data
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c){
    model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x<- x[2:length(x)]
    x <- x-mean(x)
  }

  test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)

  colnames(norm.mat.smooth) <- colnames(norm.mat)

  print("step 4: measuring baselines ...")
  if (cell.line=="yes"){
  	print("running pure cell line mode")
  	    relt <- baseline.synthetic(norm.mat=norm.mat.smooth, min.cells=10, n.cores=n.cores)
		norm.mat.relat <- relt$expr.relat
		CL <- relt$cl
        WNS <- "run with cell line mode"
    	preN <- NULL

  } else if(length(norm.cell.names)>1){

    #print(paste(length(norm.cell.names), "normal cells provided", sep=""))
    NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% norm.cell.names)])
    print(paste(NNN, " known normal cells found in dataset", sep=""))

    if (NNN==0) stop("known normal cells provided; however none existing in testing dataset")
    print("run with known normal...")

    basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)],1,median); print("baseline is from known input")

    d <- parallelDist::parDist(t(norm.mat.smooth),threads =n.cores, method="euclidean") ##use smooth and segmented data to detect intra-normal cells

    km <- 6
    fit <- hclust(d, method="ward.D2")
    CL <- cutree(fit, km)

    while(!all(table(CL)>5)){
      km <- km -1
      CL <- cutree(fit, k=km)
      if(km==2){
        break
      }
    }

    WNS <- "run with known normal"
    preN <- norm.cell.names
     ##relative expression using pred.normal cells
  	norm.mat.relat <- norm.mat.smooth-basel

  }else {
      basa <- baseline.norm.cl(norm.mat.smooth=norm.mat.smooth, min.cells=5, n.cores=n.cores)
      basel <- basa$basel
      WNS <- basa$WNS
      preN <- basa$preN
      CL <- basa$cl
      if (WNS =="unclassified.prediction"){
        Tc <- colnames(rawmat)[which(as.numeric(apply(rawmat[which(rownames(rawmat) %in% c("PTPRC", "LYZ", "PECAM")),],2, mean)) >1)]; length(Tc)
        preN <- intersect(Tc, colnames(norm.mat.smooth))

        if(length(preN)> 5){
          print("start manual mode")
          WNS <- paste("copykat failed in locating normal cells; manual adjust performed with ", length(preN), " immune cells", sep="")
          print(WNS)
          basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% preN)], 1,mean)

            }else{
                    basa <- baseline.GMM(CNA.mat=norm.mat.smooth, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99,RE.before=basa,n.cores=n.cores)
                    basel <-basa$basel
                    WNS <- basa$WNS
                    preN <- basa$preN

            }
      }
  ##relative expression using pred.normal cells
  norm.mat.relat <- norm.mat.smooth-basel

  }

  ###use a smaller set of genes to perform segmentation
  DR2 <- apply(rawmat3,1,function(x)(sum(x>0)))/ncol(rawmat3)
  ##relative expression using pred.normal cells
  norm.mat.relat <- norm.mat.relat[which(DR2>=UP.DR),]

  ###filter cells
  anno.mat2 <- anno.mat[which(DR2>=UP.DR), ]

  ToRemov3 <- NULL
  for(i in 8:ncol(anno.mat2)){
    cell <- cbind(anno.mat2$chromosome_name, anno.mat2[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    } else if(length(rle(cell[,1])$length)<23|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i<- i+1
  }

  if(length(ToRemov3)==ncol(norm.mat.relat)) stop ("all cells are filtered")

  if(length(ToRemov3)>0){
    norm.mat.relat <-norm.mat.relat[, -which(colnames(norm.mat.relat) %in% ToRemov3)]
   # print(paste("filtered out ", length(ToRemov3), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  }

  #print(paste("final segmentation: ", nrow(norm.mat.relat), " genes; ", ncol(norm.mat.relat), " cells", sep=""))

  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]

  print("step 5: segmentation...")
  results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = KS.cut, n.cores=n.cores)

  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*KS.cut, n.cores=n.cores)
  }

  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*0.5*KS.cut, n.cores=n.cores)
  }

  if(length(results$breaks)<25) stop ("too few segments; try to decrease KS.cut; or improve data")

  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)

  write.table(RNA.copycat, paste(sample.name, "CNA_raw_results_gene_by_cell.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

  print("step 6: convert to genomic bins...") ###need multi-core
  Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat=RNA.copycat, n.cores = n.cores)

  uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])

  print("step 7: adjust baseline ...")

if(cell.line=="yes"){

  mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
  write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)

  if(distance=="euclidean"){
    hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
  }else {
    hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
  }


  saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

#plot heatmap
 print("step 8: ploting heatmap ...")
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

  chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)


  if (ncol(mat.adj)< 3000){
    h <- 10
  } else {
    h <- 15
  }

  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
#library(parallelDist)

  if(distance=="euclidean"){
  jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
   heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
            ColSideColors=chr1,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

  dev.off()
  } else {
    jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
    heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                 ColSideColors=chr1,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

    dev.off()
  }
  end_time<- Sys.time()
  print(end_time -start_time)

  reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
  names(reslts) <- c("CNAmat","hclustering")
  return(reslts)
} else {
  ################removed baseline adjustment
  if(distance=="euclidean"){
  hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =n.cores, method = distance), method = "ward.D")
  }else {
  hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
  }
  hc.umap <- cutree(hcc,2)
  names(hc.umap) <- colnames(results.com)

  cl.ID <- NULL
  for(i in 1:max(hc.umap)){
    cli <- names(hc.umap)[which(hc.umap==i)]
    pid <- length(intersect(cli, preN))/length(cli)
    cl.ID <- c(cl.ID, pid)
    i<- i+1
  }

  com.pred <- names(hc.umap)
  com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
  com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "nondiploid"
  names(com.pred) <- names(hc.umap)

  ################removed baseline adjustment
  results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)
  results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
  results.com.rat.norm <- results.com.rat[,which(com.pred=="diploid")]; dim(results.com.rat.norm)

  cf.h <- apply(results.com.rat.norm, 1, sd)
  base <- apply(results.com.rat.norm, 1, mean)

  adjN <- function(j){
    a <- results.com.rat[, j]
    a[abs(a-base) <= 0.25*cf.h] <- mean(a)
    a
  }


  mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = n.cores)
  adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
  colnames(adj.results) <- colnames(results.com.rat)

  rang <- 0.5*(max(adj.results)-min(adj.results))
  mat.adj <- adj.results/rang

  print("step 8: final prediction ...")

  if(distance=="euclidean"){
    hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
  }else {
    hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
  }

  hc.umap <- cutree(hcc,2)
  names(hc.umap) <- colnames(results.com)

  saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

  cl.ID <- NULL
  for(i in 1:max(hc.umap)){
    cli <- names(hc.umap)[which(hc.umap==i)]
    pid <- length(intersect(cli, preN))/length(cli)
    cl.ID <- c(cl.ID, pid)
    i<- i+1
  }

  com.preN <- names(hc.umap)
  com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
  com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
  names(com.preN) <- names(hc.umap)

  if(WNS=="unclassified.prediction"){
    com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
    com.preN[which(com.preN == "nondiploid")] <- "c2:aneuploid:low.conf"
  }

  print("step 9: saving results...")
  res <- cbind(names(com.preN), com.preN)
  colnames(res) <- c("cell.names", "copykat.pred")

  write.table(res, paste(sample.name, "prediction.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)

  ####save copycat CNA
  write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)


  ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
  print("step 10: ploting heatmap ...")
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

  chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)

  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]

  cells <- rbind(compreN_pred,compreN_pred)

  if (ncol(mat.adj)< 3000){
    h <- 10
  } else {
    h <- 15
  }

  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

  if(distance=="euclidean"){
  jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
   heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

  legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)
  dev.off()
  } else {
    jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
    heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                 ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

    legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)

    dev.off()
  }
  end_time<- Sys.time()
  print(end_time -start_time)

  reslts <- list(res, cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
  names(reslts) <- c("prediction", "CNAmat","hclustering")
  return(reslts)
}

}

