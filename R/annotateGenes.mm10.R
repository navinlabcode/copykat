#' annotate genes with reference to hg38.
#'
#' @param mat data matrix; genes in rows; cell names in columns.
#' @param ID.type gene id type: Symbol or Ensemble.
#' @param full.anno annotation file for all known genes, automatically loaded in copycat.
#'
#' @return annotations of each genes in rows with chrom and positions.
#'
#' @examples
#' test.anno.mat <- annotateGenes.hg20(mat=matx, ID.type="ENSEMBLE_id", full.anno = full.anno)
#' @export


annotateGenes.mm10 <- function(mat, ID.type="S", full.anno=full.anno.mm10){
  print("start annotation ...")

  if(substring(ID.type,1,1) %in% c("E", "e")){
    shar <- intersect(rownames(mat), full.anno$ensembl_gene_id)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$ensembl_gene_id) %in% shar),]
    anno <- anno[!duplicated(anno$mgi_symbol),]
    anno <- anno[order(match(anno$ensembl_gene_id, rownames(mat))),]
    data <- cbind(anno, mat)

  }else if(substring(ID.type,1,1) %in% c("S", "s")) {

    shar <- intersect(rownames(mat), full.anno$mgi_symbol)
    rownames(mat)[1:10]
    full.anno$mgi_symbol[]

    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$mgi_symbol) %in% shar),]
    anno <- anno[!duplicated(anno$mgi_symbol),]
    anno <- anno[order(match(anno$mgi_symbol, rownames(mat))),]
    data <- cbind(anno, mat)
  }
}

