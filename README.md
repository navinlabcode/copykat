# CopyKAT: Inference of genomic copy number and subclonal structure of human tumors from high-throughput single cell RNAseq data

A major challenge for single cell RNA sequencing of human tumors is to distinguish cancer cells from non-malignant cell types, as well as the presence of multiple tumor subclones. CopyKAT(Copynumber Karyotyping of Tumors) is a computational tool using integrative Bayesian approaches to identify genome-wide aneuploidy at 5MB resolution in single cells to separate tumor cells from normal cells, and tumor subclones using high-throughput sc-RNAseq data. The underlying logic for calculating DNA copy number events from RNAseq data is that gene expression levels of many adjacent genes can provide depth information to infer genomic copy number in that region. CopyKAT estimated copy number profiles can achieve a high concordance (80%) with the actual DNA copy numbers obtained by whole genome DNA sequencing. The rationale for prediction tumor/normal cell states is that aneuploidy is common in human cancers (90%). Cells with extensive genome-wide copy number aberrations (aneuploidy) are considered as tumor cells, wherease stromal normal cells and immune cells often have 2N diploid or near-diploid copy number profiles.  In this vignette, we will go through examples of calculating single cell copy number profiles from 10X single cell RNA data, predicting tumor and normal cells, and inferring tumor subclones from using the {copykat} R package. 

## Step 1: installation
Installing copykat from GitHub
```{r, eval=FALSE}
library(devtools)
install_github("navinlabcode/copykat")
```

To update between versions, please remove old version with the following codes and then reinstall it with the above codes. 
```
remove.packages("copykat")
detach("package:copykat")
``` 

### The current version is V1.1.0. Updated on Sep23, 2022
### Changes in V1.1.0:
Fixed issues in outputing 'not.defined' cells in final prediction table.

### Changes in V1.0.8:
Introduced methods for calculating copy numbers from mouse scRNAseq data (to run mouse module, set genome="mm10" in the main function).  This version outputs single cell copy number results in gene by cell dimension.  Gene names are plotted in the bottom of heatmap.  Zooming into the heatmap to read gene names.

### Changes in V1.0.6.
Two coordinate errors related to new hg20 contigs were fixed. 
Added heatmap plot of single cell copy number results with gene by cell matrix; genenames are plotted in the heatmap in the PDF files. Zooming into the bottom of the heatmap to find gene names in each segment.  Due to the large file size, it could be slow. Default is: plot.genes="TRUE". Users can change it to plot.genes="FALSE" if gene names are not wanted.
Added the filtered cells back to the prediction results
Output the *.seg file, which can be loaded to IGV viewer to visualize the results directly.  default: output.seg="FALSE".  Users can change to:output.seg="TRUE". 

### Instruction continued
An example raw UMI matrix from a breast tumor sequenced by 10X 3'RNAseq protocol is included with this package named, exp.rawdata.

To test the package, simply issue this line of code in R/Rstudio:

> copykat.test <- copykat(rawmat=exp.rawdata, sam.name="test")


## Step 2: prepare the readcount input file
The only input file that you need to prepare to run copykat is the raw gene expression matrix, with gene ids in rows and cell names in columns. The gene ids can be gene symbol or ensemble id.  The matrix values are often the count of unique molecular identifier (UMI) from nowadays high througput single cell RNAseq data. The early generation of scRNAseq data may be summarized as TPM values or total read counts, which should also work. Below I provide an example of generating this UMI count matrix from 10X output.

### An example to generate the input file from 10X genomics cellranger V3 output
```{r, eval=FALSE}
  library(Seurat)
  raw <- Read10X(data.dir = data.path.to.cellranger.outs)
  raw <- CreateSeuratObject(counts = raw, project = "copycat.test", min.cells = 0, min.features = 0)
  exp.rawdata <- as.matrix(raw@assays$RNA@counts)
```

I can save the matrix for future usage.
```{r, eval=FALSE}
write.table(exp.rawdata, file="exp.rawdata.txt", sep="\t", quote = FALSE, row.names = TRUE)
```


In this vignette, we will use the example UMI count matrix, exp.rawdata to demonstrate the workflow.

## Step 3: running copykat
Now I have prepared the only one input, raw UMI count matrix, I am ready to run copykat. The default gene ids in cellranger output is gene symbol, so I put "Symbol" or "S". To filter out cells, I require at least 5 genes in each chromosome to calculate DNA copy numbers. I can tune this down to ngene.chr=1 to keep as many cells as possible, however I think using at least 5 genes to represent one chromosome is not very stringent. To filter out genes, I can tune parameters to keep only genes that are expressed in LOW.DR to UP.DR fractions of cells. I put default LOW.DR=0.05, UP.DR=0.2. I can tune down these values to keep more genes in the analysis. I need to make sure that LOW.DR is smaller than UP.DR though.  

I ask copykat to take at least 25 genes per segment. I can play around with other options ranging 15-150 genes per bin. KS.cut is the segmentation parameter, ranging from 0 to 1. Increasing KS.cut decreases sensitivity, i.e. less segments/breakpoints. Usually it works in a range of 0.05-0.15. 

Here I do parallel computation by setting n.cores = 4. Default value is 1. 

I also give a sample name by setting sam.name="test". 

One struggling  observation is that none of one clustering method could fit all datasets. In this version, I add a distance parameters for clustering that include "euclidean" distance and correlational distance, ie. 1-"pearson" and "spearman" similarity.  In general, corretional distances tend to favor noisy data, while euclidean distance tends to favor data with larger CN segments. 

I add an option to input known normal cell names as a vector object. Default is NULL.

I add a mode for cell line data that has only aneuploid or diploid cells. Setting this cell line mode by cell.line="yes". Default for tissue samples is cell.line="no". This cell line mode uses synthetic basline from the data variations, which does not represent the published algorithm. This cell line mode does not guarantee the success nor the accuracy.

I add an option to output seg file which can be directly loaded to IGV viewer for visualization. Default is output.seg="FALSE". Please change to output.seg= "TRUE" if seg file is wanted.

I add an additional plot of single cell copy number results, using gene by cell matrix.  This is by default, plot.genes="TRUE". Gene names are labelled at the bottom of heatmap.  Need to zoom in to read the tiny fonts.

Mouse scRNAseq is supported due to many requests.  I didn't have many mouse datasets to test the method yet. It is using the same method as human data.  Only difference is that the result is output in gene space instead of genomic space.  Meaning the locations of CNVs is labelled by gene names, instead of genomic positions.  Although aneuploid/diploid prediction was forced to execute, but users need to eyeball to check its accuracy.

Now run the code:

```{r, message=FALSE}
library(copykat)
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=1)
```

It might take a while to run a dataset with more than 10,000 single cells. It is suggested to run one sample at a time.  Combining different sample would pick up batch effects between samples.

After this step, copykat automatically save the calculated copy number matrix, the heatmap and tumor/normal prediction results in my working directory.  Users can also extract them from the object.  Two examples as follows:

```{r, eval=TRUE}
pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
CNA.test <- data.frame(copykat.test$CNAmat)
```

## Step 4: navigate prediction results
Now let's look at the prediction results. Predicted aneuploid cells are inferred as tumor cells; diploid cells are stromal normal cells.

This version put back all filtered cells and label them as "not defined" cells.

```{r, eval=TRUE}
head(pred.test)
     
```
The first 3 columns in the CNV matrix are the genomic coordinates. Rows are 220KB bins in genomic orders.
It also output CNVs indexed by gene names and other information such as chromosome names, start and end positions, G staining band etc.
Mouse result only output the CNV results in gene name vy cell space.

```{r, eval=TRUE}
head(CNA.test[ , 1:5])
```
Copykat generate a heatmap plot for estimated copy numbers. Rows are single cells; columns are 220kb bins in genomic order.

```{r, eval=TRUE, , echo=FALSE, fig.width=7, fig.height=7}
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

  chr <- as.numeric(CNA.test$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)

  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  com.preN <- pred.test$copykat.pred
  pred <- rbPal5(2)[as.numeric(factor(com.preN))]

  cells <- rbind(pred,pred)
  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

  legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")

```

## Step 5: define subpopulations of aneuploid tumor cells
I observed both diploid and aneuploid cells. Next step is to extract aneuploid cells that are considered as tumor cells in aneuploid tumors to define two copy number subpopulations of single tumor cells. 

```{r, eval=TRUE, fig.width=7, fig.height=6}
tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,2)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)

heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

  legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')

```
Now I have defined two subpopulations with major subclonal differences, I can move forward to compare their gene expression profiles, and evaluate gene dosage effects of subclonal copy number changes.


A final note, I also put some useful annotation data along with copykat: 
1) full.anno: the full annotation of 56051 genes (hg38) including absolute positions, chr, start, end, ensemble id, gene symbol and G band.
2) full.anno.mm10: the full annotation of 51450 genes (mm10) including absolute positions, chr, start, end, ensemble id, gene symbol and G band.
3) DNA.hg20: the coordinates in hg38 of the 220kb variable bins excluding Y chromosome.
4) cyclegenes: that are removed from copykat analysis.
5) exp.rawdata: UMI matrix from a breast tumor.

Final final note, CopyKAT had difficulty in predicting tumor and normal cells in the cases of pediatric and liquid tumors that have a few CNAs.  CopyKAT provides two ways to bypass this to give certain output instead of being dead staright: 1) input a vector of cell names of known normal cells from the same dataset 2) or try to search for T cells.


## Citation:
Gao, R., Bai, S., Henderson, Y. C., Lin, Y., Schalck, A., Yan, Y., Kumar, T., Hu, M., Sei, E., Davis, A., Wang, F., Shaitelman, S. F., Wang, J. R., Chen, K., Moulder, S., Lai, S. Y. & Navin, N. E.  (2021). Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. Nat Biotechnol.  doi:10.1038/s41587-020-00795-2.

