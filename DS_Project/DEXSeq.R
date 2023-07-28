
# library(refGenome)
setwd("/Users/jeffersonchen/programming/YuceBio/DS_Project")

# flattenedAnnotation <- flattenGTF(
#   # basic input/output options
#   "/Users/jeffersonchen/programming/YuceBio/DEXSeq/gencode.v38.annotation.gtf",
#   GTF.featureType = "exon",
#   GTF.attrType = "gene_id",
#   # the option specifying the merging algorithm
#   method = "merge")
BiocManager::install()
library(GenomicFeatures)
library(Rsubread)
library(pasillaBamSubset)
library(Rsamtools)
library(GenomicAlignments)

txdb = makeTxDbFromGFF("./Drosophila_melanogaster.BDGP5.25.62.gtf.gz")
gtf_files <- makeTxDbFromGFF("/Users/jeffersonchen/programming/YuceBio/DEXSeq/gencode.v38.for_dexseq.gff")
flattenedAnnotation = exonicParts(txdb, linked.to.single.gene.only=TRUE )

# bamFiles = c("/Users/jeffersonchen/programming/YuceBio/Bam_Files/Bam_Files")
# bamFiles =list.files(bamFiles, pattern = ".bam", full.names = T)
# class(bamFiles)
# bamFiles = BamFileList( bamFiles )


seqlevelsStyle(flattenedAnnotation) = "UCSC"
se = summarizeOverlaps(
  flattenedAnnotation, BamFileList(bamFiles), singleEnd=FALSE,
  fragments=TRUE, ignore.strand=TRUE )

se

colData(se)
colData(se)$condition = factor(c("ND", "D"))
colData(se)$libType = factor(c("paired-end"))
dxd = DEXSeqDataSetFromSE( se, design= ~ sample + exon + condition:exon )

dxd

colData(dxd)
head( counts(dxd), 60)
split( seq_len(ncol(dxd)), colData(dxd)$exon )
head( rowRanges(dxd), 3 )

dxd = estimateSizeFactors( dxd )

dxd = estimateDispersions( dxd )
plotDispEsts( dxd )

dxd = testForDEU( dxd )

dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

dxr1 = DEXSeqResults( dxd )
plotMA( dxr1, cex=0.8 )

formulaFullModel    =  ~ sample + exon + type:exon + condition:exon
formulaReducedModel =  ~ sample + exon + type:exon 

dxd = estimateDispersions( dxd, formula = formulaFullModel )
dxd = testForDEU( dxd, 
                  reducedModel = formulaReducedModel, 
                  fullModel = formulaFullModel )
dxr2 = DEXSeqResults( dxd )
plotDEXSeq( dxr2, "FBgn0010909", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

wh = (dxr2$groupID=="FBgn0010909")
stopifnot(sum(dxr2$padj[wh] < formals(plotDEXSeq)$FDR)==1)

plotDEXSeq( dxr2, "FBgn0010909", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr2, "FBgn0010909", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr2, "FBgn0010909", expression=FALSE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

DEXSeqHTML( dxr2, FDR=0.1, color=c("#FF000080", "#0000FF80") )


