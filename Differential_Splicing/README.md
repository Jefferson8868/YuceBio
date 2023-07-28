# Supplement Materials & Guiding For Differential Splicing

## MAJIQ
### Tutorial URL:
	[Reading Material](https://majiq.biociphers.org/docs/TUTORIAL%20-%20V1.1.pdf)
	(Quick Start -- Page 7)
### Required Files
	**1.** GFF3 annotation file
	**2.** Study configuration file
		<sup>
		[info]
		readlen=76
		samdir=/data/MGP/ERP000591/bam
		genome= mm10
		genome_path=/data/WASP_DATA/Genomes/goldenPath/mm10
		type=strand-specific
		[experiments]
		Hippocampus=Hippocampus1,Hippocampus2
		Liver=Liver1,Liver2
		</sup>
 
## edgeR
### Tutorial URL:
	[Reading Material](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
	(Differential Splicing -- Page 79-87)
### Required Files
	**1.** Exon Count Matrix

## DEXSeq
### Tutorial URL:
	[Reading Material](https://bioconductor.org/packages/release/bioc/manuals/DEXSeq/man/DEXSeq.pdf)
	(Additional Guiding URL For Differential Splicing)
		**1.** [From Bam/Sam Files to Visualization](https://bioconductor.riken.jp/packages/3.10/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#:~:text=Once%20you%20have%20installed%20HTSeq%2C%20you%20can%20use,the%20installation%20directory%20with%3A%20pythonScriptsDir%20%3D%20system.file%28%22python_scripts%22%2C%20package%3D%22DEXSeq%22%29)
		**2.** [From Fastq To Visualization](https://zhuanlan.zhihu.com/p/420687595)
### Required Files
	**1.** Fastq files

## limma
### Tutorial URL:
	[Reading Material](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf)
	(Differential Splicing -- Page 134-145)
### Required Files
	**1.** Grouping Information
	**2.** Reading Counts (Fastq/Counts_Result)

## JuntionSeq
### Tutorial URL:
	[Reading Material](http://hartleys.github.io/JunctionSeq/doc/JunctionSeq.pdf)
	(Require Preprocessinb Using QoRts)
	[QoRts](https://github.com/hartleys/QoRTs)
### Required Files
	**1.** Flattened GTF File
	**2.** Counts Results
	**3.** Grouping Information

## rMATS
### Tutorial URL:
	[Reading Material](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md)
	(Additional Guiding URL For Differential Splicing)
		**1.** [rMATS Introduction & Running](https://zhuanlan.zhihu.com/p/500845376)
		**2.** [Including rMATS Plotting](https://cloud.tencent.com/developer/article/1966296)
		**3.** [rMATS Installation Tutorial](https://blog.csdn.net/yaya_bioinfo/article/details/129047948)
### Required Files
	**1.** Fastq/Bam Files
	**2.** Able to Run rmats.py
	




	







