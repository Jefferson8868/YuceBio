# **Differential Splicing For DEXSeq, edgeR, Limma**

## **Preparing Data \(准备数据\)**
	**1.** Sam/Bam Files \(准备好Sam/Bam文件\)
	**2.** GTF Annotation File \(进行基因比对时使用的参考基因GTF文件\)

## **Preparing Count Data & Flattened GFF Annotation File \(准备外显子定量文件以及格式后的GFF文件\)
	**1.** Make sure "HTSeq" is installed in python packages. If not, using the following command:
		'pip install HTSeq'
	**2.** Using the following command:
		'python dexseq_prepare_annotation.py <input_gtf> <output_gff>'
	Inputing the GTF Annotation File, outputing the Flattened GFF Annotation File \(输入参考基因组文件，输出格式化后的参考文件\)
	**3.** Using the following command:
		'python dexseq_count.py [options] <flattened_gff_file> <alignment_file> <output_file>'
	Notice: alignment files is neither bam/sam files, output files should be in .txt format. 
	Important Argument: 
		- '--paired' \(if the file is read in paired or not\); 
		- '-s' \(Strand. Indicates whether the data is from a strand-specific assay \(default: yes, for indivisual datasets, suggest: no\)\) 
		- '-f' \(inputing sam/bam files\) 

## **Formatting Output .txt Files**
	> The output .txt files may conclude \" \(double quotation marks\) within the gene's name\). To address it, run the provided bash script.
	'sh clean_quotes.sh'
## **DEXSeq**
	> Using RStudio to run the following file \"DEXSeq.R\", for detailed information, see guidelines within \"DEXSeq.R\"  
	
## **edgeR**
	> Using RStudio to run the following file \"edgeR.R\", for detailed information, see guidelines within \"edgeR.R\"

## **Limma**
	> Using RStudio to run the following file \"limma.R\", for detailed information, see guidelines within \"limma.R\"

