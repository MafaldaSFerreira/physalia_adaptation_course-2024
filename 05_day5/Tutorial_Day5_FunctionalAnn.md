# Day 5. Functional annotation of candidate loci for adaptation <!-- omit from toc -->

Today, we will explore SNP annotation to see if they fall on genes of known function, and evaluate the distribution of outlier SNPs across exons, introns, regulatory regions, etc. We will also look at the annotation of repeats and transposable elements and compare that to our CNV detection.

## Table of contents <!-- omit from toc -->
- [Getting started](#getting-started)
- [5-1. SNP annotation with SnpEff](#5-1-snp-annotation-with-snpeff)
  - [Annotate the VCF file](#annotate-the-vcf-file)
- [5-2. Find the intersection between SNPs and genes with bedtools](#5-2-find-the-intersection-between-snps-and-genes-with-bedtools)
  - [In R (in a R terminal on the server)](#in-r-in-a-r-terminal-on-the-server)
  - [On the server (quit R with the command q() or quit() )](#on-the-server-quit-r-with-the-command-q-or-quit-)
  - [Have a look at your results](#have-a-look-at-your-results)
- [5-3. Gene ontology enrichment analysis with goseq](#5-3-gene-ontology-enrichment-analysis-with-goseq)
  - [Prepare a file with the correspondence between transcripts and gene ontology terms](#prepare-a-file-with-the-correspondence-between-transcripts-and-gene-ontology-terms)
  - [Enrichment analysis](#enrichment-analysis)


## Getting started 
First of all, please copy folder `05_day5` into your directory and open it. We will work from there all day.

```bash
cd
cp -r ~/Share/physalia_adaptation_course/05_day5 .
cd 05_day5
```

Now you should have everything we need for today. You can explore the different folders and files
`ls 02_data` in which you will find the vcf, `ls 03_outliers` in which I put the `.txt` files that we exported on **Day 3** when analyzing the association with temperature, and on **Day 4** when analyzing divergence between haploblocks or between sexes.

You will also find the list of all SNPs present in the VCF. You can look at each file with the command `head 03_outliers/SNP_pos.txt` for instance.

The annotated transcriptome (which is generally in `.gff` format) is located at the following path `~/Share/resources/`. A copy of the `.gff` transcriptome is inside the `02_data` folder.
You don't need to copy them as we have prepared simplified files (in R, simply selecting the relevant column) when needed. You may want to have a look to get a sense of what it looks like using `less ~/Share/resources/genome_mallotus_dummy.gff3`. Press `q` to exit the less visualization. 

## 5-1. SNP annotation with SnpEff 

SnpEff is a program that uses the `.gff` file and the position of each SNP to annotate the VCF.
If you work on a model species which already has a database, you are lucky! If not, you need to build a SnpEff database. 

Here, we have already built the database for you as this is a long process and takes some space on the server. 

>If you want to replicate this process, you can see the pipeline in [Make SNPeff DB](SNPeff_createDB.md), but please do not try to run it on the server!

If you want to re-create it, the `.gff` is inside the `02_data` while the reference genome is in the folder `01_day1/02_genome`.

If you want to, you can look at the database by doing:

```bash
java -jar ~/Share/resources/snpEff/snpEff.jar dump genome_mallotus_dummy | less
```
It may take a minute to open. To exit less, press `q`.

### Annotate the VCF file 
Now we can annotate the VCF file. We use a raw VCF file in the folder `02_data` and will write the output into the folder `04_snpEff` in which we will have all subsequent files related to the SnpEff analyses.
```bash
java -Xmx4g -jar ~/Share/resources/snpEff/snpEff.jar genome_mallotus_dummy 02_data/canada.vcf > 04_snpEff/canada_annotated.vcf
```
Let's look at the new VCF `less -S 04_snpEff/canada_annotated.vcf`. As you can see, it keeps the VCF format with its specific header, and adds annotation information for each SNP. However, this information is not easy to import into R as it is now.
To continue the analyses, we will use a few bash commands and awk to split the information by the symbol `|` and produce a new tab delimited table.

```bash
# convert the table delimiter from | to \t
cat 04_snpEff/canada_annotated.vcf | awk -F"|" '$1=$1' OFS="\t" | cut -f 1-9 > 04_snpEff/SNP_annotated_formatted.txt

# look at the first 25 lines
head -n 25 04_snpEff/SNP_annotated_formatted.txt

# look at the last lines
tail 04_snpEff/SNP_annotated_formatted.txt
```
This is what it looks like:

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
Chr1    53559   49:9:-  C       G       .       PASS    ANN=G   upstream_gene_variant
Chr1    94208   95:21:+ A       G       .       PASS    ANN=G   intergenic_region
Chr1    308478  248:57:+        T       G       .       PASS    ANN=G   downstream_gene_variant
Chr1    510235  370:36:+        G       A       .       PASS    ANN=A   intergenic_region
Chr1    586674  438:51:-        T       A       .       PASS    ANN=A   splice_region_variant&intron_variant!
```

Now we have much more information about each SNP. What can we do with that? For example, we could use this information to choose only intergenic SNPs if we want a putatively-neutral subset. Or we could ask how many SNPs are located in different annotation categories, such as regulatory regions, exons, etc.

Here is a small tutorial (optional) to look at SNP repartition:
[SNP repartition](tutorial_SNP_repartition_optional.md)


## 5-2. Find the intersection between SNPs and genes with bedtools
Bedtools is a great program to find the intersection between two files. It usually works on a specific [`.bed`](https://ensembl.org/info/website/upload/bed.html) format which has at least three columns (Chromosome, StartPosition, EndPosition) and 12 columns to the maximum.
It won't like having header. We will keeep a 4th column with SNP ID.

Inspected the previously prepared annotation file in `.bed` format:
```bash
less -S 05_bed/genome_mallotus_dummy_annotation_simplified.bed
```
(`q` to exit )

If you remember our outliers files, they were not formatted as `.bed` files. Because our SNP data does not cover all the genome, we want to know which genes are within a window of X kb around the SNP position. The size of this window should ideally be adjusted to the extent of LD decay in your organism (which can be assessed by plotting LD against physical distance in the genome with the output of plink that we generated yesterday, for example). For today, we will choose 10 kb but if you are curious you can explore different sizes.
To prepare the files, we will use R. 

### In R (in a R terminal on the server)

```R
# load file
outlier_temp_rda <- read.table("03_outliers/outlier_temp_rda.txt", header = TRUE)

# have a quick look
head(outlier_temp_rda)

# what's its dimension?
dim(outlier_temp_rda)

# which size around the SNP
window <- 10000

# add a vector with start position
outlier_temp_rda$start <- outlier_temp_rda$position - (window / 2)

# start position can't be negative! replace negative by 0
outlier_temp_rda$start[outlier_temp_rda$start < 0] <- 0 

# add a vector with stop position
outlier_temp_rda$stop <- outlier_temp_rda$position + (window / 2)

# have a look
head(outlier_temp_rda)

# which columns shoud we keep?
outlier_temp_rda_bed <- outlier_temp_rda[, c(2, 5, 6, 1)]

# save your file
write.table(outlier_temp_rda_bed, "05_bed/outlier_temp_rda.bed", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
```

Apart from generating a `.bed` file containing the position of genes around RDA outliers, we also need a `.bed` file for all SNPs in our dataset. In the gene enrichment analysis, we will ask if we have an enrichment of a specific category of genes in the RDA SNP outliers against this set of background genes in out entire SNP dataset:
```R
# all SNPs positions in our dataset
all_snps <- read.table("03_outliers/SNP_pos.txt", header = TRUE)

# create a start and end column
all_snps$start <- all_snps$position - (window / 2)
all_snps$start[all_snps$start < 0] <- 0 
all_snps$stop <- all_snps$position + (window / 2)
head(all_snps)

# output the bed file
all_snps_bed <- all_snps[, c(2, 4, 5, 1)]
write.table(all_snps_bed, "05_bed/all_snps.bed", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
```

If you have time, you can do the same for the outliers from BayPass, or the outliers of divergence on chromosome 4 or chromosome 5.

### On the server (quit R with the command q() or quit() )
If you made the `.bed` files in Rstudio, please copy back your formatted outliers and SNP `.bed` files into the `05_day5/05_bed` folder on the server. 
We will now run Bedtools. The command `intersect` will look for overlap between the file given with `-a`, and the file given with `-b`, the argument `-wb` will print the information coming from the annotation file (gene names, gene ontology, uniprot ID, etc). 
The command `>` redirects the output into the file of your choice.

With the command `wc -l` we will count the lines to see how many transcripts intersect with the outlier SNPs.

```bash
# find the intersection between the genomic location of outlier regions and the location of 
transcripts
bedtools intersect -a 05_bed/outlier_temp_rda.bed -b 05_bed/genome_mallotus_dummy_annotation_simplified.bed -wb > 05_bed/outlier_temp_rda.intersect

cat 05_bed/outlier_temp_rda.intersect | wc -l
head 05_bed/outlier_temp_rda.intersect
less -S 05_bed/outlier_temp_rda.intersect
```
For gene enrichment analyses, we will need a simpler format, so we will cut the 8th column (transcript name) with the command cut and keep it for goseq.
```bash
cat 05_bed/outlier_temp_rda.intersect | cut -f 8 > 06_go/outlier_temp_rda.transcript

head 06_go/outlier_temp_rda.transcript
```

We repeat the analyses for the set of all SNPs:

```bash
bedtools intersect -a 05_bed/all_snps.bed -b 05_bed/genome_mallotus_dummy_annotation_simplified.bed -wb > 05_bed/all_snps.intersect

cat 05_bed/all_snps.intersect | wc -l
cat 05_bed/all_snps.intersect | cut -f 8 > 06_go/all_snps.transcript
head 06_go/all_snps.transcript
```

You may now repeat the commands for other outliers files if you wish.

### Have a look at your results
Copying your `.intersect` files (and `.transcript`) back to your computer to `05_day5/05_bed`. Open them with a text editor, such as notepad ++ or equivalent to inspect them. You can see the name of the genes and their gene ontology codes. What do you think?

For information about gene ontology codes, you may want to check this [website](http://geneontology.org/docs/ontology-documentation/).

For information about proteins, you can look at https://www.uniprot.org/

## 5-3. Gene ontology enrichment analysis with goseq

**OBS! Gene enrichment analysis is mostly designed for RNAseq or whole-genome analyses in which all genes are analyzed. The next steps use our RAD-seq data just for instructional purposes and such that you can replicate the analysis with your whole-genome data, if you have it.**

Note: This script was built with the help or Dr. E. Berdan (U. Stockholm).

This part will only be in Rstudio on your computer.
Please copy *the whole folder `05_day5`* to your local computer, set your working directory in Rstudio to `05_day5`, and load the required package:
```R
# load package
library(goseq)
```

### Prepare a file with the correspondence between transcripts and gene ontology terms
Again, as we are using a non-model organism (and a dummy genome), there is no database of our transcriptome and GO terms in the library, so we will build it.

We have prepared a simplified annotation of the transcripts, and removed duplicates as much as possible. In fact, when a genome gets annotated, or when contigs from a transcriptome are aligned to a genome, one may have mutliple matches, which in turn may inflate our calculations. This also occurs when one SNP overlaps with a gene which has several transcripts - the SNP will have as many annotations as the number of transcripts. This needs to be taken into account by removing multiple matches:
```R
transcript_info = read.table("06_go/transcript_go_simplified.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(transcript_info)
row.names(transcript_info) <- transcript_info$TranscriptName
```
As you see, GO terms are side-by-side, which will not be super helpful to make a table with transcripts <-> GO.

Load the library `splitstackshape`, which is able to split a whole column, and the library `data.table`, which is super helpful to transform wide table into a long table:
```R
# load packages
library(splitstackshape)

# split the GO terms
go_split = cSplit(transcript_info, "GeneGo", sep = ";", type.convert = FALSE)

go_split$contig <- transcript_info$TranscriptName
head(go_split[, 1:10])
```

```R
# install packages
library(data.table)
library(dplyr)

# linearize the matrix
terms = colnames(select(go_split, contains("GeneGo")))

# transform table into a long format
go_long = melt(go_split, measure.vars = terms, id.vars = "contig", na.rm = TRUE)
head(go_long)

# convert to data frame.
go_ready = as.data.frame(go_long[, c(1, 3)])
head(go_ready)
```

### Enrichment analysis
First, we need to know which genes are nearby a SNP covered by a RAD locus. This will be a subset of all the genes present in the transcriptome/genome since we did a reduced-representation sequencing. This may also happen for a whole genome if, because of coverage or filtering we have genes not covered by SNPs.

```R
# upload transcript intersecting with snps
all_transcripts <- read.table("06_go/all_snps.transcript", header = FALSE)
colnames(all_transcripts)[1] <- "TranscriptName"
dim(all_transcripts) # how many?
head(all_transcripts)
```

Now we will use the information about the size of the gene to correct for the fact that SNPs are more likely to be located inside a long gene than a short gene.
We will use the `left_join()` command from dplyr to match two data.frames by transcript name, maintaining the rest of the information in the data.frames.

```R
library(dplyr)
# add size info
all_transcripts <- left_join(all_transcripts, transcript_info[, c(1,3)])
dim(all_transcripts)
head(all_transcripts)
```

We will then remove duplicated rows, as some genes are represented more than one time in our set of genes. This is because some genes may be duplicated in the genome, or several SNPs match to the same gene. We will use `distinct()` from `dplyr` to remove duplicates.

```R
# make unique
all_transcripts_unique <- all_transcripts %>% distinct(TranscriptName, .keep_all = TRUE)

# see the matrix has reduced
dim(all_transcripts_unique) 

# use transcript names as row names for goseq
row.names(all_transcripts_unique) <- all_transcripts_unique$TranscriptName

# check the data.frame
head(all_transcripts_unique)
```

Let's open one of our outlier list and format it. We will start with the outliers from the associations with temperature found by the RDA.

```R
# transcripts in outliers
outliers_transcripts_temp_rda <- read.table("06_go/outlier_temp_rda.transcript", header = FALSE)
colnames(outliers_transcripts_temp_rda) <- "TranscriptName"
head(outliers_transcripts_temp_rda)
dim(outliers_transcripts_temp_rda)
```
We will now add a column to the matrix with all the genes indicating whether this gene is found or not inside the outlier list. This will be a 0/1 vector.
```R
all_transcripts_unique$outliers_temp_rda <- as.numeric(all_transcripts_unique$TranscriptName %in% outliers_transcripts_temp_rda$TranscriptName)

head(all_transcripts_unique)
```

Now we run the goseq function `nullp()` to prepare the data and integrate gene length bias. This function requires vectors as input, so we will convert the data in data.frames to vectors.

```R
# all genes:
measured_genes = as.vector(all_transcripts_unique$TranscriptName)

# outlier genes:
outliers_genes = as.vector(all_transcripts_unique$outliers_temp_rda)

# gene length
length = as.vector(all_transcripts_unique$length)

# run nullp()
pwf_outliers = nullp(outliers_genes, bias.data = length)
row.names(pwf_outliers) <- row.names(all_transcripts_unique)

#check output
head(pwf_outliers)
```
Great! We have formatted all the files. Let's now test enrichment using our database (transcript/GO) named "go_ready" and the prepared list of genes with 0/1 info for our outliers from temperature association with RDA:
```R
# check for enrichment with goseq.
enrich_outliers = goseq(pwf_outliers, gene2cat = go_ready, use_genes_without_cat = TRUE)

# check output
head(enrich_outliers)
```
We will now correct for multiple testing with Benjamini & Hochberg correction.
```R
# Benjamini & Hochberg multiple test correction
enrich_outliers$over_represented_padjust <- p.adjust(enrich_outliers$over_represented_pvalue, method = "BH")

# check output
head(enrich_outliers)

# write output
write.table(enrich_outliers, "06_go/GO_enrich_temp_RDA.txt", sep = "\t")

```
All GO terms are presented in this matrix. We have exported it so you can look at it in an text editor.

We may want to see only the significant enrichments:
```R
enrich_outliers[which(enrich_outliers$over_represented_padjust < 0.10), ]
```

What do you think? Sometimes it may be significant but when one has low numbers (1 gene out of 4), is this really interpretable? 
Anyhow, you know how to do it! :) 

If you wish, you can repeat the analysis on the outliers from BayPass, or joined BP/rda, or on the outliers of divergence found on chr4 and chr5. 
