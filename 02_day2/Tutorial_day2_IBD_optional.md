# (Optional) Investigate structure in relation to geography <!-- omit from toc -->

>To be run on the Amazon server

## Table of contents <!-- omit from toc -->
- [1. Pairwise differentiation between populations](#1-pairwise-differentiation-between-populations)
- [2. Isolation by distance](#2-isolation-by-distance)

## 1. Pairwise differentiation between populations
We will calculate here F<sub>ST</sub> between all pairs of populations. F<sub>ST</sub> accross the genome are expected to be largely driven by neutral markers while peaks of F<sub>ST</sub> may be related to selection for local adaptation. So we will first try to get a sense of the global pattern and then look for outliers.

To get pairwise F<sub>ST</sub>, we will use the R package `StAMPP`. To save some time today, we will use a toolbox developped by Yann Dorant. You may be interested in looking at the scripts to understand how this is done.

This toolbox embeds various useful scripts in order to fastly convert a file in vcf format to common population genomics formats such as genepop, StAMPP, baypass, bayenv, among other. If you are interested to learn more about this toolbox, you will find the full description at https://gitlab.com/YDorant/Toolbox

To download the toolbox in your current working directory on the server (02_day2), use the following command line:
```bash
git clone https://gitlab.com/YDorant/Toolbox
```

To run the Toolbox, the vcf file needs to be unzipped, which can be achieved with this command:
```bash
gunzip -k populations_canada_random/populations.snps.vcf.gz
```

Ok, now we are ready to convert our VCF files to the StAMPP file format. The toolbox have an easy way to do that with a bash script. This bash script requires four arguments:
* `-v` VCF file
* `-p` population map
* `-f` output file format
* `-o` output prefix name

```bash
# make the bash script executable with chmod +x
chmod +x Toolbox/00-VCF_Reshaper.sh
# convert the file format
bash Toolbox/00-VCF_Reshaper.sh -v populations_canada_random/populations.snps.vcf -p documents/popmap_canada.txt -f StAMPP -o canada
```

Check your current directory using `ls`. You should be able to see your `*.StAMPP` input file (for a quick overview of the file use the `less -S` command).

Then, we can run the R script `StAMPP-fst.R` to compute pairwise F<sub>ST</sub> for each dataset. This script requires three arguments:
* Input StAMPP file
* Output prefix
* Number of CPU allowed (default CPU=1)

**OBS! You may see some error messages indicating that temporary files cannot be deleted, which is not important.**

>Before using the tool, you have to install it on you own server session. *If not done already*, (on the server) open an R session by typing `R` in the Terminal. Proceed to install the required packages using the following lines. Use the commands below one-by-one, and answer "yes" if asked to install additional libraries. This process may take a few minutes and print on the screen lots of text. Please be patient and don't close the Terminal:
```R
install.packages("BiocManager")
install.packages("varhandle")
# install StAMPP via BioConductor pathway
BiocManager::install("StAMPP")
# once all is installed, quit R
q()
```

Here we are ready to run StAMPP and perform pairwise F<sub>ST</sub> calculations. (On the server) type in the Terminal:
```bash
Rscript Toolbox/StAMPP-fst.R canada.StAMPP canada 1
```
Using 1 CPU, each F<sub>ST</sub> calculation should take around 5-6 minutes.

Once F<sub>ST</sub> calculations are done, you will see that four output files have been generated per dataset:
* prefix_fst_bootstrapes.txt
* prefix_fst_matrix.txt
* prefix_fst_pvalue.txt
* prefix_fst_reshape.txt

Move these files to the `FST/` directory using `mv canada_fst*.txt FST`, then go to this directory with `cd FST`.

You can explore each file with the command `less -S file.txt`. Today, we will focus on the F<sub>ST</sub> matrix.

So, now you can export the whole `02_day2` directory to your local computer.

You will have the pairwise F<sub>ST</sub> matrix (file suffix `_fst_matrix.txt`) in the subfolder `FST/`, and the info files about populations in the sub-directory `documents/`:
* documents/info_samples.csv
* documents/popmap_canada.txt


**[On your local computer]**
In Rstudio on your computer, set you working directory as `02_day2`.

We keep it simple and do a simple numeric matrix but you can imagine more fancy ways, with heatmaps or so.

First, load the required libraries:
```R
# load packages
library(dplyr)
library(magrittr)
library(tibble)
library(gplots)
library(RColorBrewer)
library(corrplot)
```
Second, use the `makeSymm()` function provided below to fill in the upper diagonal of the F<sub>ST</sub> matrix (which is a triangular matrix). This step is necessary to be able use the heatmap function for plotting. To create such function, copy/paste the commands below and press `Enter`:
```R
makeSymm <- function(m, position) {
  # add symmetrical triangle matrix (upper or lower)
  if (position == "upper") {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  if (position == "lower") {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
}
```
Now, load the data, arrange it and then plot.

```R
# load the FST matrix for all SNPs
fst.mat <- read.table("FST/canada_fst_matrix.txt")

# use the given function to fill the upper diagonal of the matrix
fst.all.mat <- fst.mat %>%
               as.matrix(.) %>%
               makeSymm(., "upper")

fst.all.mat[is.na(fst.all.mat)] <- 0 # replace NAs by 0 (NAs unaccepted for the heatmap function)
fst.all.mat[1:10, 1:10] # check the fst_matrix
              
# visualise values
corrplot(fst.all.mat, is.corr = FALSE, method = "number", addgrid.col = FALSE, diag = FALSE, type = "lower", number.digits = 3, number.cex = 0.7)

# visualize pairwise FST with a heatmap plot
gplots::heatmap.2(fst.all.mat, trace = "none",
                  col = colorRampPalette(brewer.pal(9, "Reds"))(15),
                  key.xlab = "FST")
```
Note that values go up to F<sub>ST</sub> = 0.02! But most of them are very low this value.
![img_fst_all](../images/Fst_heatmap_all_snps.png)

>What do you notice? Is it heterogeneous? Do some population look more differentiated than others?

>Why do you think A and J are so different? 

## 2. Isolation by distance
To explore whether the observed pattern relates to the geographic distance between populations, we will perform an Isolation-by-distance (IBD) test:
```R
# load packages
library(reshape2)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)

# import information about populations
info_pop <- read.table("documents/info_pop_geo_eco.txt", header = TRUE)
head(info_pop)
# calculate geographic (euclidian) distances between all pairs of populations
distance <- geosphere::distm(info_pop[, c(3, 4)], fun = distGeo) %>%
  as.matrix(.)
# change it from meters to km
distance <- distance/1000

# set the colnames and rownames of the distance matrix
dimnames(distance) <- list(info_pop$pop, info_pop$pop)
distance

# prepare datasets
# linearize the distance matrix
dist.melt <- reshape2::melt(distance) %>%
  set_colnames(., c("pop1", "pop2", "distance"))
head(dist.melt)

# linearize the fst matrix
fst.melt <- reshape2::melt(fst.all.mat) %>%
  set_colnames(., c("pop1", "pop2", "FST"))

# join the distance and fst
IBD.df <- left_join(dist.melt, fst.melt, by = c("pop1", "pop2")) %>%
  filter(., distance > 0)
head(IBD.df)

# test association with FST
cor.test(log(IBD.df$distance), IBD.df$FST / (1 - IBD.df$FST))

# plot IBD
ggplot(IBD.df) + aes(x = log(distance), y = FST / (1 - FST)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) +
  theme_bw()
```
![img_IBD](../images/IBD_plot_all_snps.png)

The results of this analysis indicate that isolation-by-distance is not significant, and it does not seem that geography can explain the genetic distances very well with this full dataset.

>If we come back to our heatmap, we can notice the cluster of populations C, F, I. How do you interpret it?

Let's look at the sex ratio in our data:
```R
info_ind <- read.table("documents/info_samples_canada.txt", header = TRUE)
head(info_ind)
table(info_ind$pop, info_ind$sex)
```
```R
    F  M
  A  0 20
  B 10 10
  C  0 20
  D 15  5
  E 10 10
  F  5 15
  G 10 10
  H 10 10
  I  0 20
  J 18  2
  K 10 10
  L 10 10
 ```
It seems that the field sampling has not been very good at balancing sex-ratio between populations. *We should be worried about sex-linked markers driving the clustering pattern!*. This may be one of the reasons why A and J are the most diferrentiated lineages.

**OBS! We intentionaly subsetted the dataset to create this bias for learning purposes, as this issue may easily happen for some species or low sample sizes (the bias was controlled for in the publication).**

In the PCA analysis of the Canadian populations, we saw that a region on chromosome 4 and sex-linked markers on chromosome 5 were driving the population structure pattern. Will that influence our pairwise FST estimates? Possibly.

Let's re-run the steps above on the vcf in which we removed the chr5 and chr4 (skip that if you are late).

Let's look at the pairwise F<sub>ST</sub> matrix and the IBD stats
![img_fst_all](../images/Fst_heatmap_no_chr4-5.png)

>What do you see now?
>Pay attention to the absolute value of F<sub>ST</sub>, to the clustering and IBD (or the absence of).
