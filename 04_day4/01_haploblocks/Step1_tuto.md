## Step 1. Genotyping the individuals
Thanks to our local PCA exploration, we know that there are non-recombining haploblocks which may be explained by a chromosomal inversion in chromosome 4.
We located the breakpoints approximately from 4.8MB to 16.6MB.

Now we want to classify our individuals according to the cluster that we observed on the PCA.
To do so, we will extract the relevant regions from the VCF file and then perform a PCA and clustering approach.

#### on the server : make a VCF with Chr4:4.8-16.6MB
In folder`02_data`, you will find the VCF for Canadian populations with all SNPs (1 SNP randomly selected by locus) and all individuals.
You also have a file with information about the individuals `info_samples_canada.txt`.

We will make a VCF with Chr 4 4.8-16.6MB and export it as `.012`as we did on day 2 to do the PCA.
```bash
# unzip de vcf file
gunzip 02_data/canada.vcf.gz

# extract SNPs in Chr 4 between coordinates 4800000 and 16600000
vcftools --vcf 02_data/canada.vcf --chr Chr4 --from-bp 4800000 --to-bp 16600000 --recode --out 02_data/canada.chr4inv
vcftools --vcf 02_data/canada.chr4inv.recode.vcf --012 --out 02_data/canada.chr4inv
```
Now you can copy the generated files into  on your local computer following the same architecture (`04_day4/01_haploblocks/02_data`) 

#### On your computer in R studio
Perform a PCA :
```R
# read the geno 012 file
geno <- read.table("02_data/canada.chr4inv.012")[, -1]

# read the individual information
indv <- read.table("02_data/canada.chr4inv.012.indv")

# attribute individual IDs as row names to the geno file
rownames(geno) <- indv[, 1]

# check the geno matrix
geno[1:6, 1:6] 

# run the pca
geno.pca <- prcomp(geno) 

# plot the pca
plot(geno.pca$x[, 1], geno.pca$x[, 2]) 
```
We tend to see three groups although the smallest one is not well defined, possibly because there are less individuals belonging to this group.
To cluster those individuals into three groups, we will use kmeans.

```R
# cluster individuals according to their pca loading using kmeans
geno_kmean <- kmeans(geno.pca$x[, 1], c(min(geno.pca$x[, 1]), (min(geno.pca$x[, 1]) + max(geno.pca$x[, 1])) / 2, max(geno.pca$x[, 1])))

# inspect the object
geno_kmean

# plot pca again coloring individuals according to their cluster
plot(geno.pca$x[, 1],geno.pca$x[, 2], col = geno_kmean$cluster, pch = 20)
```
![pca](06_images/pca_cluster.png)

Notice the ratio between the "between sum of squares" over the "total sum of square" in `geno_kmean` object, which is indicative of the fit of the clusters.
Here the fit is around 92%, which is not perfect but still meaningful, since we still get three clusters likely corresponding to the three inversion genotypes (AA, AB and BB).

Let's use these three clusters for subsequent analyses. In real life, you may want to be more refined. For instance, you could perform the PCA using the most strongly differentiated windows or SNPs between individuals, or excluding dubious intermediate individuals.
To split the VCF between our three groups we will need a list of individual IDs in each category.

```R
# save the individual and kmeans results into a new data.frame
info_cluster_inv <- cbind(indv, geno_kmean$cluster)
# add colnames 
colnames(info_cluster_inv) <- c("id_inv", "cluster_inv")
head(info_cluster_inv)

# prepare a list per group 
AA_ind <- info_cluster_inv[info_cluster_inv$cluster_inv == "1", 1]
AB_ind <- info_cluster_inv[info_cluster_inv$cluster_inv == "2", 1]
BB_ind <- info_cluster_inv[info_cluster_inv$cluster_inv == "3", 1]
  
# export files
write.table(info_cluster_inv, "02_data/info_cluster_inv.txt", quote = FALSE, row.names = FALSE)
write.table(AA_ind, "02_data/AA.list", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(AB_ind, "02_data/AB.list", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(BB_ind, "02_data/BB.list", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
