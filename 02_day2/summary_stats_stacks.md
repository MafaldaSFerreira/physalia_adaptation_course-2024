Stacks generates files summarizing several summary statistics, including F<sub>ST</sub>. Let's look at one of these files in more detail:
```bash
cd ~/stacks/populations_2lin_random
head -n 30 populations.fst_NWA-GRE.tsv
```
F<sub>ST</sub> for each locus is at the 8<sup>th</sup> column. We can look at the 20 most differentiated SNPs with this command:
```bash
grep -v "^#" populations.fst_NWA-GRE.tsv | sort -rk8,8 | head -n 20
```
The command `grep -v` excludes lines that start with a `#`, then you sort the file by values of F<sub>ST</sub> from the highest to the lowest (`-r`), and print to the screen the first 20 lines.

If you have more than two populations, Stacks will print out `populations.fst_POP1-POP2.tsv` for each unique combination of populations. Look into the `populations_canada_random` folder. 
```bash
cat populations.fst_*.tsv | grep -v "^#" - | sort -rk8,8 | head -n 20
```

>Have a look at this list. Are these SNPs randomly distributed across the 5 chromosomes?

Here we can create a list of the 100 most differentiated SNPs among all populations with this one-liner:
```bash
cat populations.fst_*.tsv | grep -v "^#" - | sort -rk8,8 | head -n 1000 | sort -un -k1,1 -k2,2 -s | sort -rk8,8 | cut -f 1,6 | tail -n +2 | head -n 100 | sort -n > high_fst.whitelist.tsv
```
Note that it is mostly based on commands that we've used above, with the difference that we merge all files containing population pairwise information and print only the SNP information. Note that even though the data are mapped to a reference genome, Stacks has built a catalog of loci and uses locus ID and position within the locus (stored in columns 1 and 6 in this file) to retrieve information on these SNPs (i.e. Stacks doesn't use chromosome and position to identify SNPs).

Now that we have created a "whitelist" listing the 100 most differentiated markers, we can rerun populations and generate F-statistics and input files for these SNPs only:
```bash
cd ~/stacks
mkdir populations_canada_random_highfst
populations -t 16 -P ~/stacks/gstacks/ -M ~/scripts/popmap_canada.txt -O populations_canada_random_highfst --fstats --vcf --genepop --structure -W ~/stacks/populations_canada_random/high_fst.whitelist.tsv
```

Whitelists and blacklists are very useful to select or discard loci for analysis. For quick data exploration, or for analyses that don't require or can't handle large datasets, you can use a similar code to create a whitelist with a random subset of SNPs. For 1000 SNPs for example, one can use:
```bash
grep -v "^#" populations.sumstats.tsv | cut -f 1,4 | sort | uniq | shuf | head -n 1000 | sort -n > 1000snps_whitelist.tsv
```
Note that Stacks takes into account the whole locus, rather than just the random SNP selected in each locus. You can examine this in the file called `populations.sumstats_summary.tsv`. For easier visualization, you can download the file and import it in excel. Additionally, Stacks provides many haplotype-based statistics. Although we won't delve into these today, these statistics are important to obtain estimates of D<sub>xy</sub> (a measure of absolute divergence, gene and haplotype diversity) and to phase SNPs, which can provide additional important information for both evolutionary questions and conservation applications.

We recommend to take a look at the following review on the use of haplotype information in conservation genomics
>Leitwein, M., Duranton, M., Rougemont, Q., Gagnaire, P.A. and Bernatchez, L., 2020. Using haplotype information for conservation genomics. Trends in Ecology & Evolution, 35(3), pp.245-258.

Also, note that statistics will differ depending on whether you look at all SNPs in one locus or only one random SNP per locus. For example, statistics from `~/stacks/populations_2lin/populations.log`, including all SNPs:
```bash
Population summary statistics (more detail in populations.sumstats_summary.tsv):
  NWA: 39.172 samples per locus; pi: 0.054733; all/variant/polymorphic sites: 770730/45238/25848; private alleles: 18011
  GRE: 39.138 samples per locus; pi: 0.054095; all/variant/polymorphic sites: 770730/45238/27188; private alleles: 19351

Population pair divergence statistics (more in populations.fst_summary.tsv and populations.phistats_summary.tsv):
  NWA-GRE: mean Fst: 0.024875; mean Phi_st: 0.077792; mean Fst': 0.074473; mean Dxy: 0.0049455
```
Statistics from `~/stacks/populations_2lin_random/populations.log`, including only one random SNP per locus:
```bash
Population summary statistics (more detail in populations.sumstats_summary.tsv):
  NWA: 39.161 samples per locus; pi: 0.051743; all/variant/polymorphic sites: 770730/8530/4634; private alleles: 3345
  GRE: 39.12 samples per locus; pi: 0.054156; all/variant/polymorphic sites: 770730/8530/5171; private alleles: 3882

Population pair divergence statistics (more in populations.fst_summary.tsv and populations.phistats_summary.tsv):
  NWA-GRE: mean Fst: 0.028809; mean Phi_st: 0.033052; mean Fst': 0.020139; mean Dxy: 0.00092926
```

Here, the value of π are not very different because they are averaged over only the polymorphic sites (which is NOT how we should estimate π anyway). Same for F<sub>ST</sub>. Do you see how different the haplotype-based estimates are? They are all lower in the random-SNP dataset because Stacks computes these statistics as if there was only one polymorphic site in each of these loci.


