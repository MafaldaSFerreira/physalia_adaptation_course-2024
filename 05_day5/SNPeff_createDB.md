#### Prepare the database for SNPeff
First download [snpEff](http://pcingola.github.io/SnpEff/download/), unzip it, and go to its respective directory. 

Then add a line in the configuration file `snpEff.config` with the name of the new genome. For this, you can use the the text editor `nano`:
```bash
# save in an environmental variable the path on the server to the executable (*.jar) file of SnpEff
PATH_TO_SNPEFF="/home/ubuntu/src/conda/envs/adaptg/share/snpeff-5.2-1"

# open the config file with nano
nano $PATH_TO_SNPEFF/snpEff.config

# scroll down using the arrow keys until the section

#-------------------------------------------------------------------------------
# Databases & Genomes
#
# One entry per genome version. 
#
# For genome version 'ZZZ' the entries look like
#       ZZZ.genome              : Real name for ZZZ (e.g. 'Human')
#       ZZZ.reference           : [Optional] Comma separated list of URL to site/s where information for building ZZZ database was ex$
#       ZZZ.chrName.codonTable  : [Optional] Define codon table used for chromosome 'chrName' (Default: 'codon.Standard')
#
#-------------------------------------------------------------------------------

#---
# Non-standard Databases
#---

```
After the `# Mouse` genome, add the lines below:
```bash
# capelin genome
genome_mallotus_dummy.genome : capelin
```
To save the changes, press `Ctrl + x`, `y`, `Enter`.

Examine the changes with `less $PATH_TO_SNPEFF/snpEff.config`.

Now, add the `*.gff` file into a data folder. and the genome as a fasta (.fa and .genome) in a genomes folder:




```bash
PATH_TO_FASTA="/home/ubuntu/Share/resources"
PATH_TO_GFF="/home/ubuntu/Share/resources"

# go to the directory where the executable file of SnpEff is located
cd $PATH_TO_SNPEFF

# create some directories
mkdir data
cd data
mkdir genome_mallotus_dummy
mkdir genomes
cd $PATH_TO_SNPEFF

# copy the GFF file
cp $PATH_TO_GFF/genome_mallotus_dummy.gff3 data/genome_mallotus_dummy/genes.gff

# copy the FASTA file
cp $PATH_TO_FASTA/genome_mallotus_dummy.fasta data/genome_mallotus_dummy/sequences.fa

```

Now you can build the database for the new genome using:
```bash
java -jar $PATH_TO_SNPEFF/snpEff.jar build -noCheckCds -noCheckProtein -gff3 -v genome_mallotus_dummy

```
If you want to, you can look at the database using:
```bash
java -jar $PATH_TO_SNPEFF/snpEff.jar dump genome_mallotus_dummy | less
```
To exit `less` simply press `q`.