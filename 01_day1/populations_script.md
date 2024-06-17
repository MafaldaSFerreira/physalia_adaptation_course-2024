Here is the script for the Canada dataset:
```bash
cd ~/scripts
### The command below is to count the number of different populations present in the dataset to set the -p command in the program populations
cut -f 2 popmap_canada_day1.txt | sort | uniq | wc -l
cd ~/stacks
mkdir populations_canada_random
populations -P ~/stacks/gstacks/ -M ~/scripts/00_documents/popmap_canada_day1.txt -O populations_canada_random -t 4 -p 12 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp

```

Here is the script for the complete dataset:
```bash
cd ~/scripts
cut -f 2 popmap_all_day1.txt | sort | uniq | wc -l
cd ~/stacks
mkdir populations_all_random
populations -P ~/stacks/gstacks/ -M ~/scripts/00_documents/popmap_all_day1.txt -O populations_all_random -t 4 -p 14 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp

```
