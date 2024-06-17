interactive -A naiss2024-22-819 -n2 -t 2:00:00

```bash
ml load bioinfo-tools sratools/3.0.7
prefetch --option-file sralist.txt
fasterq-dump --split-files SRR11730659.sra
```

```bash
#!/bin/bash
#SBATCH -A naiss2024-22-819
#SBATCH -p core
#SBATCH -t 1
#SBATCH -t 1:00:00
#SBATCH --array=0-280:1
#SBATCH -J fetch
#SBATCH -e fetch_%A_%a.err            # File to which STDERR will be written
#SBATCH -o fetch_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

ml load bioinfo-tools sratools

#cp */*.sra fastqs/

sra_files=(*.sra)

# obtain the file
file="${sra_files[$SLURM_ARRAY_TASK_ID]}"

fasterq-dump --split-files $file

mv ${file/.sra/.fastq} ${file/.sra/.fq}
do gzip ${file/.sra/.fq}

#while read -r SRA ID POP; do mv "$SRA".fq.gz "$ID".fq.gz; done < ../../00_metadata/populations_to_target.txt 
```

while read -r SRA ID POP; do cp "$SRA".fq.gz "$ID".fq.gz; done < ../../00_metadata/populations_to_target.txt 

```bash
ml load bioinfo-tools samtools
samtools faidx genome_mallotus_dummy.fasta 
```

mv MID_13.fastq MID_13.fq

```bash
# Global variables
GENOMEFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/02_genome"
GENOME="genome_mallotus_dummy.fasta"
DATAFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/03_raw_reads/fastqs"
ALIGNEDFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/04_aligned_files"
NCPU=2
```


```bash
#!/bin/bash
#SBATCH -A naiss2024-22-819
#SBATCH -p core
#SBATCH -t 5
#SBATCH -t 2:00:00
#SBATCH -J mapping
#SBATCH -e mapping_%A_%a.err            # File to which STDERR will be written
#SBATCH -o mapping_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

# Global variables
GENOMEFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/02_genome"
GENOME="genome_mallotus_dummy.fasta"
DATAFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/03_raw_reads/fastqs"
ALIGNEDFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/04_aligned_files"
NCPU=5

ml load bioinfo-tools samtools bwa

# Index genome if not alread done
#bwa index -p "$GENOMEFOLDER"/"$GENOME" "$GENOMEFOLDER"/"$GENOME"

for file in $(ls -1 "$DATAFOLDER"/*.fq.gz)
do
    # Name of uncompressed file
    echo "Aligning file $file"

    name=$(basename "$file")
    ID="@RG\tID:ind\tSM:ind\tPL:IonProton"

    # Align reads 1 step
    bwa mem -t "$NCPU" -k 19 -c 500 -T 0 -E 2,2 -R "$ID" "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" 2> /dev/null | samtools view -Sb -q 1 -F 4 -F 256 -F 2048 > "$ALIGNEDFOLDER"/"${name%.fq.gz}".bam

    # Samtools sort
    samtools sort --threads "$NCPU" -o "$ALIGNEDFOLDER"/"${name%.fq.gz}".sorted.bam "$ALIGNEDFOLDER"/"${name%.fq.gz}".bam

    # Cleanup
    rm "$ALIGNEDFOLDER"/"${name%.fq.gz}".bam
done
```

```bash
#!/bin/bash
#SBATCH -A naiss2024-22-819
#SBATCH -p core
#SBATCH -t 4
#SBATCH -t 2:00:00
#SBATCH --array=0-280:1
#SBATCH -J mapping
#SBATCH -e mapping_%A_%a.err            # File to which STDERR will be written
#SBATCH -o mapping_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

# Global variables
GENOMEFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/02_genome"
GENOME="genome_mallotus_dummy.fasta"
DATAFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/03_raw_reads/fastqs"
ALIGNEDFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/04_aligned_files"
NCPU=4

ml load bioinfo-tools samtools bwa

# Index genome if not alread done
#bwa index -p "$GENOMEFOLDER"/"$GENOME" "$GENOMEFOLDER"/"$GENOME"

fq_files=("$DATAFOLDER"/*.fq.gz)
# obtain the file
file="${fq_files[$SLURM_ARRAY_TASK_ID]}"

echo "Aligning file $file"

name=$(basename "$file")
ID="@RG\tID:ind\tSM:ind\tPL:IonProton"

# Align reads 1 step
bwa mem -t "$NCPU" -k 19 -c 500 -T 0 -E 2,2 -R "$ID" "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" 2> /dev/null | samtools view -Sb -q 1 -F 4 -F 256 -F 2048 > "$ALIGNEDFOLDER"/"${name%.fq.gz}".bam

# Samtools sort
samtools sort --threads "$NCPU" -o "$ALIGNEDFOLDER"/"${name%.fq.gz}".sorted.bam "$ALIGNEDFOLDER"/"${name%.fq.gz}".bam

# Cleanup
rm "$ALIGNEDFOLDER"/"${name%.fq.gz}".bam
```

## Which locations did they use:
```bash
cat 02_day2/documents/info_samples.csv | cut -f6,7 -d";" | tr ";" "\t" | awk '{printf "%s N %s W\n", $1, $2}' | sort | uniq | head -n12 > coordinates.txt
```

# To run stacks:
```bash
../software/stacks-2.66/gstacks -I bamfiles -M 00_documents/dummy_pop.txt -O stacks/gstacks -t 3

../software/stacks-2.66/populations -P stacks/gstacks -M 00_documents/dummy_pop.txt -O stacks/dummy_popmap -t 4 -p 2 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp
```


```bash
while read -r line; do sort --random-sort $line.txt | head -n 20 | cut -f1,43 -d"," | tr "," "\t" > $line.pop.txt; done < populations_letters.txt

while read -r line; do awk -v filename=$line '{printf "%s\t%s\t%s\n", $1, $2, filename}' $line.pop.txt; done < populations_letters.txt
```



# run stacks:

```bash
#!/bin/bash
#SBATCH -A naiss2024-22-819
#SBATCH -p core
#SBATCH -t 3
#SBATCH -t 5:00:00
#SBATCH -J stacks
#SBATCH -e stacks_%A_%a.err            # File to which STDERR will be written
#SBATCH -o stacks_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

# Global variables
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/scripts"
#ALIGNEDFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/04_aligned_files"

ml load bioinfo-tools Stacks/2.66

cd $WD
# gstacks:
gstacks -I bamfiles -M 00_metadata/popmap_all_day1.txt -O stacks/gstacks -t 3

```


run_stacks_2lin.sh

#!/bin/bash
###stacks_gstacks_2lin.sh
cd
mkdir -p stacks
cd stacks
mkdir -p gstacks_2lin
gstacks -I ~/bamfiles -M ~/scripts/popmap_2lin_day1.txt -O gstacks_2lin -t 3


```bash
#!/bin/bash
#SBATCH -A naiss2024-22-819
#SBATCH -p core
#SBATCH -t 3
#SBATCH -t 5:00:00
#SBATCH -J stacks
#SBATCH -e stacks_%A_%a.err            # File to which STDERR will be written
#SBATCH -o stacks_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

# Global variables
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/scripts"
#ALIGNEDFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/04_aligned_files"

ml load bioinfo-tools Stacks/2.66

cd $WD
mkdir stacks/gstacks_2lin
# gstacks:
gstacks -I bamfiles -M 00_metadata/popmap_2lin_day1.txt -O stacks/gstacks_2lin -t 3

mkdir stacks/populations_2lin_random

populations -P stacks/gstacks_2lin/ -M 00_metadata/popmap_2lin_day1.txt -O stacks/populations_2lin_random -t 4 -p 2 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp
```