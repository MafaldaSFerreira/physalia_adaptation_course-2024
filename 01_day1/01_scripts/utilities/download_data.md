interactive -A naiss2024-22-819 -n2 -t 2:00:00

~~~bash
ml load bioinfo-tools sratools/3.0.7
prefetch --option-file sralist.txt
fasterq-dump --split-files SRR11730659.sra
~~~

~~~bash
#!/bin/bash
#SBATCH -A naiss2024-22-819
#SBATCH -p core
#SBATCH -t 1
#SBATCH -t 24:00:00
#SBATCH -J fetch
#SBATCH -e fetch_%A_%a.err            # File to which STDERR will be written
#SBATCH -o fetch_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

ml load bioinfo-tools sratools

#cp */*.sra fastqs/

while read -r line line2 line3; do fasterq-dump --split-files $line.sra; done < ../../00_metadata/populations_to_target.txt 

for i in $(ls *.fastq); do mv $i ${i/.fastq/.fq}; done
for i in $(ls *.fq); do gunzip $i; done

while read -r SRA ID POP; do mv "$SRA".fq.gz "$ID".fq.gz; done < ../../00_metadata/populations_to_target.txt 
~~~

~~~bash
ml load bioinfo-tools samtools
samtools faidx genome_mallotus_dummy.fasta 
~~~

mv MID_13.fastq MID_13.fq

~~~bash
# Global variables
GENOMEFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/02_genome"
GENOME="genome_mallotus_dummy.fasta"
DATAFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/03_raw_reads/fastqs"
ALIGNEDFOLDER="/proj/snic2020-2-19/private/herring/users/mafalda/Physalia/04_aligned_files"
NCPU=2
~~~


~~~bash
#!/bin/bash
#SBATCH -A naiss2024-22-819
#SBATCH -p core
#SBATCH -t 5
#SBATCH -t 24:00:00
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
~~~

## Which locations did they use:
cat 02_day2/documents/info_samples.csv | cut -f6,7 -d";" | tr ";" "\t" | awk '{printf "%s N %s W\n", $1, $2}' | sort | uniq | head -n12 > coordinates.txt


# To run stacks:
../software/stacks-2.66/gstacks -I bamfiles -M 00_documents/dummy_pop.txt -O stacks/gstacks -t 3

../software/stacks-2.66/populations -P stacks/gstacks -M 00_documents/dummy_pop.txt -O stacks/dummy_popmap -t 4 -p 2 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp




while read -r line; do sort --random-sort $line.txt | head -n 20 | cut -f1,43 -d"," | tr "," "\t" > $line.pop.txt; done < populations_letters.txt

while read -r line; do awk -v filename=$line '{printf "%s\t%s\t%s\n", $1, $2, filename}' $line.pop.txt; done < populations_letters.txt