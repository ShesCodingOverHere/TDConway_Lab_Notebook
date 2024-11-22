# Assembly to Annotation
### Date last edited: 10.22.2024

This is a step-by-step pipeline to assembling and annotation genome assemblies.

## The strains used for PacBio hifi sequencing ~100 individuals each from male and female
*D. affinis* - 141.02 (dark)

*D. azteca* - "*azteca*"

*D. algonquin* - NH2

*D. narragansett* - GL120

*D. bifasciata* male (we sequenced) -

*D. bifasciata* female (Ryan Bracewell) - 14012-0181.02

## Step 1: Assemble long reads via hifiasm on the Cluster

```
#!/bin/bash
#SBATCH --job-name=hifiasm  # Job name
#SBATCH --partition=kucg      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL,BEGIN     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tconway@ku.edu   # Where to send mail
#SBATCH --ntasks=1          # Run on a single CPU
#SBATCH --mem=240gb           # Job memory request
#SBATCH --time=7-24:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=hifiasm_M%j.log  # Standard output and error log


module load conda
conda activate hifiasm

hifiasm -o /output/path.asm -t 32 /path/to/longreads.fastq.gz
```

## Step 2: Blast to remove contaminants locally
```
blastn -task megablast -query foo.fa -remote -db nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -out foo.megablast
```
- When this is finished, go through the file and remove bacteria/virus using NCBI

## Step 3: Stop and check the BUSCO score and quast on the linux

```
# busco
nohup busco -i /hdd/Taylor/data/foo.fa -o /hdd/Taylor/data/BUSCO.foo -l /hdd/Taylor/data/diptera_odb10 -m genome --auto-lineage-euk -f &

# quast
python3 ../software/quast-5.2.0/quast.py /hdd/Taylor/data/foo.fa
```
# Stats from unmasked genomes sans bacteria

| Species/Sex | length (Mb) | # of scaffolds | N50 | BUSCO (c) | BUSCO (s) | BUSCO (d) |
|:-----------:|:------:|:--------------:|:---:|:---------:|:---------:|:---------:|
| *D. affinis* male | 209Mb | 54 | 16.8Mb |  |  |  |
| *D. affinis* female | 201Mb | 57 | 25.9Mb |  |  |  |
| *D. algonquin* male | 235Mb | 42 | 10.9Mb |  |  |  |
| *D. algonquin* female | 196Mb | 34 | 44.2Mb |  |  |  |
| *D. azteca* male | 221Mb | 38 | 33.9Mb |  |  |  |
| *D. azteca* female | 192Mb | 30 | 34.4Mb |  |  |  |


## Step 4: Align fasta to reference via mummer to rename scaffolds locally

```
nucmer --maxgap=500 -mincluster=100 reference.fasta query.fasta
delta-filter -q -r out.delta > foo.filter
show-coords -B foo.filter > foo.coords
```

- Using the Rscript "chrom_mapping.R", check PID for each alignment and rename accordingly using:

```
sed 's/scaffold_to_be_renamed/rename_it_here/g' foo.fa >temp1
sed 's/scaffold_to_be_renamed/rename_it_here/g' temp1 >temp2
sed 's/scaffold_to_be_renamed/rename_it_here/g' temp2 >temp1

and so on...
```

## Step 5: Locate and mask repeats with repeatmodelor and repeatmasker on the cluster

```
#!/bin/bash
#SBATCH --job-name=RMaffM  # Job name
#SBATCH --partition=kucg      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL,BEGIN     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tconway@ku.edu   # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=4-00:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=RMaffM_%j.log  # Standard output and error log

module load repeatmodeler
module load repeatmasker/4.0.9

#usage: sbatch RepeatMasker.args.job <fasta> <prefix>

cd $SCRATCH
mkdir RMaffinisM_pilon2

echo "STARTING"
cd RMaffinisM_pilon2
cp $HOME/$1 .

BuildDatabase -name $2 -engine ncbi $1

RepeatModeler -engine ncbi -pa 8 -database $2

RepeatMasker -pa 8 -gff -lib $2-families.fa -dir MaskerOutput$2 $1

echo done
```

- With this data, you can look at Y-linked repeat families.

## Step 6: Annotate with helixer

- Go to https://www.plabipd.de/helixer_main.html
- Input fasta
- Change "Select Lineage-specific mode" to invertebrate
- Enter GFF label name and email address
- Submit job and wait
- grep gene foo.gff > genes.txt
- Import genes.txt into spreadsheet
- Convert gff to fasta using gffread (see below for code)
- blastx Y_transcripts.fa
- look up each gene on flybase and fill out spreadsheet

```
# gffread
gffread your_transcripts.gff -g genomic_reference.fasta -w your_transcripts.fasta​
```

## Step 7: Confirm Y genes with expression data in IGV

- Make bams from RNA seq reads using hisat2 (including mismatches) and bowtie (excluding mismatches) on the clusty (see below for code)
- Import genome fasta, gff file, and the RNAseq bams into IGV and look through Y-linked genes to determined what genes are expressed, or not, in each sex/tissue

```
#!/bin/bash
#SBATCH --job-name=makin_bams  # Job name
#SBATCH --partition=kucg      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL,BEGIN     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tconway@ku.edu   # Where to send mail
#SBATCH --ntasks=1          # Run on a single CPU
#SBATCH --mem=64gb           # Job memory request
#SBATCH --time=24:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=bams_%j.log  # Standard output and error log

module load bwa
module load samtools
module load hisat2
module load bowtie2

### Using hisat2 with allowance for mismatches
# usage: sbatch coverage.job <fasta> <fastq> <output>
#Change to ouput directory
mkdir scratch/hisat_$3
cd scratch/hisat_$3

# Build index
hisat2-build $1 reference_genome_index

# Align reads
hisat2 -p 8 --dta -x reference_genome_index -U $2 --mp 0,0 -S $3.sam

# Convert to BAM
samtools view -Sb $3.sam > $3.bam

# Sort BAM
samtools sort $3.bam -o sorted_$3.bam

# Index BAM
samtools index sorted_$3.bam


### Bowtie. No mismatches
# usage: sbatch coverage.job <fasta> <fastq> <output>

cd scratch
mkdir $3_bowtie
cd $3_bowtie

# Reference index
bowtie2-build $1 reference_index

# align
bowtie2 -x reference_index -U $2 --no-unal -N 0 -p 4 -S $3.sam

# sam to bam
samtools view -Sb $3.sam > $3.bam

# sort bam
samtools sort $3.bam -o sorted_$3.bam

# index bam
samtools index sorted_$3.bam
```