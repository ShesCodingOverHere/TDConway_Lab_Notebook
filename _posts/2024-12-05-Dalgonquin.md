# *D. algonquin* Strain: NH2


## Assembly
```
#!/bin/bash
#SBATCH --job-name=hifiasm  # Job name
#SBATCH --partition=kucg      # Partition Name (Required)
#SBATCH --mail-type=END,FAIL,BEGIN     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tconway@ku.edu   # Where to send mail
#SBATCH --ntasks=1          # Run on a single CPU
#SBATCH --mem=100gb           # Job memory request
#SBATCH --time=7-24:00:00        # Time limit days-hrs:min:sec
#SBATCH --output=hifiasm_M%j.log  # Standard output and error log


module load conda
conda activate hifiasm

hifiasm -o /home/t043c581/scratch/data/Dalg_downsampled.asm -t 32 /home/t043c581/scratch/data/DalgM_hifi.downsampled.fastq
```

![alt text](image-16.png)

## Remove bacterial contamination

```
blastn -task megablast -query foo.fa -remote -db nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -out foo.megablast
```

- When this is finished, go through the file and remove bacteria/virus using NCBI

### quast, busco, and other data I will undoubtedly be asked for...

```
# busco
nohup busco -i /hdd/Taylor/data/foo.fa -o /hdd/Taylor/data/BUSCO.foo -l /hdd/Taylor/data/diptera_odb10 -m genome --auto-lineage-euk -f &

# quast
python3 ../software/quast-5.2.0/quast.py /hdd/Taylor/data/foo.fa
```

| Species (sex)| Length| # of Scaffolds| N50 | BUSCO (complete)| BUSCO (single)| BUSCO (dup)|
|-------------:|:-----:|:-------------:|:---:|:---------------:|:-------------:|:----------:|
| _D. algonquin_(male)| ~235Mb| 42| 10.9Mb| 98.6%| 91.1%| 7.5%|
| _D. algonquin_(female)| ~196Mb| 34| 44.1Mb| 99.3%| 98.4%| 0.9%|


## Align fasta to reference via mummer to rename scaffolds locally

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

## Coverage
### Method 1

Make .bam files with minimap and plot coverage using chromosome quotient method

#### Code:
```
# on the linux
# align female longreads to male reference
minimap2 -t 8 -ax map-hifi /hdd/Taylor/data/DalgM_hifi_nobact_v2.fa /hdd/Taylor/data/D-algonquin_F_HiFi.fastq.gz |samtools view -bS > /hdd/Taylor/data/DalgF_hifi.bam

# align male longreads to male reference
minimap2 -t 8 -ax map-hifi /hdd/Taylor/data/DalgM_hifi_nobact_v2.fa /hdd/Taylor/data/D-algonquin_M_HiFi.fastq.gz |samtools view -bS > /hdd/Taylor/data/DaztM_hifi.bam

# sort bams and get coverage
samtools sort DalgM_hifi.bam >DalgM_hifi_sorted.bam
samtools coverage DalgM_hifi_sorted.bam >DalgM_hifi.cov
samtools sort DalgF_hifi.bam >DalgF_hifi_sorted.bam
samtools coverage DalgF_hifi_sorted.bam >DalgF_hifi.cov
```

Using RStudio:
```
# Load required libraries
library(ggplot2)
library(cowplot)
library(ggrepel)

# Load data
aztmale <- read.delim("/Users/conway/Desktop/CurrentWorkingDatasets/DalgM_hifi.cov", header = TRUE)
aztfem <- read.delim("/Users/conway/Desktop/CurrentWorkingDatasets/DalgF_hifi.cov", header = TRUE)

# Rename columns for clarity
colnames(algmale) <- c("scaffold", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")
colnames(algfem) <- c("scaffold", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")

# Normalize data
algmale$normalized <- algmale$numreads / sum(algmale$numreads)
algfem$normalized <- algfem$numreads / sum(algfem$numreads)

# Calculate scaffold sizes
algmale$size <- algmale$endpos - algmale$startpos
algfem$size <- algfem$endpos - algfem$startpos

# Compute log2 ratio of female/male normalized reads
algmale$log2fem_male <- log2(algfem$normalized / algmale$normalized)

# List of scaffolds to highlight
#highlight_scaffolds <- c("ptg000009l", "ptg000018l", "ptg000019l",
#                         "ptg000022l", "ptg000028l", "ptg000029l",
#                         "ptg000038l", "ptg000046l", "ptg000098l")

# Add a column to indicate if a scaffold should be highlighted
#algmale$highlight <- ifelse(algmale$scaffold %in% highlight_scaffolds, "yes", "no")

ggplot(algmale[algmale$size > 10000, ], aes(x = size, y = log2fem_male)) +
  # Add points with conditional coloring
  geom_point(aes(color = highlight), size = 3, alpha = 0.7) +
  # Use a log10 scale for the x-axis
  scale_x_continuous(trans = "log10", labels = scales::comma_format()) +
  # Add horizontal reference lines
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_hline(yintercept = log2(0.3), color = "red", linetype = "dashed") +
  geom_hline(yintercept = log2(2), color = "green", linetype = "dashed") +
  # Label only points below the red line
  geom_text_repel(aes(label = ifelse(log2fem_male < log2(0.3), scaffold, "")),
                  size = 5, 
                  box.padding = 0.4, 
                  point.padding = 0.4, 
                  max.overlaps = Inf) +
  # Customize color scale
  scale_color_manual(values = c("no" = "blue", "yes" = "blue"), guide = "none") +
  # Axis labels
  labs(
    x = "Scaffold Size (log10 scale)",
    y = "Log2(Female/Male) Normalized Reads",
    title = "Algonquin Normalized Coverage Using Chrom Quotient Concept"
  ) +
  # Set y-axis limits
  ylim(-7.5, 2.5) +
  # Apply clean theme
  theme_cowplot(font_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Center and bold title
    axis.title = element_text(size = 16),  # Larger axis titles
    axis.text = element_text(size = 14),   # Larger axis text
    panel.grid.major = element_line(color = "grey90", size = 0.5),  # Subtle grid lines
    legend.position = "none"               # Remove legend
  )
```
![alt text](image-37.png)
Points aligning around 1 should be X-linked, and points aligning around 0 should be autosomal. Anything under the red line is putative Y-linked. This is because females have 2 Xs when males have 1, and males have 1 Y while females have zero. You expect (when log2 tranformed) for the autosomal reads to cancel out and end up around 0.

### Method 2

Indexcov

#### Code:
```
goleft indexcov --directory ../indexcov_Dalg *.bam
```

## Putative Y scaffolds
![alt text](image-17.png) ![alt text](image-18.png)
![alt text](image-19.png) ![alt text](image-20.png)
![alt text](image-21.png) ![alt text](image-22.png)
![alt text](image-23.png) ![alt text](image-24.png)
![alt text](image-25.png) ![alt text](image-26.png)
![alt text](image-27.png) ![alt text](image-28.png)
![alt text](image-29.png) ![alt text](image-30.png)
![alt text](image-31.png) ![alt text](image-32.png)
![alt text](image-33.png) ![alt text](image-34.png)
![alt text](image-35.png) ![alt text](image-36.png)

### Confirmation of putative Y scaffolds
| Scaffold| Length| # of genes (unmasked)| # of genes (masked)| samtools coverage| indexcov| PCR|
|-------------:|:-----:|:-------------:|:---:|:---------------:|:-------------:|:----------:|
| ptg00050l| | | | y| n| |
| ptg00004l| | | | y| y| |
| ptg00032l| | | | y| y| |
| ptg00023l| | | | y| y| |
| ptg00027l| | | | y| y| |
| ptg00008l| | | | y| y| |
| ptg00017l| | | | y| y| |
| ptg00007l| | | | y| y| |
| ptg00018l| | | | y| y| |
| ptg00019l| | | | y| y| |
| ptg00033l| | | | n| y| |

## Locate and mask repeats with repeatmodelor and repeatmasker on the cluster

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

## Annotate with helixer

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

## Renaming Transcripts
### Using Nilanjan's method

Step 1: On the cluster

```
#!/bin/bash
#SBATCH --job-name=best_hits         # Job name
#SBATCH --output=best_hits.out	  # Standard output log
#SBATCH --error=best_hits.err        # Standard error log
#SBATCH --time=36:00:00               # Time limit (hh:mm:ss)
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # CPUs per task (adjust as needed)
#SBATCH --mem=32G                     # Memory allocation (adjust as needed)
#SBATCH --partition=kucg          # Partition name (adjust as needed)
#SBATCH --mail-type=END,FAIL,BEGIN     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tconway@ku.edu   # Where to send mail

# Load necessary modules (adjust to your cluster's configuration)
module load blast+
module load seqtk

# make a database
#makeblastdb -in Dmel_translation_clean.fasta -dbtype prot -out dmel_protein_database

# Run blastx
blastx -query $1_transcripts.fa \
       -db dmel_protein_database \
       -outfmt 6 \
       -evalue 1e-5 \
       -max_target_seqs 1 \
       -num_threads 4 \
       -out /home/t043c581/scratch/data/blast_$1_transcript.txt

# Extract best hit protein IDs
cut -f2 blast_$1_transcript.txt | sort | uniq > $1_best_hit_proteins.list

# retrieve protein sequences
seqtk subseq Dmel_translation_clean.fasta $1_best_hit_proteins.list > $1_best_hit_proteins.fa

# Reciprocal TBLASTN search
makeblastdb -in $1_transcripts.fa -dbtype nucl -out $1_transcripts_db

# Identify reciprocal best hits
tblastn -query $1_best_hit_proteins.fa -db $1_transcripts_db -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -out $1_blast_reciprocal.txt

# Identify reciprocal best hits
awk '{print $1"\t"$2}' $1_blast_reciprocal.txt > $1_forward_hits.txt
awk '{print $2"\t"$1}' $1_blast_reciprocal.txt > $1_reciprocal_hits.txt
sort forward_hits.txt reciprocal_hits.txt | sed  's/-P[ABCDEFGHIJKLMNOPQRSTUVWXYZ]//g'  | uniq > reciprocal_best_hits.txt

# Add a subscript if duplicate gene found
awk '{if(a[$2]++){print $1"\t"$2"."a[$2]}else{print $0}}' $1_reciprocal_best_hits.txt > $1_RBH.txt
awk '{print $0 ".1"}' $1_RBH.txt > $1_RBH2.txt

# Replace in the file intended such as .gtf file
awk 'NR==FNR { mapping[$1] = $2; next } { for (key in mapping) gsub(key, mapping[key]) } 1' $1_RBH2.txt $1_unmasked.gff > $1_temp.gff
awk 'NR==FNR { mapping[$1] = $2; next } { for (key in mapping) gsub(key, mapping[key]) } 1' $1_RBH.txt $1_temp.gff > $1_renamed.gff
```
