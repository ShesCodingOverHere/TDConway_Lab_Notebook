# _D. affinis_ Y Project
We are finding the elusive Y-Chromosome.


## Part 1: Aligning the female assembly to _D. pseudoobscura_, then the male assembly to the female
This can be found via this path: /home/runcklesslab/Taylor/AffinisAlignPseudoO
### Align female assembly to D. pseudoobscura with MUMMER (or other - could look into whole genome aligners)
We imported the fasta file from _D. pseudoobscura_ from NCBI and used MUMmer software to align with the female _D. affinis_ assembly.

#### Coding Steps:
- nucmer --maxgap=500 -mincluster=100 GCF_009870125.1_UCI_Dpse_MV25_genomic.fna Daffinis_Female_HiC.20220721.fa
- delta-filter -q -r Pseudo_Affinis.delta > Pseudo_Affinis.filter
- show-coords -r Pseudo_Affinis.filter > Pseudo_Affinis.coords

### Align male assembly to female assembly
Similar to using MUMmer above, I used the program to align the male _D. affinis_ assembly to the female _D. affinis_ assembly.

#### Coding Steps:
- nucmer --maxgap=500 -mincluster=100 Daffinis_Female_HiC.20220721.fa Daffinis_Male_HiC.20221222.fa
- delta-filter -q -r Female_Male_Affinis.delta > Female_Male_Affinis.filter
- show-coords -r Female_Male_Affinis.filter > Female_Male_Affinis.coords

### Males to D Pseudo
- nucmer –maxgap=500 -mincluster=100 /home/runcklesslab/Taylor/Short_Read_Mapping/GCF_009870125.1_UCI_Dpse_MV25_genomic.fna /home/runcklesslab/Taylor/Short_Read_Mapping/Daffinis_Male_HiC.20221222.fa
- delta-filter -q -r out.delta > Pseudo_Male.filter
- show-coords -r Pseudo_Male.filter > Pseudo_Male.coords

#### Output files
- Female_Male_Affinis.coords
- Pseudo_Affinis.coords

## Part 2: Map short reads from males and females to the male assembly
This can be found via this directory: /home/runcklesslab/Taylor/Short_Read_Mapping

### Mapping and sorting

#### Files used
- /mnt/data/Daffinis/RawReads/Fragments/affinis_female_frag_R1.fastq.gz
- /mnt/data/Daffinis/RawReads/Fragments/affinis_female_frag_R2.fastq.gz
- /mnt/data/Daffinis/RawReads/Fragments/affinis_male_frag_R1.fastq.gz
- /mnt/data/Daffinis/RawReads/Fragments/affinis_male_frag_R2.fastq.gz

#### Coding Steps
- Inside of this script --> MaleFemaleMapping.sh
	- bwa index Daffinis_Male_HiC.20221222.fa
	- bwa mem Daffinis_Male_HiC.20221222.fa affinis_male_frag_R1.fastq.gz affinis_male_frag_R2.fastq.gz | samtools view -hb -F 4 - | samtools sort - > maleShortReads2MaleAssembly.bam; samtools index maleShortReads2MaleAssembly.bam
	- bwa mem Daffinis_Male_HiC.20221222.fa affinis_female_frag_R1.fastq.gz affinis_female_frag_R2.fastq.gz | samtools view -hb -F 4 - | samtools sort - > femaleShortReads2MaleAssembly.bam; samtools index femaleShortReads2MaleAssembly.bam

- Ran as: nohup sh MaleFemaleMapping.sh &
### samtools coverage
#### Coding Steps
- This script --> SamtoolsCoverage.sh
	- samtools coverage /home/runcklesslab/Taylor/Short_Read_Mapping/femaleShortReads2MaleAssembly.bam >femaleShortReadsCoverage.out
	- samtools coverage /home/runcklesslab/Taylor/Short_Read_Mapping/maleShortReads2MaleAssembly.bam >maleShortReadsCoverage.out
- Ran as: nohup sh SamtoolsCoverage.sh &

#### Output files

- maleShortReadsCoverage.out
- femaleShortReadsCoverage.out




## So I made some ggplots



![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/MeanDepthPlotMalesToFemales.png)

- This is just a plot of the mean depth with the female on the Y-axis and males on the X-axis
- If it falls on the line, we assume it's an autosome
- Above the line it is likely associated with the X-chromosome
- Below the line- potentially Y-chromosome
- But there is a better way to do this plot

![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/BetterFitLineMeanDepthMaleToFemale.png)

- I ran this: both$f.normreads=both$f.numberreads/sum(both$f.numberreads)
both$m.normreads=both$m.numberreads/sum(both$m.numberreads)
to get normalized reads for a more accurate line. 
- Now what doesn't fall on the line are likely chromosomal

![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/RatioNormalizedReads.png)

- This plot is fun! It shows us that we have some outliers above(which we used blastn to determine they were bacteria contamination)
- The points aligning the 0 line are likely autosomal
- The points around the 1 line are most likely X related
- We are interested in the points way below everything because they could be Y related

![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/PotentialYCandidates.png)


- This is looking at the points below that we are interested in. 


### Looking at outliers from the GGplot and seeing what to keep


*both[log2(both$f.normreads/both$m.normreads)>(2.5),]*

*both[log2(both$f.normreads/both$m.normreads)<(-.8),]*

#### Get rid of because they aren't in drosophila:
#### All Bacteria
-  Consensus_Consensus_Consensus_disjointig_490_pilon_pilon_pilon_pilon_pilon__unscaffolded
- PGA_scaffold30__23_contigs__length_1144930
- Consensus_Consensus_Consensus_disjointig_1063_pilon_pilon_pilon_pilon_pilon__unscaffolded
- Consensus_Consensus_Consensus_disjointig_1539_pilon_pilon_pilon_pilon_pilon__unscaffolded
- Consensus_Consensus_Consensus_disjointig_1555_pilon_pilon_pilon_pilon_pilon__unscaffolded
- Consensus_Consensus_Consensus_disjointig_1633_pilon_pilon_pilon_pilon_pilon__unscaffolded
- Consensus_Consensus_Consensus_disjointig_490_pilon_pilon_pilon_pilon_pilon__unscaffolded
- Consensus_Consensus_Consensus_disjointig_861_pilon_pilon_pilon_pilon_pilon__unscaffolded

#### Throw out... It's yeast
- Consensus_Consensus_Consensus_disjointig_353_pilon_pilon_pilon_pilon_pilon__unscaffolded
- Consensus_Consensus_Consensus_disjointig_560_pilon_pilon_pilon_pilon_pilon__unscaffolded


### Keep These:
These are all on the negative side of the plot.
#### In Drosophila
- PGA_scaffold11__133_contigs__length_7237977
- PGA_scaffold31__9_contigs__length_634783
- PGA_scaffold33__196_contigs__length_9964445
- Consensus_Consensus_Consensus_disjointig_37_pilon_pilon_pilon_pilon_pilon__unscaffolded
- Consensus_Consensus_Consensus_disjointig_566_pilon_pilon_pilon_pilon_pilon__unscaffolded