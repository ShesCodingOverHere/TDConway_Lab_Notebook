# Finding Primers for Potential Ys

I plotted the 4 potential Ys against the female to see where they land using R
Script:C:/Users/tdcon/Desktop/KUResearch/PutativeYtoX.R

### These are the potentials with their plots
 - PGA_scaffold11__133_contigs__length_7237977 
 
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/PotentialY1.png)
 
 Without mtDNA:
 
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/PotYWithoutmtDNA.png)
 
 Scaffold 11 mapped mtDNA. We will look in the sequence for primers where there is missing mappage.
 
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/MTDNASCAFF11.png)
 - PGA_scaffold33__196_contigs__length_9964445
 
 ![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/PotentialY2.png)

 - PGA_scaffold31__9_contigs__length_634783
 
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/PotentialY3.png)

- Consensus_Consensus_Consensus_disjointig_37_pilon_pilon_pilon_pilon_pilon__unscaffolded

![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/PotentialY4.png)

### Command line work:

In the first potential, we found a big segment that partially mapped to the mtDNA but seemingly didn't map to anything for it's first 1000000ish bases. We used bioawk to pull out 10000 bases from that scaffold in the male HIC data that we assume is in that missing region to use to find primers.

`bioawk -c fastx '$name=="PGA_scaffold11__133_contigs__length_7237977" {print ">sc11.1000000\n" substr($seq,1000000,10000)}' /home/runcklesslab/Taylor/Short_Read_Mapping/Daffinis_Male_HiC.20221222.fa >MaleScaff11.1000000to1010000.fa`

Next we made a blast database for the female HiC data:

`makeblastdb -in /home/runcklesslab/Taylor/Short_Read_Mapping/Daffinis_Female_HiC.20220721.fa -out FemaleHiC -title FemaleHiC -dbtype nucl`

Used blast to compare the male scaffold to the female HiC blast file.

`blastn -query MaleScaff11.1000000to1010000.fa -db FemaleHiC -outfmt 6 >blastfile.out`

Made a file to hold the forward and reverse primer region to use to find a primer that fits what we need.

`MaleScaff1.1000000.primers`

`blastn -query MaleScaff1.1000000.primers -db FemaleHiC -outfmt 6 -task blastn-short`

Using this method to find primers isn't as effecient as we could be, so we took a different approach.

Tried to use this to unzip some tar files, but that did not work...

`tar -xvf /home/runcklesslab/Taylor/Transcriptome/affinisTranscriptomes.tar.gz`

Used this older transcriptome data to make a female blast file:

`blastn -query affinis_transcriptome.20120529.fasta -db ../Short_Read_Mapping/FemaleHiC2 -outfmt "6 std qlen" -max_target_seqs 1 > transcriptome.blast &`

This is a fun bit of code that let's you watch as a file creates itself.

`watch tail transcriptome.blast`

Found the length of the fasta file to compare to mapped bases we found in R:

`grep ">" affinis_transcriptome.20120529.fasta |wc -l`

We also made a male blast file:

`blastn -query affinis_transcriptome.20120529.fasta -db ../Short_Read_Mapping/MaleHiC -outfmt "6 std qlen" -max_target_seqs 1 > transcriptome.blast.male &`

This next 3 lines of code we used to look for our scaffold in question in the female and male blast files, and then in the male blast file with a percent ID greater than 99% and matched mapped bases greater than 80%

`grep PGA_scaffold11__133_contigs__length_7237977 transcriptome.blast`

`grep PGA_scaffold11__133_contigs__length_7237977 transcriptome.blast.male`

`grep PGA_scaffold11__133_contigs__length_7237977 transcriptome.blast.male |awk '$3>99 && $4/$13>0.8'`

From our findings there, we chose a match with a larger length and made it into its own fasta file:

`bioawk -cfastx '$name=="comp698032_c0_seq1" {print ">comp698032_c0_seq1\n" $seq}' affinis_transcriptome.20120529.fasta >comp698032_c0_seq1.fa`

Used blast query to check this new fasta's presence in both male and female:

`blastn -query comp698032_c0_seq1.fa -db ../Short_Read_Mapping/FemaleHiC2 -outfmt "6 std qlen" -max_target_seqs 1`

Redid those last few steps a couple of times, and will continue to do this until I find some good potential primers. I will eventually plug these shorter transcriptome bases into primer3.

`bioawk -cfastx '$name=="comp29667_c0_seq1" {print ">comp29667_c0_seq1\n" $seq}' affinis_transcriptome.20120529.fasta >comp29667_c0_seq1.fa`

`blastn -query comp29667_c0_seq1.fa -db ../Short_Read_Mapping/FemaleHiC2 -outfmt "6 std qlen" -max_target_seqs 1`

`blastn -query comp29667_c0_seq1.fa -db ../Short_Read_Mapping/MaleHiC -outfmt "6 std qlen" -max_target_seqs 1`

Similar to above, we grepped out a sequence within some parameters

`grep PGA_scaffold11__133_contigs__length_7237977 transcriptome.blast.male |awk '$3>99 && $4/$13>0.8 && $4>800'`

Transcriptome segments tried in those repeated steps:
 
| transcriptome |scaffold|
|--|--|
|comp29667_c0_seq1| 11 |
| comp29693_c0_seq1 | 11 |
| comp698032_c0_seq1 | 11 |
| comp30262_c0_seq1| 11|
| comp28193_c0_seq1| 31|
| comp174474_c0_seq1| 31|

Not much luck here...

Started looking at scaffold 31:

`grep PGA_scaffold31__9_contigs__length_634783 transcriptome.blast.male`

## A different approach
We made a table of unique hits in the male that don't exist in the female

The code is as follows:

    cut -f1 transcriptome.blast |sort |uniq >transcriptome.femalehits.uniq
    cut -f1 transcriptome.blast.male |sort |uniq >transcriptome.malehits.uniq
    while read p ; do echo -n $p " " ; grep -c $p transcriptome.femalehits.uniq ; done < transcriptome.malehits.uniq |awk '$2==0 {print $1}' >MaleOnlyTranscriptomeHits.out
    while read p ; do grep $p transcriptome.blast.male ; done < MaleOnlyTranscriptomeHits.out |awk '{c[$2]++}END{for (i in c) printf("%s\t%s\n",i,c[i])}'


|Scaffold |Matches |
|--|--|
|PGA_scaffold10__20_contigs__length_963180 |87 |
|PGA_scaffold1__21_contigs__length_49891873 |23 |
|PGA_scaffold11__133_contigs__length_7237977 |49 |
|PGA_scaffold6__31_contigs__length_4556253 |1274 |
|PGA_scaffold2__103_contigs__length_47279294 |56 |
|PGA_scaffold33__196_contigs__length_9964445 |183 |
|PGA_scaffold30__23_contigs__length_1144930 |29 |

---

|Transcript |Scaffold |Primer |Forward |Reverse |PCR amplification |
|--|--|--|--|--|--|
|comp1083250_c0_seq1|6|Yes|ATTCCATGATCAGGTTATTGCCATT|ATATCATCAATGTCTGTTGTTTCGC|  |
|comp1326310_c0_seq1|11|Yes|AAGCACCGTGTACTTCTAAGCAG|GCTGAAGTTGTTGATGCTGTTG|  |
|comp29437_c0_seq1|10|Yes|CCACACTTTCATGGTTTGTATTCAC|CATTTTCCACGTTCTAGCATTCAAG|  |
