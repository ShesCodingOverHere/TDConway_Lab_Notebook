Path on the linux machine:
/home/runcklesslab/Taylor/DiscoverY

Using this as a guide:
[DiscoverY Github Repository](https://github.com/makovalab-psu/DiscoverY/blob/master/README.md)

This is the paper:
[DiscoverY Paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5996-3)

- I deleted their test data and named the male and female HiC files male.fasta and female.fasta respectively
- To make kmer files, I ran:
	- Females

    cd dependency
    ln -s /home/runcklesslab/Taylor/Short_Read_Mapping/Daffinis_Female_HiC_FixedScaffolds.fa
    ./run_dsk_Linux.sh /home/runcklesslab/Taylor/Short_Read_Mapping/Daffinis_Female_HiC_FixedScaffolds.fa 25
    
	- Males



Run with Melanogaster data:
Melanogaster data:
SRA: SRP007888
female-->SRR332298
male--> SRR332299

Turn the fastq into fasta

    sed -n '1~4s/^@/>/p;2~4p' /home/runcklesslab/Taylor/ChromosomeQuotient/SRR332299_male.fastq > male_contigs.fasta

    sed -n '1~4s/^@/>/p;2~4p' /home/runcklesslab/Taylor/ChromosomeQuotient/SRR332298_female.fastq > female.fasta


make the kmer files:

Female:

    ./run_dsk_Linux.sh /home/runcklesslab/Taylor/ChromosomeQuotient/SRR332298_female.fastq 25
Male:



### Run with their data...
Their data built into the program is not a full dataset.

I'm running it with the gorilla data by using the 4 fastq files from the paper and using a full genome of the western gorilla and concatenating the Y assembly to the rest of the assembly.

https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=555244

Their gorilla data will not finish for whatever reason. I have zero idea.
