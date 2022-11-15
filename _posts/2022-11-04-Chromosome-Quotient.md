https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-273#Sec18

I am running this without the use of the code in the paper.
I'll use both melanogaster(to check for accuracy) and affinis data.

 
Melanogaster data:
SRA: SRP007888
female-->SRR332298
male--> SRR332299
Go to SRA on NCBI and put in SRP007888.
Click "Send results to run selector"

To download:
fasterq-dump --split-files SRR332298
fasterq-dump --split-files SRR332299

Affinis data:
We used `wc -l` to determine how long our read files were.

Made a subset of our data using this code:

    seqtk sample -s 1234 affinis_female_frag_R2.fastq.gz 10000000 >affinis_female_R2_10mil.fastq
    seqtk sample -s 1234 affinis_male_frag_R2.fastq.gz 10000000 >affinis_male_R2_10mil.fastq
    seqtk sample -s 1234 affinis_male_frag_R1.fastq 10000000 >affinis_male_10mil.fastq
    seqtk sample -s 1234 affinis_female_frag_R1.fastq 10000000 >affinis_female_10mil.fastq


Shell script used to align to the reference:
/home/runcklesslab/Taylor/Short_Read_Mapping/MaleFemaleMapping.sh
