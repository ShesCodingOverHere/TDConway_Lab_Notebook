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


