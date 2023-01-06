
Directory: Y20221214

bwa index Daffinis_Female_plus_scaf6.33.11.fa  
bwa mem Daffinis_Female_plus_scaf6.33.11.fa ../Transcriptome/affinis_transcriptome.20120529.fasta | samtools view -hb - | samtools sort - > transcriptome.v.femaleplus.bam  
samtools index transcriptome.v.femaleplus.bam  
samtools tview transcriptome.v.femaleplus.bam Daffinis_Female_plus_scaf6.33.11.fa -p PGA_scaffold33__196_contigs__length_9964445:3143599

