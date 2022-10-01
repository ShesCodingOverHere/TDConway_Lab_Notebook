# Plotting Depth for the Putative Y Scaffolds

We used python (code at the bottom of this post) to get the depth for each scaffold for male and female. Then we used R to plot those points.

## Plots for Each Scaffold

### Scaffold 11
#### Males minus females where male depth is greater than 10
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff11DepthMaleMinusFemaleGT10.png)

#### Males minus females where male depth is greater than 0
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff11DepthmaleminusfemaleGT0.png)

#### Male Depth > 10
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff11DepthMaleDividedByFemaleGT10.png)

#### Female Depth > 0
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff11DepthmaleminusfemaleGT0.png)

### Scaffold 31
#### Males minus females where male depth is greater than 10
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff31DepthMaleminusFemaleGT10.png)

#### Males minus females where male depth is greater than 0
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff31MaleMinusFemaleGT0.png)

#### Male Depth > 10
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff31DepthMaleDividedByFemaleGT10.png)

#### Scaffold 33
#### Male Depth > 10
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff33Depthmaleminusfemale.png)

#### Male Depth > 0
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff33SmoothDepthGT0.png)

#### Male Depth > 10
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff33Depthmaledividedbyfemale.png)

#### Scaffold 37(It's not really a scaffold, it's a disjoint??)
#### Male Depth > 10
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff37DepthmaleminusfemaleGT10.png)

#### Male Depth > 0
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff37DepthMaleMinusFemaleGT0.png)

#### Male Depth > 10
![](https://raw.githubusercontent.com/ShesCodingOverHere/TDConway_Lab_Notebook/master/images/Scaff37DepthMaleDividedByFemaleGT10.png)

### Python code:
``samtools depth -aa -r PGA_scaffold33__196_contigs__length_9964445 /home/runcklesslab/Taylor/Short_Read_Mapping/femaleShortReads2MaleAssembly.bam >Scaff33.female.depth``

``samtools depth -aa -r PGA_scaffold33__196_contigs__length_9964445 /home/runcklesslab/Taylor/Short_Read_Mapping/maleShortReads2MaleAssembly.bam >Scaff33.male.depth``

`wc -l *depth`

`paste *depth | cut -f 1,2,3,6 >Scaff33.combined.depth`

`samtools depth -aa -r PGA_scaffold11__133_contigs__length_7237977 /home/runcklesslab/Taylor/Short_Read_Mapping/femaleShortReads2MaleAssembly.bam >Scaff11.female.depth`

`samtools depth -aa -r PGA_scaffold11__133_contigs__length_7237977 /home/runcklesslab/Taylor/Short_Read_Mapping/maleShortReads2MaleAssembly.bam >Scaff11.male.depth`

`paste Scaff11* | cut -f 1,2,3,6 >Scaff11.combined.depth`

`samtools depth -aa -r PGA_scaffold31__9_contigs__length_634783 /home/runcklesslab/Taylor/Short_Read_Mapping/maleShortReads2MaleAssembly.bam >Scaff31.male.depth`

`samtools depth -aa -r PGA_scaffold31__9_contigs__length_634783 /home/runcklesslab/Taylor/Short_Read_Mapping/femaleShortReads2MaleAssembly.bam >Scaff31.female.depth`

`paste Scaff31* | cut -f 1,2,3,6 >Scaff31.combined.depth`

`samtools depth -aa -r Consensus_Consensus_Consensus_disjointig_37_pilon_pilon_pilon_pilon_pilon__unscaffolded /home/runcklesslab/Taylor/Short_Read_Mapping/maleShortReads2MaleAssembly.bam >Scaff37.male.depth`

`samtools depth -aa -r Consensus_Consensus_Consensus_disjointig_37_pilon_pilon_pilon_pilon_pilon__unscaffolded /home/runcklesslab/Taylor/Short_Read_Mapping/femaleShortReads2MaleAssembly.bam >Scaff37.female.depth`

`paste Scaff37* | cut -f 1,2,3,6 >Scaff37.combined.depth`
