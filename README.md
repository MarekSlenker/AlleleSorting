# AlleleSorting


The R scripts for automated sorting of alleles to distinct homeologs as done in [Šlenker et al. 2021](https://www.frontiersin.org/articles/). Currently, we are able to sort alleles of suspected allotetraploid. However, future developments will focus on resolving genomes of higher ploidy levels as well as autopolyploids.


## Identifying homeologs inherited from different parents

#### Find two pairs of alleles (sequences), in which the alleles are closest to each other within the pairs while more distant between the pairs
Interallelic distances are estimated from the branch lengths of the corresponding ML tree (computed by `cophenetic` function of R base package `stats`). Distances between alleles are computed as an average length of parh connecting all possible allele pairs or trios. As next, we compare calculated distances. If an average distance between alleles within any two pairs is more than `between_homeolog_distance` (desired treshold)-time shorter than the average distance between alleles within any other possible arrangements, these pairs of alleles are considered to be unequivocally different and attributable to different homeologs.  

&nbsp;  

| distance | formula | tree A | tree B | tree C | tree D | tree E |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|  12_34  | **average** [ (a1_a2) + (a3_a4) ]  | 2  | 10  |  14  | 10  |  6  |
|  13_24  | **average** [ (a1_a3) + (a2_a4) ]  | 20  | 20  | 20  |  16 |  8  |
|  14_23  | **average** [ (a1_a4) + (a2_a3) ]  | 20  | 20  | 20  |  16 |  8  |
|  123  | **average** [ (a1_a2) + (a1_a3) + (a2_a3) ]  | 14  | 19.33  | 19.33  | 14  |  2.66  |
|  124  | **average** [ (a1_a2) + (a1_a4) + (a2_a4) ]  | 14  | 19.33  | 19.33  | 14  |  8.66  |
|  134  | **average** [ (a1_a3) + (a1_a4) + (a3_a4) ]  | 14  | 14  | 16.66  | 11.33  |  8.66  |
|  234  | **average** [ (a2_a3) + (a2_a4) + (a3_a4) ]  | 14  | 14  | 16.66  | 16.66  |   8.66 |



![trees](/trees.png)

lets set `between_homeolog_distance = 4`  
**Tree A:** &nbsp;&nbsp; The shortest distance (2) was found in arrangement as a1 and a2 as a a one pair and a3 and a4 as a another pair. The second shortest distance was found in trios (14). This distance is more than 4-time longer (between_homeolog_distance), thus we consider pairs a1-a2 and a3-a4 as unequivocally different and attributable to different homeologs (in the following steps).  
**Trees B, C, D:** &nbsp;&nbsp; The shortest distances are not 4-time longer than any other possible distance. These alleles can not be equivocally sorted to homeologs.  
**Tree E:** &nbsp;&nbsp; The shortest distance was found in trio, what interrupts any further attempt to find two pairs of distinct alleles.  



&nbsp;  &nbsp; 

#### Whwn sorting is succesful (uneqivocal)

&nbsp;  &nbsp;  


#### Whwn sorting is NOT succesful


&nbsp;  &nbsp;  &nbsp;  &nbsp;  


### Najdi rodicov

pouzijeme **prvy skript**, dostane sekvencie roztriedene na A1 A2 B1 B2 a LOG file, ktory popisuje co sa s kazdou vzorkou stalo.
ak sme zvolili REMOVE,  vysledne sekvencie obsahuju len tie, ktore presli cez TRESHOLD  --> exon based analysis
ak sme zvolili MASK, sekvencie su makovane N. vysledne fasta subory obsahuju vsetky sekvencie, tie ktore nepresli treshold su maskovane - identicke. oznacenie je vsak zmenene ako A1-B2, co umoznije konkatenaciu napr cez AMAS


&nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  







Following arangements of alleles can be found: 





The optimal threshold for unequivocal allele sorting was set to 4. This means that if an average distance between alleles within the proposed two pairs was more than four-time shorter than the average distance between alleles within any other possible arrangements, these pairs of alleles were considered to be unequivocally different and attributable to different homeologs. If the allele sorting did not pass the desired threshold, interallelic SNPs were masked by using Ns on those positions (for further concatenation to gene alignments, see below) or the sample was removed from the alignment (for exon-based analyses). As next, the allele pairs were attributed to different homeologs and labelled by calculating their distances to the alleles of all diploid species. The allele pair that was closer to C. amara (proposed as the maternal parent according to the plastid phylogeny, see below) was marked as homoelog ‘A’ and the other as homeolog ‘B’. Gene alignments were also assembled, in which the respective exons were concatenated to genes to obtain longer alignments with potentially higher phylogenetic signal. The concatenated exons included those with successfully sorted alleles into ‘A’ and ‘B’ homeologs, and those for which allele sorting was equivocal, with masked interallelic SNPs. After exon concatenation, the allele sorting into two homeologs was verified for each gene, with the same threshold as set for the exons above, to confirm unambiguity, or to remove the equivocal sample from the gene alignments. Both exon-based and gene-based alignments were used for species tree inference in ASTRAL-III. The labelled homeologs, representing the two subgenomes within C. barbaraeoides, were treated as independent accessions. Used scripts are available at http://github.com/MarekSlenker/AlleleSorting...(?). 


### Software Dependencies
R: https://www.r-project.org/  


Following packages are required: `ape`, `seqinr`, and `filelock`. Within `R` command line use e.g. command
```R
install.packages(pkgs=c("ape", "seqinr", "filelock"), repos="https://mirrors.nic.cz/R/", dependencies="Imports")
```
to install needed packages.

### Installation
You can download this repository zipped (button on the right-hand side of the screen) or, if you have git installed on your system, clone it with:

```bash
git clone https://github.com/MarekSlenker/AlleleSorting.git
```

### How to cite

**Use following formula:** Allele sequences were sorted to homeologs using the scripts and following the workflow available at https://github.com/MarekSlenker/AlleleSorting, described in detail by [Šlenker et al. (2021)](https://www.frontiersin.org/articles/).

# END



##################################################################################xx


zarovnane sakvencie
ocakavany format, alely su oznacene ako <sample_name>-h1, <sample_name>-h2, <sample_name>-h3, <sample_name>-h4




pouzijeme **prvy skript**, dostane sekvencie roztriedene na A1 A2 B1 B2 a LOG file, ktory popisuje co sa s kazdou vzorkou stalo.
ak sme zvolili REMOVE,  vysledne sekvencie obsahuju len tie, ktore presli cez TRESHOLD  --> exon based analysis
ak sme zvolili MASK, sekvencie su makovane N. vysledne fasta subory obsahuju vsetky sekvencie, tie ktore nepresli treshold su maskovane - identicke. oznacenie je vsak zmenene ako A1-B2, co umoznije konkatenaciu napr cez AMAS


konkatenovane exony do genov treba odtestovat, ci splnaju TRESHOLD   **SKRIPT_2**   - daj 2 moznosti, no defaultne vyhod tie co nesplnaju.
vypise LOG

ak vyhadzujeme, dostavame parcialne datasety. Vsetko popisane v LOG ---> gene based alignments



geny








| distance | formula | tree A | tree B | tree C | tree D | tree E |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|  dist12_34  | **average** [ (a1_a2) + (a3_a4) ]  | 2  | 10  |  14  | 10  |  6  |
|  dist13_24  | **average** [ (a1_a3) + (a2_a4) ]  | 20  | 20  | 20  |  16 |  8  |
|  dist14_23  | **average** [ (a1_a4) + (a2_a3) ]  | 20  | 20  | 20  |  16 |  8  |
|  dist123  | **average** [ (a1_a2) + (a1_a3) + (a2_a3) ]  | 2+20+20=14  | 18+20+20=19.33  | 19.33  | 10+12+20=14  |  2.66  |
|  dist124  | **average** [ (a1_a2) + (a1_a4) + (a2_a4) ]  | 14  | 19.33  | 19.33  | 14  |  8.66  |
|  dist134  | **average** [ (a1_a3) + (a1_a4) + (a3_a4) ]  | 14  | 14  | 20+20+10=16.66  | 12+12+10=11.33  |  8.66  |
|  dist234  | **average** [ (a2_a3) + (a2_a4) + (a3_a4) ]  | 14  | 14  | 16.66  | 16.66  |   8.66 |



| distance | formula | tree A | tree B | tree C |
| ------------- | ------------- | ------------- | ------------- | ------------- |
|  dist12_34  | **average** [ (a1_a2) + (a3_a4) ]  | 2  | 11  | 20  |
|  dist13_24  | **average** [ (a1_a3) + (a2_a4) ]  | 22  | 22  | 22  |
|  dist14_23  | **average** [ (a1_a4) + (a2_a3) ]  | 22  | 22  | 22  |

|  dist123  | **average** [ (a1_a2) + (a1_a3) + (a2_a3) ]  | 2+22+22=15.33  | 20+22+22=21.33  | 20+22+22=21.33  |
|  dist124  | **average** [ (a1_a2) + (a1_a4) + (a2_a4) ]  | 15.33  | 21.33  | 21.33  |
|  dist134  | **average** [ (a1_a3) + (a1_a4) + (a3_a4) ]  | 15.33  | 21.33 | 21.33  |
|  dist234  | **average** [ (a2_a3) + (a2_a4) + (a3_a4) ]  | 15.33  | 15.33  | 21.33  |



