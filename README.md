# AlleleSorting


The R scripts for automated sorting of alleles to distinct homeologs as done in [Šlenker et al. 2021](https://www.frontiersin.org/articles/). Currently, we are able to sort alleles of suspected allotetraploid. However, future developments will focus on resolving genomes of higher ploidy levels as well as autopolyploids.


## Identifying homeologs inherited from different parents

#### Find two pairs of alleles (sequences), in which the alleles are closest to each other within the pairs while more distant between the pairs
Interallelic distances are estimated from the branch lengths of the corresponding ML tree (computed by `cophenetic` function of R base package `stats`). Distances between alleles are computed as an average length of path connecting all possible allele pairs or trios. As next, we compare calculated distances. If an average distance between alleles within any two pairs is more than `between_homeolog_distance` (desired treshold, see below)-time shorter than the average distance between alleles within any other possible arrangements, these pairs of alleles are considered to be unequivocally different and attributable to different homeologs. The treatment of each sample is witen into logfile. 

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

#### Sorting is succesful (uneqivocal) -> attribute to one parent
If the sorting is succesful, we have two pairs. We dont kow where they belong, only think we know is, that they are sufficiently similar within pairs and different between pairs.  
The next step is to calculate distance of each pair to a potential parent (`Parent_A`). The parir closer to `Parent_A` is labelled as homeolog A, and its alleles as A1 and A2 (in a random order). Pair more distant is labeled as homeolog B and its alleles as B1 and B2. The "`*** conditions are met`" is writen to log file.  
Arguments `use_between_parents_distance` (logical) and `between_parents_distance` (numeric) can be used for more precise filtering. If `use_between_parents_distance=TRUE`, homeolog A has to be more than `between_parents_distance` closer to `Parent_A` than to homeolog B, to consider them as uneqivocally sorted.  

&nbsp;  &nbsp;  

#### Sorting is NOT succesful  - mask/remove
If the sorting is not succesful, sample alleles can be masked by using Ns on the positions of intra allelic SNPs, or removed. Removing unsorted alleles is required for coalescent-based species tree estimation. Hovever, if we are able to concat sequences to longer blocks (to genes, for example), partially masked sequences can be concatenated with sorted sequences to prolong allignment length. In this case, partially masked sequences will not disrupt phylogenetic analysis. On the contrary, these sequences can bring better support on the speces level.  After concatenation, the allele sorting into two homeologs has to be verified again.  

&nbsp;  &nbsp;  

As a **result**, you will get allignment(s) where suffixes of sequences of sorted samples are renamed to A1, A2, B1 and B2, indicating their belonging to parental genomes. Unsorted samples will be removed if `if_bellow_treshold = "remove"` is used. If you chose to mask instead fo removing, unsorted samples will be renamed too, to allow following concatention. The treatment of each sample is witen into logfile.

The "masking" approach homogenize intra allelic variation, hovever, if we are able to concat sequences to longer blocks (to genes, for example), even masked sequences can provide additional support for sample/ species level differenciation.
After concatenation, the allele sorting into two homeologs was verified for each gene, to confirm unambiguity, or to remove the equivocal sample from the gene alignments.



&nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  




## The pipeline


#### Input files:
1) Aigned sequences in fasta format  
**tetraploid** sampes, which should be sorted to two homeologs have to be named as follow:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample>-a1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample>-a2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample>-a3  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample>-a4  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103-a1`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103-a2`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103-a3`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103-a4`  
where <4xSample> have to be listed in a sample file (see bullet 3).    
&nbsp;  
other sequences, or at leas sequences of sampes belonging to **diploid** Parent A have to be named as:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <2xSample>-a1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <2xSample>-a2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `ambC014_101-a1`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `ambC014_101-a2`  


2) Phylogenetic trees in newick format, bulit from alligned sequences.  

3) A sample file, containing a list of **tetraploid** samples, which alleles have to be sorted to homeologs. Each sample in a separate line.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample1>  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample2>  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC010_102`  


4) A species mapping file containing a list **diploid** species and individuals (samples). 
Not all samples need to be listed here. To sort allels of tetraploids, only samples of Parent A are required. Hovever, if you want to use other tools, e.g. to calculate distances of identified homeologs of tetraploid to all other diploid species, those other species have to be listed here. 
This mapping file should have one line per species, and each line needs to be in following format:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; species1:<2xSample1>,<2xSample2>  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `amara:ambC014_101,ambC024_103,ammC029_102`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `matthioli:matGRM9`  

5) A exon / gene list to process. Aligned sequence and phylogenetic tree built from particular sequence have to have the same basis of the file name. This is the list of those names. Script will iterate through this list and take sequence and corresponding phylogenetic tree.   


#### Pipeline Output:
1) Alligned sequences in fasta format. Alleles of tetraploids are either sorted or removed.
2) log file, containing a treatment of each sample for each sequence.

#### Running the pipeline

Use the `Main.R` script as a reference and edit it according your needs. Set arguments and load functions from `functions.R` into environment. 
```R
source('./functions.R')
```

Scripts depend on the following list of variables. You have to set them in the begining of the `Main.R` script. Please ensure they are set up correctly.  
* `samples_2x_File`: path to a species mapping file containing a list **diploid** species and samples (see above).  
* `samples_4x_File`: path to a file containing a list of **tetraploid** samples, which alleles have to be sorted to homeologs (see above).  
* `Parent_A`: A one of potential parent. Species name has to be included in `2xSamplesFile`. Can be represented by multiple samples.  
* `between_homeolog_distance`: two allele pairs are considered to be different homeologs inherited from different parents, if distance between this arrangement of alleles is `between_homeolog_distance`-shorter than distance of any other possible arrangement.  
* `use_between_parents_distance` and `between_parents_distance`: can be used for more precise filtering. If `use_between_parents_distance=TRUE`, homeolog A has to be more than `between_parents_distance` closer to Parent_A than to homeolog B, to consider them as uneqivocally sorted.  

* `treeDir`: path to a directory containing phylogenetic trees.  
* `seqDir`: path to a directory containing alligned sequences.  
* `patterns_exonsFile`, `patterns_genesFile`: A exon / gene list to process.  
* `if_bellow_treshold`: two options - "remove" or "mask", see above.  
* `ResDir`: a results dir for sorted sequences and log file.  
* `logfile`: log file, containing a treatment of each sample for each sequence.  


Alleles are sorted to A and B homeologs using a `sortAllels` function. Prefered approach is to use multiple threads, to speed up computation (using a `mclapply` function). A critical section in code is protected by exclusive locking. Alretnatively, use a for loop to do all in a single thread.  
&nbsp;  
If sorted sequences are concatenated, phylogenetic trees have to be computed and sorting have to be reevaluated again, using `confirmSorting` function.  

Distances of each homeolog pair to diploid samples can be calculated using `distanceTable` function.  





&nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  
&nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  
&nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  



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



