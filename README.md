# AlleleSorting


The R scripts for automated sorting of alleles obtained from polyploid genomes to distinct homeologs as done in [Šlenker et al. 2021](https://www.frontiersin.org/articles/10.3389/fpls.2021.659275). Currently, we are able to sort alleles of suspected allotetraploids. However, future developments will focus on resolving genomes of higher ploidy levels as well as autopolyploids.


## Identifying homeologs inherited from different parents

#### Find two pairs of alleles (sequences), in which the alleles are closest to each other within the pairs while more distant between the pairs
First, interallelic distances are computed from the branch lengths of the corresponding ML tree (computed by `cophenetic` function of R base package `stats`). Distances between alleles are computed as an average length of path connecting all possible allele pairs or trios. As next, we compare the computed distances. If an average distance between alleles within any two pairs is more than `between_homeolog_distance` (the desired threshold, see below)-time shorter than the average distance between alleles within any other possible arrangements, these pairs of alleles are considered to be unequivocally different and attributable to different homeologs. The treatment of each sample is written into a logfile. 

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



![trees](/trees3.png)

lets set `between_homeolog_distance = 4`  
**Tree A:** &nbsp;&nbsp; The shortest distance (2) was found in the arrangement as a1 and a2 as a one pair and a3 and a4 as another pair. The second shortest distance was found in trios (14). This distance is more than 4-time longer (between_homeolog_distance), thus we consider pairs a1-a2 and a3-a4 as unequivocally different and attributable to different homeologs (in the following steps).  
**Trees B, C, D:** &nbsp;&nbsp; The shortest distance is not a 4-time shorter than the distance of any other possible arrangement. These alleles cannot be equivocally sorted to homeologs.  
**Tree E:** &nbsp;&nbsp; The shortest distance was found in a trio, what interrupts any further attempt to find two pairs of distinct alleles.  

**For more details about searching for the optimal `between_homeolog_distance` value, see Supplementary Text 1 (https://www.frontiersin.org/articles/10.3389/fpls.2021.659275/full#supplementary-material).**


&nbsp;  &nbsp; 

#### Sorting is successful (unequivocal) -> attribute to one parent
If the sorting is successful, we have two pairs. We don't know where they belong (to which subgenomes/parental genomes), the only thing we know is, that they are sufficiently similar within pairs and different between pairs.  
The next step is to calculate the distance of each pair to a potential parent (`Parent_A`), e.g. the maternal parent suggested from cpDNA phylogeny. The pair closer to `Parent_A` is labelled as homeolog A, and its alleles as A1 and A2 (in random order). The pair more distant is labelled as homeolog B and its alleles as B1 and B2. The "`*** conditions are met`" is written to logfile.  
Arguments `use_between_parents_distance` (logical) and `between_parents_distance` (numeric) can be used for more precise filtering. If `use_between_parents_distance=TRUE`, homeolog A has to be more than `between_parents_distance` closer to `Parent_A` than to homeolog B, to consider them as unequivocally  sorted.  

&nbsp;  &nbsp;  

#### Sorting is NOT successful  - mask/remove
If the sorting is not successful, sample alleles can be masked by using Ns on the positions of intra-allelic SNPs, or removed. Removing unsorted alleles is required for coalescent-based species tree estimation. However, if we are able to concatenate sequences to longer blocks (to genes, for example), partially masked sequences can be concatenated with sorted sequences to prolong alignment length. In this case, partially masked sequences will not disrupt phylogenetic analysis. On the contrary, these sequences can bring better support on the species level.  After concatenation, the allele sorting into two homeologs has to be verified again.  

&nbsp;  &nbsp;  

As a **result**, you will get alignment(s) where suffixes of sequences of sorted samples are renamed to A1, A2, B1 and B2, indicating their assignment to parental genomes. Unsorted samples will be removed if `if_below_treshold = "remove"` is used. If you chose masking instead of removing, unsorted samples will be masked but renamed as A1, A2, B1 and B2 too, to allow following concatenation. The treatment of each sample is written into logfile.
The "masking" approach homogenizes inter-allelic variation, however, if we are able to concatenate sequences to longer blocks (to genes, for example), even masked sequences can provide additional support for sample/species-level differentiation.
The sorting into homeologs has to be verified after concatenation, to confirm unambiguity, or to remove the equivocal sample from the gene alignments.



&nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  




## The pipeline


#### Input files:
1) Aligned sequences in fasta format  
**tetraploid** samples, which should be sorted to two homeologs have to be named as follow:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample>-a1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample>-a2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample>-a3  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample>-a4  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103-a1`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103-a2`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103-a3`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103-a4`  
where <4xSample>s have to be listed in a sample file (see bullet 3).    
&nbsp;  
other sequences, or at least sequences of samples belonging to a **diploid** Parent A has to be named as:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <2xSample>-a1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <2xSample>-a2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `ambC014_101-a1`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `ambC014_101-a2`  


2) Phylogenetic trees in newick format, built from aligned sequences.  

3) A sample file, containing a list of **tetraploid** samples, which alleles have to be sorted to homeologs. Each sample in a separate line.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample1>  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <4xSample2>  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC007_103`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `barC010_102`  


4) A species mapping file containing a list of **diploid** species and individuals (samples). 
Not all samples need to be listed here. To sort alleles of tetraploids, only samples of Parent A are required. However, if you want to use other tools, e.g. to calculate distances between identified homeologs of tetraploids to other diploid species, these other species have to be listed here. 
This mapping file should have one line per species, and each line needs to be in the following format:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; species1:<2xSample1>,<2xSample2>  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `amara:ambC014_101,ambC024_103,ammC029_102`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `matthioli:matGRM9`  

5) An exon/gene list to process. Aligned sequences and the corresponding phylogenetic trees must have the same basis of the file name. This is the list of those names. The script will iterate through this list and take a sequence and corresponding phylogenetic tree.   


#### Pipeline Output:
1) Aligned sequences in fasta format. Alleles of tetraploids are either sorted or removed.
2) log file, containing a treatment of each sample for each sequence.

#### Running the pipeline

Use the `Main.R` script as a reference and edit it according to your needs. Set arguments and load functions from `functions.R` into the environment. 
```R
source('./functions.R')
```

Scripts depend on the following list of variables. You have to set them in the `Main.R` script. Please ensure they are set up correctly.  
* `samples_2x_File`: path to a species mapping file containing a list of **diploid** species and samples (see above).  
* `samples_4x_File`: path to a file containing a list of **tetraploid** samples, which alleles have to be sorted to homeologs (see above).  
* `Parent_A`: one of a potential parent. Species name has to be included in `samples_2x_File`. Can be represented by multiple samples (individuals) as specified in the species mapping file.  
* `between_homeolog_distance`: two allele pairs are considered to be different homeologs inherited from different parents if the distance between this arrangement of alleles is `between_homeolog_distance`-shorter than the distance of any other possible arrangement.  
* `use_between_parents_distance` and `between_parents_distance`: can be used for more precise filtering. If `use_between_parents_distance=TRUE`, homeolog A has to be more than `between_parents_distance` closer to Parent_A than to homeolog B, to consider them as unequivocally sorted.  

* `treeDir`: path to a directory containing phylogenetic trees.  
* `seqDir`: path to a directory containing aligned sequences.  
* `patterns_exonsFile`, `patterns_genesFile`: An exon/gene list to process.  
* `if_below_treshold`: two options - "remove" or "mask", see above.  
* `ResDir`: a results directory for sorted sequences and log file.  
* `logfile`: the log file, containing a treatment of each sample for each sequence.  


Alleles are sorted to A and B homeologs using a `sortAllels()` function. If sorted sequences are concatenated, phylogenetic trees have to be computed and sorting has to be reevaluated again, using `confirmSorting()` function. Distances of each homeolog pair to diploid samples can be calculated using `distanceTable()` function.  

The functions are implemented in a thread-safe way. A critical sections in code are protected by exclusive locking. Thus, multithreading (using a `mclapply()` function) is a preferred approach as it significantly speeds up computation time. Alternatively, use a `for()` loop to do all in a single thread.
&nbsp;  



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

**Use the following formula:** Allele sequences were sorted to homeologs using the scripts and following the workflow available at https://github.com/MarekSlenker/AlleleSorting, described in detail by [Šlenker et al. (2021)](https://www.frontiersin.org/articles/10.3389/fpls.2021.659275).

### Other questions not covered here and reporting problems
If you have a question or you encounter a problem, feel free to email me at marek.slenker@savba.sk, or use [issues](https://github.com/MarekSlenker/AlleleSorting/issues). I will do my best to help you.

