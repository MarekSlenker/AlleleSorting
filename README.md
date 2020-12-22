# AlleleSorting


The R scripts for automated sorting of alleles to distinct homeologs as done in [Šlenker et al. 2021](https://www.frontiersin.org/articles/). Currently, we are able to sort alleles of suspected allotetraploid. However, future developments will focus on resolving genomes of higher ploidy levels as well as autopolyploids.


### Software Dependencies
R: https://www.r-project.org/  


Following packages are required: `ape`, `seqinr`, and `filelock`. Within `R` command line started in `AlleleSorting` directory use e.g. command
```R
install.packages(pkgs=c("ape", "seqinr", "filelock"),
lib="rpackages", repos="https://mirrors.nic.cz/R/", dependencies="Imports")
```
to install needed packages.

### Installation
You can download this repository zipped (button on the right-hand side of the screen) or, if you have git installed on your system, clone it with:

```bash
git clone https://github.com/MarekSlenker/AlleleSorting.git
```



zarovnane sakvencie
ocakavany format, alely su oznacene ako <sample_name>-h1, <sample_name>-h2, <sample_name>-h3, <sample_name>-h4




pouzijeme **prvy skript**, dostane sekvencie roztriedene na A1 A2 B1 B2 a LOG file, ktory popisuje co sa s kazdou vzorkou stalo.
ak sme zvolili REMOVE,  vysledne sekvencie obsahuju len tie, ktore presli cez TRESHOLD  --> exon based analysis
ak sme zvolili MASK, sekvencie su makovane N. vysledne fasta subory obsahuju vsetky sekvencie, tie ktore nepresli treshold su maskovane - identicke. oznacenie je vsak zmenene ako A1-B2, co umoznije konkatenaciu napr cez AMAS


konkatenovane exony do genov treba odtestovat, ci splnaju TRESHOLD   **SKRIPT_2**   - daj 2 moznosti, no defaultne vyhod tie co nesplnaju.
vypise LOG

ak vyhadzujeme, dostavame parcialne datasety. Vsetko popisane v LOG ---> gene based alignments



geny
