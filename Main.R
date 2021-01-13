source('./functions.R')
options(scipen=999)


# GLOBAL atguments
samples_2x_File="test_dataset/inputData/2xSamples"
samples_4x_File="test_dataset/inputData/4xSamples"
Parent_A="amara"
between_homeolog_distance = 4
use_between_parents_distance=FALSE
between_parents_distance = 2

# read and parse parents, samples, and patterns
pp = read.delim(samples_2x_File, header = F, sep = ":", row.names = 1, stringsAsFactors = FALSE)
samples2x=list()
for (p in 1:length(row.names(pp))) {
  
  if (length(strsplit(pp$V2[p], split = ",")[[1]]) == 1) {
    eval( parse (text= paste("samples2x$", row.names(pp)[p], " = c(\"", pp$V2[p], "\")",sep = ""   )))
    next
  }
  eval( parse (text= paste("samples2x$", row.names(pp)[p], " = ", paste(strsplit(pp$V2[p], split = ","), sep = "\", \""),  sep = "") ))
}
rm (pp)
samples4x = read.delim(samples_4x_File, sep = ",", header = F, stringsAsFactors = FALSE)
# eval(parse( text =  paste("parent_A = samples2x$",Parent_A, sep = "")))



library(ape)
library(seqinr)
library(filelock)
library(parallel)




#  
# sortAllels 
#

# set local atguments
seqDir="test_dataset/inputData/seq"
treeDir="test_dataset/inputData/trees"
patterns_exonsFile="test_dataset/inputData/patterns_exons"
if_bellow_treshold = "mask"  # "mask" or "remove"
ResDir="test_dataset/1_AB_Homeologs_exons_seq"
logfile = "test_dataset/1_AB_Homeologs_exons_seq/1_AB_Homeologs_exons_seq.log"
dir.create(ResDir)


patterns_exons = read.delim(patterns_exonsFile, header = FALSE,  stringsAsFactors = FALSE)
cat("", samples4x$V1, sep = "\t", file = logfile)

# using multiple threads
mclapply(patterns_exons$V1, sortAllels, treeDir=treeDir, seqDir = seqDir, resdir=ResDir, samples2x=samples2x, samples4x=samples4x, Parent_A=Parent_A, if_bellow_treshold=if_bellow_treshold, between_homeolog_distance=between_homeolog_distance, logfile=logfile, mc.cores = 4)

# same in single thread
for (pattern in patterns_exons$V1) {
  cat(pattern, "\n")
  sortAllels(pattern, treeDir=treeDir, seqDir = seqDir, resdir=ResDir, samples2x=samples2x, samples4x=samples4x, Parent_A=Parent_A, if_bellow_treshold=if_bellow_treshold, between_homeolog_distance=between_homeolog_distance, logfile=logfile )
}


#
# Concatenation using AMAS (https://github.com/marekborowiec/AMAS)
# mkdir test_dataset/2_AB_Homeologs_genes_seq
# while read r; do AMAS.py concat -i test_dataset/1_AB_Homeologs_exons_seq/$r* -f fasta -d dna -t test_dataset/2_AB_Homeologs_genes_seq/$r.fasta; done < test_dataset/inputData/patterns_genes 
#


#  
# confirmSorting 
#

# set local atguments
patterns_genesFile="test_dataset/inputData/patterns_genes"
seqDir="test_dataset/2_AB_Homeologs_genes_seq"
treeDir="test_dataset/3_AB_Homeologs_genes_trees"
if_bellow_treshold = "remove"
ResDir="test_dataset/4_AB_Homeologs_genes_seq_confirmed"
logfile = "test_dataset/4_AB_Homeologs_genes_seq_confirmed/4_AB_Homeologs_genes_seq_confirmed.log"
dir.create(ResDir)



patterns_genes = read.delim(patterns_genesFile, header = FALSE,  stringsAsFactors = FALSE)
cat("", samples4x$V1, sep = "\t", file = logfile)

# using multiple threads
mclapply(patterns_genes$V1, confirmSorting, treeDir=treeDir, seqDir = seqDir, resdir=ResDir, samples2x=samples2x, samples4x=samples4x, Parent_A=Parent_A, if_bellow_treshold=if_bellow_treshold, between_homeolog_distance=between_homeolog_distance, logfile=logfile, mc.cores = 4)

# same in single thread
for (pattern in patterns_genes$V1) {
  cat(pattern, "\n")
  confirmSorting(pattern, treeDir=treeDir, seqDir = seqDir, resdir=ResDir, samples2x=samples2x, samples4x=samples4x, Parent_A=Parent_A, if_bellow_treshold=if_bellow_treshold, between_homeolog_distance=between_homeolog_distance, logfile=logfile )
}


#  
# distanceTable 
#

# set local atguments
patterns_genesFile="test_dataset/inputData/patterns_genes"
treeDir="test_dataset/5_AB_Homeologs_genes_seq_confirmed_trees"
tableFile = "test_dataset/table.txt"

patterns_genes = read.delim(patterns_genesFile, header = FALSE,  stringsAsFactors = FALSE)


cat("\t", paste(rep(samples4x$V1, each = length(names(samples2x))*2), collapse = "\t"), "\n", file = tableFile)
cat("sequence\t", rep(paste(paste(names(samples2x), collapse = "_to_A1A2\t"), "_to_A1A2\t", paste(names(samples2x), collapse = "_to_B1B2\t"), "_to_B1B2\t", sep = ""), nrow(samples4x) ), file = tableFile, append = TRUE)

# using multiple threads
mclapply(patterns_genes$V1, distanceTable, treeDir=treeDir, samples2x=samples2x, samples4x=samples4x, logfile=tableFile, mc.cores = 4)

# same in single thread
for (pattern in patterns_genes$V1) {
  cat(pattern, "\n")
  distanceTable(pattern, treeDir=treeDir, samples2x=samples2x, samples4x=samples4x, logfile=tableFile )
}



# END  -------------------------------

  