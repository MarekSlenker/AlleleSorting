


treeDir="inputData/trees/"
seqDir="inputData/seq/"
parentsFile="inputData/parents"
samplesFile="inputData/samples"
Parent_A="amara"
patterns_exonsFile="inputData/patterns_exons"
patterns_genesFile="inputData/patterns_genes"
logfile = "xxx.log"
between_homeolog_distance = 4

use_between_parents_distance=FALSE
between_parents_distance = 2

if_bellow_treshold = "remove"  # "mask" or "remove"

AB_HomeologsDir="AB_Homeologs"
dir.create(AB_HomeologsDir)



library(ape)
library(seqinr)
library(filelock)
library(parallel)

options(scipen=999)



# read and parse parents, samples, and patterns
pp = read.delim(parentsFile, header = F, sep = ":", row.names = 1)
parents=list()
for (p in 1:length(row.names(pp))) {
  
  if (length(strsplit(pp$V2[p], split = ",")[[1]]) == 1) {
    eval( parse (text= paste("parents$", row.names(pp)[p], " = c(\"", pp$V2[p], "\")",sep = ""   )))
    next
  }
  eval( parse (text= paste("parents$", row.names(pp)[p], " = ", paste(strsplit(pp$V2[p], split = ","), sep = "\", \""),  sep = "") ))
}
rm (pp)

samples = read.delim(samplesFile, sep = ",", header = F)
patterns_exons = read.delim(patterns_exonsFile, header = FALSE,  stringsAsFactors = FALSE)
patterns_genes = read.delim(patterns_genesFile, header = FALSE,  stringsAsFactors = FALSE)
eval(parse( text =  paste("parent_A = parents$",Parent_A, sep = "")))





#############x
cat("", samples$V1, sep = "\t", file = logfile )
mclapply(patterns_exons$V1, identitiHomeologs_by_1parent_and_mask_by_N, treeDir=treeDir, seqDir = seqDir, resdir=AB_HomeologsDir, samples=samples, parent_A=parent_A, if_bellow_treshold=if_bellow_treshold, logfile=logfile, mc.cores = 4)




----------------------------------------------------------------

pattern=patterns_exons$V1[1]








mclapply(2:4, function(i,j,k) c(i,j,k), i=1, k=5)



for (n in names(sample_distance_table)) {
  write.table(t(sample_distance_table[[n]]),file = paste("sample_distance_table.",n, ".txt", sep = "")) 
}










for (pattern in patterns$V1) {
   cat("processing", pattern, "\n")
  
  
  
  tree = read.tree(file = list.files(args[1], pattern = pattern, full.names = TRUE))
  seq = read.dna(file = list.files(args[2], pattern = pattern, full.names = TRUE), format = "fasta", as.character = T)
  #seq=read.fasta(file = list.files(args[2], pattern = pattern, full.names = TRUE), seqtype = "DNA", as.string = T)
  

    distanceTable1 = calculateDistances(tree, samples)

    str(distanceTable1)

}


calculateDistances <- function(tree, samples) {
  
  dist = cophenetic(tree)  #teraz je to naozaj vzdialenost, tam po spolocny uzol i naspat
  #distTable_1 = array(dim = c(7 , 7, length(samples$V1)))
  distTable_1 = list()
  
  for (sample in samples$V1) {
    dist12_34 = (dist[paste(sample, "-h1", sep = ""),paste(sample, "-h2", sep = "")] + dist[paste(sample, "-h3", sep = ""),paste(sample, "-h4", sep = "")])/2
    dist13_24 = (dist[paste(sample, "-h1", sep = ""),paste(sample, "-h3", sep = "")] + dist[paste(sample, "-h2", sep = ""),paste(sample, "-h4", sep = "")])/2
    dist14_23 = (dist[paste(sample, "-h1", sep = ""),paste(sample, "-h4", sep = "")] + dist[paste(sample, "-h2", sep = ""),paste(sample, "-h3", sep = "")])/2
    dist123 = (dist[paste(sample, "-h1", sep = ""),paste(sample, "-h2", sep = "")] + dist[paste(sample, "-h1", sep = ""),paste(sample, "-h3", sep = "")] + dist[paste(sample, "-h2", sep = ""),paste(sample, "-h3", sep = "")])/3 
    dist124 = (dist[paste(sample, "-h1", sep = ""),paste(sample, "-h2", sep = "")] + dist[paste(sample, "-h1", sep = ""),paste(sample, "-h4", sep = "")] + dist[paste(sample, "-h2", sep = ""),paste(sample, "-h4", sep = "")])/3 
    dist134 = (dist[paste(sample, "-h1", sep = ""),paste(sample, "-h3", sep = "")] + dist[paste(sample, "-h1", sep = ""),paste(sample, "-h4", sep = "")] + dist[paste(sample, "-h4", sep = ""),paste(sample, "-h3", sep = "")])/3 
    dist234 = (dist[paste(sample, "-h3", sep = ""),paste(sample, "-h2", sep = "")] + dist[paste(sample, "-h4", sep = ""),paste(sample, "-h3", sep = "")] + dist[paste(sample, "-h2", sep = ""),paste(sample, "-h4", sep = "")])/3 
    
    distTable_sample = data.frame(names=c("dist12_34", "dist13_24", "dist14_23", "dist123", "dist124", "dist134", "dist234"),
                            dist = c(dist12_34, dist13_24, dist14_23, dist123, dist124, dist134, dist234),
                            twoOrThree = c("2pairs", "2pairs", "2pairs", "trio", "trio", "trio", "trio"),
                            parent_1.1 = c("-h1", "-h1", "-h1", NA, NA, NA, NA),
                            parent_1.2 = c("-h2", "-h3", "-h4", NA, NA, NA, NA),
                            parent_2.1 = c("-h3", "-h2", "-h2", NA, NA, NA, NA),
                            parent_2.2 = c("-h4", "-h4", "-h3", NA, NA, NA, NA)
                            
                       )
  
    
    #distTable_1[which(sample == samples$V1)] = list(distTable_sample)
    distTable_1[sample] = list( distTable_sample[order(distTable_sample$dist),])
    
  }
  
  return(distTable_1)
}












seq[]

seq["ulUD4-h1"]

rownames(seq)


  
treeFiles = list.files(path = args[1], full.names = T, pattern = ".support")
cat("tree\tsample\tshortestDist\t2ndShortestDist\tdistRatio\t\tsample\tshortestDist\t2ndShortestDist\tdistRatio\t\tsample\tshortestDist\t2ndShortestDist\tdistRatio\t\tsample\tshortestDist\t2ndShortestDist\tdistRatio\t\t", file=args[2], append = F)

for (t in treeFiles) {
  cat(t, "\n")
  tree  = read.tree(file = t)   #strom ulozeny z mega, dlzky etiev = div. time <-- stary koment, spred X rokov :)))
  dist = cophenetic(tree)  #teraz je to naozaj vzdialenost, tam po spolocny uzol i naspat

  cat(paste("\n",t, "\t", sep = ""), file=args[2], append=TRUE)
  for (sample in 4:(length(args))) {
    dist12_34 = (dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h2", sep = "")] + dist[paste(args[sample], "-h3", sep = ""),paste(args[sample], "-h4", sep = "")])/2
    dist13_24 = (dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h3", sep = "")] + dist[paste(args[sample], "-h2", sep = ""),paste(args[sample], "-h4", sep = "")])/2
    dist14_23 = (dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h4", sep = "")] + dist[paste(args[sample], "-h2", sep = ""),paste(args[sample], "-h3", sep = "")])/2
  
    dist123 = (dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h2", sep = "")] + 
               dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h3", sep = "")] +
               dist[paste(args[sample], "-h2", sep = ""),paste(args[sample], "-h3", sep = "")])/3 
    dist124 = (dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h2", sep = "")] + 
               dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h4", sep = "")] +
               dist[paste(args[sample], "-h2", sep = ""),paste(args[sample], "-h4", sep = "")])/3 
    dist134 = (dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h3", sep = "")] + 
               dist[paste(args[sample], "-h1", sep = ""),paste(args[sample], "-h4", sep = "")] +
               dist[paste(args[sample], "-h4", sep = ""),paste(args[sample], "-h3", sep = "")])/3 
    dist234 = (dist[paste(args[sample], "-h3", sep = ""),paste(args[sample], "-h2", sep = "")] + 
               dist[paste(args[sample], "-h4", sep = ""),paste(args[sample], "-h3", sep = "")] +
               dist[paste(args[sample], "-h2", sep = ""),paste(args[sample], "-h4", sep = "")])/3 
  
    dists = data.frame(d = c(dist12_34, dist13_24, dist14_23, dist123, dist124, dist134, dist234),
                     names=c("dist12_34", "dist13_24", "dist14_23", "dist123", "dist124", "dist134", "dist234"))
    dists = dists[order(dists$d),]
    #dists
    #minPos = which(dists$d %in% min(dists$d))
    
    # WRITE: 
    cat(paste(args[sample], 
            #length(minPos), 
            dists$names[1],
            dists$names[2],
            dists$d[2]/dists$d[1],
            #paste(dists$d[-c(1,2)], collapse = "\t"),
            "", "",
            sep = "\t" ), file=args[2], append=TRUE)
  
  }

  #cat("\n", file=args[2], append=TRUE)

}

