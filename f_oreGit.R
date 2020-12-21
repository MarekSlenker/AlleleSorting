


identitiHomeologs_by_1parent_and_mask_by_N <- function(pattern, treeDir, seqDir, resdir, samples, parent_A, if_bellow_treshold, logfile) {
  tocat = character();  tocat = c("\n", pattern, "\t")
  
  tree = read.tree(file = list.files(treeDir, pattern = pattern, full.names = TRUE))
  seq = read.dna(file = list.files(seqDir, pattern = pattern, full.names = TRUE), format = "fasta", as.character = T)
  
  # let's calculate all possible distances among alleles
  distanceTable1 = calculate_h1..h4_Distances(tree, samples)
  
  for (sample in samples$V1) {
    
    # are there any two pairs of alleles in which the alleles are closest to each other within the pairs while more threshold_distance distance between the pairs?
    if (  (distanceTable1[[sample]]$dist[2] / distanceTable1[[sample]]$dist[1] ) >  between_homeolog_distance &
          distanceTable1[[sample]]$twoOrThree[1] ==  "2pairs" &
          distanceTable1[[sample]]$twoOrThree[2] ==  "trio"
    ) {  
      
      # YES, they are
      
      homeolog_A.1 = paste(sample, distanceTable1[[sample]]$homeolog_A.1[1], sep = ""); homeolog_A.2 = paste(sample, distanceTable1[[sample]]$homeolog_A.2[1], sep = ""); homeolog_B.1 = paste(sample, distanceTable1[[sample]]$homeolog_B.1[1], sep = "");  homeolog_B.2 = paste(sample, distanceTable1[[sample]]$homeolog_B.2[1], sep = "")
      
      # let's calculate distances of two pairs to a potential parent
      dt_2 = distancesToParents(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, parent_A)
      
      if ( (dt_2[which(colnames(dt_2) == "dist_A_to_parent")] < dt_2[which(colnames(dt_2) == "dist_B_to_parent")]) ) {
        
        # use_between_parents_distance?
        if (use_between_parents_distance & (dt_2[which(colnames(dt_2) =="ratio")] < between_parents_distance)){
          tocat = c(tocat, "the_ratio_of_distances_of_homeologs_to_parents_is_lower_than_between_parents_distance=interallelic_SNPs_are_masked/removed", "\t")
          
          seq = do_if_bellow_treshold(if_bellow_treshold, sample, seq)
          # it doesn't matter how, but sequences have to be renamed
          rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
          
          next
        }
        
        colnames(dt_2) = c("A_homeolog", "B_homeolog", "ratio")
        
        # rename fasta sequence
        rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
        
      } else if ( dt_2[which(colnames(dt_2) == "dist_A_to_parent")] > dt_2[which(colnames(dt_2) == "dist_B_to_parent")] ){
        
        # use_between_parents_distance ?
        if (use_between_parents_distance & (dt_2[which(colnames(dt_2) =="ratio")] < between_parents_distance)){
          tocat = c(tocat, "the_ratio_of_distances_of_homeologs_to_parents_is_lower_than_between_parents_distance=interallelic_SNPs_are_masked", "\t")
          
          seq = do_if_bellow_treshold(if_bellow_treshold, sample, seq)
          
          # it doesn't matter how, but sequences have to be renamed
          rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
          
          next
        }
        colnames(dt_2) = c("B_homeolog", "A_homeolog", "ratio") 
        
        # rename fasta sequence
        rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-B2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-A2", sep = "")
        
      } else  { 
        
        tocat = c(tocat, "the_ratio_of_sample_distances_to_parents_is_equal=interallelic_SNPs_are_masked", "\t")
        
        seq = do_if_bellow_treshold(if_bellow_treshold, sample, seq)
        
        # it doesn't matter how, but sequences have to be renamed
        rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
        
        next
      }
      
      tocat = c(tocat, "*** conditions are met", "\t")
      
      
    } else {  # NO, it is not OK
      tocat = c(tocat, "the_ratio_of_distances_of_potential_homeologs_is_lower_than_between_homeolog_distance=interallelic_SNPs_are_masked", "\t")
      
      seq = do_if_bellow_treshold(if_bellow_treshold, sample, seq)
      
      # it doesn't matter how, but sequences have to be renamed
      
      homeolog_A.1 = paste(sample, "-A1", sep = ""); homeolog_A.2 = paste(sample, "-A2", sep = ""); homeolog_B.1 = paste(sample, "-B1", sep = "");  homeolog_B.2 = paste(sample, "-B2", sep = "")
      rownames(seq)[grep(paste(sample, "-h1 ", sep = ""), rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(paste(sample, "-h2 ", sep = ""), rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(paste(sample, "-h3 ", sep = ""), rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(paste(sample, "-h4 ", sep = ""), rownames(seq))] = paste(sample, "-B2", sep = "")
      
    }
    
    
  }
  
  write.dna(seq, file = paste(resdir, "/", pattern, ".AB_homeologs.fasta", sep = ""), format = "fasta", nbcol = -1, colsep = "")
  
  lck = lock(logfile, exclusive = TRUE, timeout = Inf)
  cat(tocat, file = logfile, append = T)
  unlock(lck)
}

do_if_bellow_treshold <- function(if_bellow_treshold, sample, seq) {
  switch (if_bellow_treshold,
          mask = {
            seq = mask_by_N(sample, seq)
            }
          , remove = {
            seq = seq[- grep(sample, rownames(seq)),]
          }  
  )
  
  return(seq)
}


calculate_h1..h4_Distances <- function(tree, samples) {
  
  dist = cophenetic(tree)  
  
  distTable = list()
  
  for (sample in samples$V1) {
    sh1 = paste(sample, "-h1", sep = "");   sh2 = paste(sample, "-h2", sep = "");   sh3 = paste(sample, "-h3", sep = "");   sh4 = paste(sample, "-h4", sep = "")
    
    dist12_34 = (dist[sh1,sh2] + dist[sh3,sh4])/2
    dist13_24 = (dist[sh1,sh3] + dist[sh2,sh4])/2
    dist14_23 = (dist[sh1,sh4] + dist[sh2,sh3])/2
    dist123 = (dist[sh1,sh2] + dist[sh1,sh3] + dist[sh2,sh3])/3 
    dist124 = (dist[sh1,sh2] + dist[sh1,sh4] + dist[sh2,sh4])/3 
    dist134 = (dist[sh1,sh3] + dist[sh1,sh4] + dist[sh4,sh3])/3 
    dist234 = (dist[sh3,sh2] + dist[sh4,sh3] + dist[sh2,sh4])/3 
    
    distTable_sample = data.frame(names=c("dist12_34", "dist13_24", "dist14_23", "dist123", "dist124", "dist134", "dist234"),
                                  dist = c(dist12_34, dist13_24, dist14_23, dist123, dist124, dist134, dist234),
                                  twoOrThree = c("2pairs", "2pairs", "2pairs", "trio", "trio", "trio", "trio"),
                                  homeolog_A.1 = c("-h1", "-h1", "-h1", NA, NA, NA, NA),
                                  homeolog_A.2 = c("-h2", "-h3", "-h4", NA, NA, NA, NA),
                                  homeolog_B.1 = c("-h3", "-h2", "-h2", NA, NA, NA, NA),
                                  homeolog_B.2 = c("-h4", "-h4", "-h3", NA, NA, NA, NA)
                                  
    )
    
    distTable[sample] = list( distTable_sample[order(distTable_sample$dist),])
  }
  
  return(distTable)
}


distancesToParents <- function(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, parent_A) {
  
  # tree$tip.label = sub("_$", "", tree$tip.label)
  
  dist = cophenetic(tree)  #teraz je to naozaj vzdialenost, tam po spolocny uzol i naspat
  distTable = data.frame(dist_A_to_parent = numeric(length = length(parent_A)), dist_B_to_parent = numeric(length = length(parent_A)), ratio = numeric(length = length(parent_A)))
  rownames(distTable) = parent_A
  
  
  for (parent in parent_A) {
    #break
    dist_A1_parent_h1 = dist[homeolog_A.1 , paste(parent, "-h1", sep = "")]
    dist_A1_parent_h2 = dist[homeolog_A.1 , paste(parent, "-h2", sep = "")]
    dist_A2_parent_h1 = dist[homeolog_A.2 , paste(parent, "-h1", sep = "")]
    dist_A2_parent_h2 = dist[homeolog_A.2 , paste(parent, "-h2", sep = "")]
      
    dist_B1_parent_h1 = dist[homeolog_B.1 , paste(parent, "-h1", sep = "")]
    dist_B1_parent_h2 = dist[homeolog_B.1 , paste(parent, "-h2", sep = "")]
    dist_B2_parent_h1 = dist[homeolog_B.2 , paste(parent, "-h1", sep = "")]
    dist_B2_parent_h2 = dist[homeolog_B.2 , paste(parent, "-h2", sep = "")]
      
    distTable[parent,]$dist_A_to_parent = mean(dist_A1_parent_h1, dist_A1_parent_h2, dist_A2_parent_h1, dist_A2_parent_h2)
    distTable[parent,]$dist_B_to_parent = mean(dist_B1_parent_h1, dist_B1_parent_h2, dist_B2_parent_h1, dist_B2_parent_h2)
    distTable[parent,]$ratio = max(distTable[parent,]$dist_A_to_parent, distTable[parent,]$dist_B_to_parent)/min(distTable[parent,]$dist_A_to_parent, distTable[parent,]$dist_B_to_parent)  
  }
  
  dt_2 = apply(distTable, 2, mean)
  dt_2 = t(dt_2)
  return(dt_2)
}

mask_by_N <- function(sample, seq) {
  
  samplePos = grep(sample, rownames(seq))
  
  x = apply(seq[samplePos,], 2, function(x) {
    length(unique(x))
  })
  
  seq[samplePos, which(x > 1)] = "n"
  
  return(seq)
}

