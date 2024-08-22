


sortAllels <- function(pattern, treeDir, seqDir, resdir, samples2x, samples4x, Parent_A, if_below_treshold, between_homeolog_distance, logfile) {
  tocat = character();  tocat = c("\n", pattern, "\t")
  
  tree = read.tree(file = list.files(treeDir, pattern = pattern, full.names = TRUE))
  seq = read.dna(file = list.files(seqDir, pattern = pattern, full.names = TRUE), format = "fasta", as.character = T)
  
  eval( parse (text= paste("Parent_A_list = list(", Parent_A, " = samples2x$", Parent_A, ")",  sep = "") ))
  

  # let's calculate all possible distances among alleles
  distanceTable1 = calculate_a1..a4_Distances(tree, samples4x)
  
  for (sample in samples4x$V1) {
    if (any(grepl(sample, names(distanceTable1)))) {

      homeolog_A.1 = paste(sample, distanceTable1[[sample]]$homeolog_A.1[1], sep = ""); homeolog_A.2 = paste(sample, distanceTable1[[sample]]$homeolog_A.2[1], sep = ""); homeolog_B.1 = paste(sample, distanceTable1[[sample]]$homeolog_B.1[1], sep = "");  homeolog_B.2 = paste(sample, distanceTable1[[sample]]$homeolog_B.2[1], sep = "")
      
      
      if (evaluateDistance (distanceTable1, sample, between_homeolog_distance)) {
        
        dt_2 = distancesToParentsList(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, Parent_A_list)
        
        toSolve = TRUE
        if ( dt_2$dist_A_to_parent < dt_2$dist_B_to_parent ) {
          # use_between_parents_distance?
          if (use_between_parents_distance & (dt_2$ratio < between_parents_distance)){
            tocat = c(tocat, "masked/removed", "\t")
            
            seq = do_if_below_treshold(if_below_treshold, sample, seq)
            # it doesn't matter how, but sequences have to be renamed
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
            
            toSolve = FALSE
          }
          if (toSolve) {
            tocat = c(tocat, "*** conditions are met", "\t")
            
            colnames(dt_2) = c("A_homeolog", "B_homeolog", "ratio")
            # rename fasta sequence
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
          }
        } else if (dt_2$dist_A_to_parent > dt_2$dist_B_to_parent ){
          # use_between_parents_distance ?
          if (use_between_parents_distance & (dt_2$ratio < between_parents_distance)){
            tocat = c(tocat, "masked/removed", "\t")
            
            seq = do_if_below_treshold(if_below_treshold, sample, seq)
            
            # it doesn't matter how, but sequences have to be renamed
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-B2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-A2", sep = "")
            
            toSolve = FALSE
          }
          if (toSolve) {
            tocat = c(tocat, "*** conditions are met", "\t")
            
            colnames(dt_2) = c("B_homeolog", "A_homeolog", "ratio") 
            # rename fasta sequence
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-B2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-A2", sep = "")
          }
        
        } else {
          tocat = c(tocat, "masked/removed", "\t")
          seq = do_if_below_treshold(if_below_treshold, sample, seq)
          
          rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
        }
      
      } else { # evaluateDistance = FALSE
        
        tocat = c(tocat, "masked/removed", "\t")
        
        seq = do_if_below_treshold(if_below_treshold, sample, seq)
        rownames(seq)[grep(paste(sample, "-a1", sep = ""), rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(paste(sample, "-a2", sep = ""), rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(paste(sample, "-a3", sep = ""), rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(paste(sample, "-a4", sep = ""), rownames(seq))] = paste(sample, "-B2", sep = "")
      }
    } else {
      tocat = c(tocat, "missing", "\t")
    }  
        
    ###########x
#    returnedList = evaluateDistance_sortSequence(distanceTable1, sample, between_homeolog_distance, tree, Parent_A_list, use_between_parents_distance, between_parents_distance, tocat, if_below_treshold, seq, homeolog_A.1, homeolog_A.2, homeolog_B.1,  homeolog_B.2)
#    seq = returnedList[[1]]
#    tocat = returnedList[[2]]
  }
  
  write.dna(seq, file = paste(resdir, "/", pattern, ".AB_homeologs.fasta", sep = ""), format = "fasta", nbcol = -1, colsep = "")
  
  lck = lock(logfile, exclusive = TRUE, timeout = Inf)
  cat(tocat, file = logfile, append = T)
  unlock(lck)
}

confirmSorting <- function(pattern, treeDir, seqDir, resdir, samples2x, samples4x, Parent_A, if_below_treshold, between_homeolog_distance, logfile) {
  tocat = character();  tocat = c("\n", pattern, "\t")
  
  tree = read.tree(file = list.files(treeDir, pattern = pattern, full.names = TRUE))
  seq = read.dna(file = list.files(seqDir, pattern = pattern, full.names = TRUE), format = "fasta", as.character = T)
  
  eval( parse (text= paste("Parent_A_list = list(", Parent_A, " = samples2x$", Parent_A, ")",  sep = "") ))
  
  # let's calculate all possible distances among alleles
  distanceTable1 = calculate_AB_Distances(tree, samples4x)
  
  for (sample in samples4x$V1) {
    
    if (length(grep(sample, tree$tip.label)) == 0) {
      next
    }
    
    # homeolog_A.1 = paste(sample, distanceTable1[[sample]]$homeolog_A.1[1], sep = ""); homeolog_A.2 = paste(sample, distanceTable1[[sample]]$homeolog_A.2[1], sep = ""); homeolog_B.1 = paste(sample, distanceTable1[[sample]]$homeolog_B.1[1], sep = "");  homeolog_B.2 = paste(sample, distanceTable1[[sample]]$homeolog_B.2[1], sep = "")
    # homeolog_A.1 = NA; homeolog_A.2 = NA; homeolog_B.1 = NA; homeolog_B.2 = NA
    
    if (evaluateDistance (distanceTable1, sample, between_homeolog_distance) ){
      
      tocat = c(tocat, "*** conditions are met", "\t")
      
    } else {
      tocat = c(tocat, "masked/removed", "\t")
      
      seq = do_if_below_treshold(if_below_treshold, sample, seq)
    }
    
    
    # I need to do it like this, as I need to return 2 objects
    # returnedList = evaluateDistance_sortSequence(distanceTable1, sample, between_homeolog_distance, tree, Parent_A_list, use_between_parents_distance, between_parents_distance, tocat, if_below_treshold, seq, homeolog_A.1, homeolog_A.2, homeolog_B.1,  homeolog_B.2)
    # seq = returnedList[[1]]
    # tocat = returnedList[[2]]
  }
  
  write.dna(seq, file = paste(resdir, "/", pattern, ".AB_homeologs.fasta", sep = ""), format = "fasta", nbcol = -1, colsep = "")
  
  lck = lock(logfile, exclusive = TRUE, timeout = Inf)
  cat(tocat, file = logfile, append = T)
  unlock(lck)
}

  
distanceTable <- function(pattern, treeDir, samples2x, samples4x, logfile) {
  distCat = character();  distCat = c("\n", pattern, "\t")
  
  tree = read.tree(file = list.files(treeDir, pattern = pattern, full.names = TRUE))
  
  
  for (sample in samples4x$V1) {
    
    # sample is missing in tree
    if (length(grep(sample, tree$tip.label)) == 0) {
      distCat = c(distCat, paste(rep("NA", length(samples2x)), collapse = "\t" ), "\t", paste(rep("NA", length(samples2x)), collapse = "\t" ), "\t")
      next
    }
    # is it OK? sure, you checked it earlier
    
    homeolog_A.1 = paste(sample, "-A1", sep = ""); homeolog_A.2 = paste(sample, "-A2", sep = ""); homeolog_B.1 = paste(sample, "-B1", sep = "");  homeolog_B.2 = paste(sample, "-B2", sep = "")
    
    dt_2 = distancesToParentsList(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, samples2x)
    
    distCat = c(distCat, paste(dt_2$dist_A_to_parent, collapse = "\t" ), "\t", paste(dt_2$dist_B_to_parent, collapse = "\t" ), "\t")
    
    
  }
  
  lck2 = lock(logfile, exclusive = TRUE, timeout = Inf)
  cat(distCat, file = logfile, append = T)
  unlock(lck2)
}

calculate_AB_Distances <- function(tree, samples4x) {
  
  dist = cophenetic(tree)  
  distTable_1 = list()
  
  
  for (sample in samples4x$V1) {
    if (length(grep(sample, tree$tip.label)) == 0) {
      next
    }
    
    
    A1 = paste(sample, "-A1", sep = "");   A2 = paste(sample, "-A2", sep = "");   B1 = paste(sample, "-B1", sep = "");   B2 = paste(sample, "-B2", sep = "")
    
    distA1A2_B1B2 = (dist[A1,A2] + dist[B1,B2])/2
    distA1B1_A2B2 = (dist[A1,B1] + dist[A2,B2])/2
    distA1B2_A2B1 = (dist[A1,B2] + dist[A2,B1])/2
    distA1A2B1 = (dist[A1,A2] + dist[A1,B1] + dist[A2,B1])/3 
    distA1A2B2 = (dist[A1,A2] + dist[A1,B2] + dist[A2,B2])/3 
    distA1B1B2 = (dist[A1,B1] + dist[A1,B2] + dist[B2,B1])/3 
    distA2B1B2 = (dist[B1,A2] + dist[B2,B1] + dist[A2,B2])/3 
    
    distTable_sample = data.frame(names=c("distA1A2_B1B2", "distA1B1_A2B2", "distA1B2_A2B1", "distA1A2B1", "distA1A2B2", "distA1B1B2", "distA2B1B2"),
                                  dist = c(distA1A2_B1B2, distA1B1_A2B2, distA1B2_A2B1, distA1A2B1, distA1A2B2, distA1B1B2, distA2B1B2),
                                  twoOrThree = c("2pairs", "2pairs", "2pairs", "trio", "trio", "trio", "trio"),
                                  homeolog_A.1 = c("-A1", "-A1", "-A1", NA, NA, NA, NA),
                                  homeolog_A.2 = c("-A2", "-B1", "-B2", NA, NA, NA, NA),
                                  homeolog_B.1 = c("-B1", "-A2", "-A2", NA, NA, NA, NA),
                                  homeolog_B.2 = c("-B2", "-B2", "-B1", NA, NA, NA, NA)
                                  
    )
    
    
    #distTable_1[which(sample == samples4x$V1)] = list(distTable_sample)
    distTable_1[sample] = list( distTable_sample[order(distTable_sample$dist),])
  }
  
  return(distTable_1)
}

do_if_below_treshold <- function(if_below_treshold, sample, seq) {
  switch (if_below_treshold,
          mask = {
            seq = mask_by_N(sample, seq)
            }
          , remove = {
            seq = seq[- grep(sample, rownames(seq)),]
          }  
  )
  
  return(seq)
}

calculate_a1..a4_Distances <- function(tree, samples4x) {
  
  dist = cophenetic(tree)  
  
  distTable = list()
  
  for (sample in samples4x$V1) {
    sh1 = paste(sample, "-a1", sep = "");   sh2 = paste(sample, "-a2", sep = "");   sh3 = paste(sample, "-a3", sep = "");   sh4 = paste(sample, "-a4", sep = "")
    if (any(grepl(sample, colnames(dist))))
    {
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
                                      homeolog_A.1 = c("-a1", "-a1", "-a1", NA, NA, NA, NA),
                                      homeolog_A.2 = c("-a2", "-a3", "-a4", NA, NA, NA, NA),
                                      homeolog_B.1 = c("-a3", "-a2", "-a2", NA, NA, NA, NA),
                                      homeolog_B.2 = c("-a4", "-a4", "-a3", NA, NA, NA, NA)
                                      
        )
    
        distTable[sample] = list( distTable_sample[order(distTable_sample$dist),])
      }
  }
  
  return(distTable)
}

distancesToParentsList <- function(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, Parent_A_list) {
  
  dist = cophenetic(tree)  
  #distTable = data.frame(dist_A_to_parent = numeric(length = length(Parent_A_list)), dist_B_to_parent = numeric(length = length(Parent_A_list)), ratio = numeric(length = length(Parent_A_list)))
  #rownames(distTable) = Parent_A_list
  
  
  distTable_2 = data.frame(dist_A_to_parent = numeric(length = length(Parent_A_list)), dist_B_to_parent = numeric(length = length(Parent_A_list)), ratio = numeric(length = length(Parent_A_list)))
  rownames(distTable_2) = names(Parent_A_list)
  
  
  
  for (parent in names(Parent_A_list)) {
    #break

    eval( parse (text= paste("parentSamples = Parent_A_list$", parent,  sep = "") ))
    
    parentSamplesDistTable = data.frame(dist_A_to_parent = numeric(length = length(parentSamples)), 
                                        dist_B_to_parent = numeric(length = length(parentSamples)), 
                                        ratio = numeric(length = length(parentSamples)))
    rownames(parentSamplesDistTable) = parentSamples
    
    for (parentSample in parentSamples) {
      if (any(grepl(parentSample, colnames(dist)))){

          dist_A1_parent_a1 = dist[homeolog_A.1 , paste(parentSample, "-a1", sep = "")]
          dist_A1_parent_a2 = dist[homeolog_A.1 , paste(parentSample, "-a2", sep = "")]
          dist_A2_parent_a1 = dist[homeolog_A.2 , paste(parentSample, "-a1", sep = "")]
          dist_A2_parent_a2 = dist[homeolog_A.2 , paste(parentSample, "-a2", sep = "")]
          
          dist_B1_parent_a1 = dist[homeolog_B.1 , paste(parentSample, "-a1", sep = "")]
          dist_B1_parent_a2 = dist[homeolog_B.1 , paste(parentSample, "-a2", sep = "")]
          dist_B2_parent_a1 = dist[homeolog_B.2 , paste(parentSample, "-a1", sep = "")]
          dist_B2_parent_a2 = dist[homeolog_B.2 , paste(parentSample, "-a2", sep = "")]
          
          parentSamplesDistTable[which(parentSample == parentSamples), 1] = mean( c( dist_A1_parent_a1, dist_A1_parent_a2, dist_A2_parent_a1, dist_A2_parent_a2)) 
          parentSamplesDistTable[which(parentSample == parentSamples), 2] = mean( c( dist_B1_parent_a1, dist_B1_parent_a2, dist_B2_parent_a1, dist_B2_parent_a2)) 
          parentSamplesDistTable[which(parentSample == parentSamples), 3] = parentSamplesDistTable[which(parentSample == parentSamples), 2]/parentSamplesDistTable[which(parentSample == parentSamples), 1]
        } else
        {
        parentSamplesDistTable[which(parentSample == parentSamples), 1] = NaN
        parentSamplesDistTable[which(parentSample == parentSamples), 2] = NaN
        parentSamplesDistTable[which(parentSample == parentSamples), 3] = NaN
        }
      }
    distTable_2[parent,]$dist_A_to_parent = mean(parentSamplesDistTable$dist_A_to_parent, na.rm = T)
    distTable_2[parent,]$dist_B_to_parent = mean(parentSamplesDistTable$dist_B_to_parent, na.rm = T)  
    distTable_2[parent,]$ratio = mean(parentSamplesDistTable$ratio, na.rm = T)  
  }
  
  return(distTable_2)
  
  
  
  
  #-------------------------------------
  for (parent in Parent_A) {
    #break
    dist_A1_parent_h1 = dist[homeolog_A.1 , paste(parent, "-a1", sep = "")]
    dist_A1_parent_h2 = dist[homeolog_A.1 , paste(parent, "-a2", sep = "")]
    dist_A2_parent_h1 = dist[homeolog_A.2 , paste(parent, "-a1", sep = "")]
    dist_A2_parent_h2 = dist[homeolog_A.2 , paste(parent, "-a2", sep = "")]
    
    dist_B1_parent_h1 = dist[homeolog_B.1 , paste(parent, "-a1", sep = "")]
    dist_B1_parent_h2 = dist[homeolog_B.1 , paste(parent, "-a2", sep = "")]
    dist_B2_parent_h1 = dist[homeolog_B.2 , paste(parent, "-a1", sep = "")]
    dist_B2_parent_h2 = dist[homeolog_B.2 , paste(parent, "-a2", sep = "")]
    
    distTable[parent,]$dist_A_to_parent = mean(dist_A1_parent_h1, dist_A1_parent_h2, dist_A2_parent_h1, dist_A2_parent_h2)
    distTable[parent,]$dist_B_to_parent = mean(dist_B1_parent_h1, dist_B1_parent_h2, dist_B2_parent_h1, dist_B2_parent_h2)
    distTable[parent,]$ratio = max(distTable[parent,]$dist_A_to_parent, distTable[parent,]$dist_B_to_parent)/min(distTable[parent,]$dist_A_to_parent, distTable[parent,]$dist_B_to_parent)  
  }
  
  dt_2 = apply(distTable, 2, mean)
  dt_2 = t(dt_2)
  return(dt_2)
}

distancesToParents <- function(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, Parent_A) {
  
  # tree$tip.label = sub("_$", "", tree$tip.label)
  
  dist = cophenetic(tree)  
  distTable = data.frame(dist_A_to_parent = numeric(length = length(Parent_A)), dist_B_to_parent = numeric(length = length(Parent_A)), ratio = numeric(length = length(Parent_A)))
  rownames(distTable) = Parent_A
  
  
  for (parent in Parent_A) {
    #break
    dist_A1_parent_h1 = dist[homeolog_A.1 , paste(parent, "-a1", sep = "")]
    dist_A1_parent_h2 = dist[homeolog_A.1 , paste(parent, "-a2", sep = "")]
    dist_A2_parent_h1 = dist[homeolog_A.2 , paste(parent, "-a1", sep = "")]
    dist_A2_parent_h2 = dist[homeolog_A.2 , paste(parent, "-a2", sep = "")]
      
    dist_B1_parent_h1 = dist[homeolog_B.1 , paste(parent, "-a1", sep = "")]
    dist_B1_parent_h2 = dist[homeolog_B.1 , paste(parent, "-a2", sep = "")]
    dist_B2_parent_h1 = dist[homeolog_B.2 , paste(parent, "-a1", sep = "")]
    dist_B2_parent_h2 = dist[homeolog_B.2 , paste(parent, "-a2", sep = "")]
      
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

evaluateDistance <- function(distanceTable1, sample, between_homeolog_distance) {
  if (  (distanceTable1[[sample]]$dist[2] / distanceTable1[[sample]]$dist[1] ) >  between_homeolog_distance &
        distanceTable1[[sample]]$twoOrThree[1] ==  "2pairs" &
        distanceTable1[[sample]]$twoOrThree[2] ==  "trio"
  ) {  
    return(TRUE)
    
  } else {  # NO, it is not OK
    
    return(FALSE)
  }
  
}

evaluateDistance_sortSequence <- function(distanceTable1, sample, between_homeolog_distance, tree, Parent_A_list, use_between_parents_distance, between_parents_distance, tocat, if_below_treshold, seq, homeolog_A.1, homeolog_A.2, homeolog_B.1,  homeolog_B.2) {
      toSolve = TRUE
      if (  (distanceTable1[[sample]]$dist[2] / distanceTable1[[sample]]$dist[1] ) >  between_homeolog_distance &
            distanceTable1[[sample]]$twoOrThree[1] ==  "2pairs" &
            distanceTable1[[sample]]$twoOrThree[2] ==  "trio"
      ) {  
        
        # YES, they are
        # let's calculate distances of two pairs to a potential parent
        dt_2 = distancesToParentsList(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, Parent_A_list)
        
        #********************************************************
        if ( isTRUE(dt_2$dist_A_to_parent < dt_2$dist_B_to_parent) ) {
          
          # use_between_parents_distance?
          if (use_between_parents_distance & (dt_2$ratio < between_parents_distance)){
            tocat = c(tocat, "masked/removed", "\t")
            
            seq = do_if_below_treshold(if_below_treshold, sample, seq)
            # it doesn't matter how, but sequences have to be renamed
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
            
            toSolve = FALSE
          }
          
          if (toSolve) {
            
          colnames(dt_2) = c("A_homeolog", "B_homeolog", "ratio")
          
          # rename fasta sequence
          rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
          
          }
          
        } else if (isTRUE(dt_2$dist_A_to_parent > dt_2$dist_B_to_parent )){
          
          # use_between_parents_distance ?
          if (use_between_parents_distance & (dt_2$ratio < between_parents_distance)){
            tocat = c(tocat, "masked/removed", "\t")
            
            seq = do_if_below_treshold(if_below_treshold, sample, seq)
            
            # it doesn't matter how, but sequences have to be renamed
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
            
            toSolve = FALSE
          }
          
          if (toSolve) {
            
          colnames(dt_2) = c("B_homeolog", "A_homeolog", "ratio") 
          
          # rename fasta sequence
          rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-B2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-A2", sep = "")
          
          }
          
        } else  { 
          
          tocat = c(tocat, "masked/removed", "\t")
          
          seq = do_if_below_treshold(if_below_treshold, sample, seq)
          
          # it doesn't matter how, but sequences have to be renamed
          rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
          
          toSolve = FALSE
        }
        
        if (toSolve) {
          tocat = c(tocat, "*** conditions are met", "\t")
          
          toSolve = FALSE
        }
        
        
      } else {  # NO, it is not OK
        if (toSolve) {
          
        tocat = c(tocat, "masked/removed", "\t")
        
        seq = do_if_below_treshold(if_below_treshold, sample, seq)
        
        # it doesn't matter how, but sequences have to be renamed
        #homeolog_A.1 = paste(sample, "-A1", sep = ""); homeolog_A.2 = paste(sample, "-A2", sep = ""); homeolog_B.1 = paste(sample, "-B1", sep = "");  homeolog_B.2 = paste(sample, "-B2", sep = "")
        #homeolog_A.1 = paste(sample, "-a1", sep = ""); homeolog_A.2 = paste(sample, "-a2", sep = ""); homeolog_B.1 = paste(sample, "-a3", sep = "");  homeolog_B.2 = paste(sample, "-a4", sep = "")
        
        rownames(seq)[grep(paste(sample, "-a1", sep = ""), rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(paste(sample, "-a2", sep = ""), rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(paste(sample, "-a3", sep = ""), rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(paste(sample, "-a4", sep = ""), rownames(seq))] = paste(sample, "-B2", sep = "")
        # if samples ends with a1 ... a4, they have to be renamed. if not, they are named as A1 A2 B1 B2 from previous run.
        #rownames(seq)[grep(paste(sample, "-A1", sep = ""), rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(paste(sample, "-A2", sep = ""), rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(paste(sample, "-B1", sep = ""), rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(paste(sample, "-B2", sep = ""), rownames(seq))] = paste(sample, "-B2", sep = "")
      # rownames(seq)[grep(paste(sample, "-a1 ", sep = ""), rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(paste(sample, "-a2 ", sep = ""), rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(paste(sample, "-a3 ", sep = ""), rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(paste(sample, "-a4 ", sep = ""), rownames(seq))] = paste(sample, "-B2", sep = "")
        
        }
      }
      
      return(list(seq, tocat))
}
