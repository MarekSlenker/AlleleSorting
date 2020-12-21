# hladanieParo2_funkcie

calculate_h1..H4_Distances <- function(tree, samples) {
  
  dist = cophenetic(tree)  #teraz je to naozaj vzdialenost, tam po spolocny uzol i naspat
  #colnames(dist)
  #distTable_1 = array(dim = c(7 , 7, length(samples$V1)))
  distTable_1 = list()
  
  
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
    
    
    #distTable_1[which(sample == samples$V1)] = list(distTable_sample)
    distTable_1[sample] = list( distTable_sample[order(distTable_sample$dist),])
  }
  
  return(distTable_1)
}

calculate_AB_Distances <- function(tree, samples) {
  
  dist = cophenetic(tree)  #teraz je to naozaj vzdialenost, tam po spolocny uzol i naspat
  #distTable_1 = array(dim = c(7 , 7, length(samples$V1)))
  distTable_1 = list()
  
  
  for (sample in samples$V1) {
    if (! (paste(sample, "-A1", sep = "") %in% tree$tip.label)) {
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
    
    
    #distTable_1[which(sample == samples$V1)] = list(distTable_sample)
    distTable_1[sample] = list( distTable_sample[order(distTable_sample$dist),])
  }
  
  return(distTable_1)
}

distancesToParents <- function(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, parents) {
  
  tree$tip.label = sub("_$", "", tree$tip.label)
  
  dist = cophenetic(tree)  #teraz je to naozaj vzdialenost, tam po spolocny uzol i naspat
  distTable_2 = data.frame(dist_A_to_parent = numeric(length = length(parents)), dist_B_to_parent = numeric(length = length(parents)), ratio = numeric(length = length(parents)))
  rownames(distTable_2) = parents
  
  
  for (parent in parents) {
    
    parentSamples = read.table(parent, header = FALSE,  stringsAsFactors = FALSE)
    
    parentSamplesDistTable = data.frame(dist_A_to_parent = numeric(length = length(parentSamples$V1)), 
                                        dist_B_to_parent = numeric(length = length(parentSamples$V1)), 
                                        ratio = numeric(length = length(parentSamples$V1)))
    
    for (parentSample in parentSamples$V1) {
      dist_A1_parent_h1 = dist[homeolog_A.1 , paste(parentSample, "-h1", sep = "")]
      dist_A1_parent_h2 = dist[homeolog_A.1 , paste(parentSample, "-h2", sep = "")]
      dist_A2_parent_h1 = dist[homeolog_A.2 , paste(parentSample, "-h1", sep = "")]
      dist_A2_parent_h2 = dist[homeolog_A.2 , paste(parentSample, "-h2", sep = "")]
      
      dist_B1_parent_h1 = dist[homeolog_B.1 , paste(parentSample, "-h1", sep = "")]
      dist_B1_parent_h2 = dist[homeolog_B.1 , paste(parentSample, "-h2", sep = "")]
      dist_B2_parent_h1 = dist[homeolog_B.2 , paste(parentSample, "-h1", sep = "")]
      dist_B2_parent_h2 = dist[homeolog_B.2 , paste(parentSample, "-h2", sep = "")]
      
      parentSamplesDistTable[which(parentSample == parentSamples$V1), 1] = mean( c( dist_A1_parent_h1, dist_A1_parent_h2, dist_A2_parent_h1, dist_A2_parent_h2)) 
      parentSamplesDistTable[which(parentSample == parentSamples$V1), 2] = mean( c( dist_B1_parent_h1, dist_B1_parent_h2, dist_B2_parent_h1, dist_B2_parent_h2)) 
      parentSamplesDistTable[which(parentSample == parentSamples$V1), 3] = max(parentSamplesDistTable[which(parentSample == parentSamples$V1), 2], parentSamplesDistTable[which(parentSample == parentSamples$V1), 1])/min(parentSamplesDistTable[which(parentSample == parentSamples$V1), 2], parentSamplesDistTable[which(parentSample == parentSamples$V1), 1])    
    }
    
    distTable_2[parent,]$dist_A_to_parent = mean(parentSamplesDistTable$dist_A_to_parent)
    distTable_2[parent,]$dist_B_to_parent = mean(parentSamplesDistTable$dist_B_to_parent)  
    distTable_2[parent,]$ratio = mean(parentSamplesDistTable$ratio)  
  }
  
  return(distTable_2)
}

maskujN <- function(sample, seq) {
  
  samplePos = grep(sample, rownames(seq))
  
  x = apply(seq[samplePos,], 2, function(x) {
    length(unique(x))
  })
  
  seq[samplePos, which(x > 1)] = "n"
  
  return(seq)
}

identitiHomeologs_by_1parent_mask_by_N <- function(pattern) {
  tocat = character();  tocat = c("\n", pattern, "\t")
  
  tree = read.tree(file = list.files(treeDir, pattern = pattern, full.names = TRUE))
  seq = read.dna(file = list.files(seqDir, pattern = pattern, full.names = TRUE), format = "fasta", as.character = T)
  
  distanceTable1 = calculate_h1..H4_Distances(tree, samples)
  
  for (sample in samples$V1) {
    
    # JE OK? = splna parametre?
    if (  (distanceTable1[[sample]]$dist[2] / distanceTable1[[sample]]$dist[1] ) >  treshold_distance & 
          distanceTable1[[sample]]$twoOrThree[1] ==  "2pairs" &
          distanceTable1[[sample]]$twoOrThree[2] ==  "trio"
    ) {  # je OK
      homeolog_A.1 = paste(sample, distanceTable1[[sample]]$homeolog_A.1[1], sep = ""); homeolog_A.2 = paste(sample, distanceTable1[[sample]]$homeolog_A.2[1], sep = ""); homeolog_B.1 = paste(sample, distanceTable1[[sample]]$homeolog_B.1[1], sep = "");  homeolog_B.2 = paste(sample, distanceTable1[[sample]]$homeolog_B.2[1], sep = "")
      
      # vyrataj vzdialenosti aliel ku rodicom
      dt_2 = distancesToParents(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, parents)
      
          if (
            dt_2[parents[1],]$dist_A_to_parent < dt_2[parents[1],]$dist_B_to_parent #  & dt_2[parents[1],]$ratio > treshold_distance
          ) {
            colnames(dt_2) = c("amara_homeolog", "acris_homeolog")
        
            # zmenim v strome ; presisem seq
            # tree$tip.label[which(tree$tip.label == homeolog_A.1)] = paste(sample, "-A1", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_A.2)] = paste(sample, "-A2", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_B.1)] = paste(sample, "-B1", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_B.2)] = paste(sample, "-B2", sep = "")
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
        
          } else if (
            dt_2[parents[1],]$dist_A_to_parent > dt_2[parents[1],]$dist_B_to_parent # & dt_2[parents[1],]$ratio > treshold_distance
          ){
            colnames(dt_2) = c("acris_homeolog", "amara_homeolog") 
        
            # tree$tip.label[which(tree$tip.label == homeolog_A.1)] = paste(sample, "-B1", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_A.2)] = paste(sample, "-B2", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_B.1)] = paste(sample, "-A1", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_B.2)] = paste(sample, "-A2", sep = "")
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-B2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-A2", sep = "")
        
          } else  { # vyhodim
            
            tocat = c(tocat, "ratio_ku_rodicovi_je_rovnake_pre_AaB alebo nesplna treshold_distance = Maskujem_N", "\t")
            
            seq = maskujN(sample, seq)
            # musis ich premenovat, je jedno ako, lebo sekvencie su rovnake
            rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
            next
          }
      
      tocat = c(tocat, "** zapisem distances", "\t")
      #sample_distance_table[[sample]] = cbind(sample_distance_table[[sample]], dt_2)
      
    } else {  # nie je OK
      tocat = c(tocat, "___mesplna treshold_distance_maskujem_N", "\t")
      
      seq = maskujN(sample, seq)
      # musis ich premenovat, je jedno ako, lebo sekvencie su rovnake
      homeolog_A.1 = paste(sample, "-A1", sep = ""); homeolog_A.2 = paste(sample, "-A2", sep = ""); homeolog_B.1 = paste(sample, "-B1", sep = "");  homeolog_B.2 = paste(sample, "-B2", sep = "")
      rownames(seq)[grep(paste(sample, "-h1 ", sep = ""), rownames(seq))] = paste(sample, "-A1", sep = ""); 
      rownames(seq)[grep(paste(sample, "-h2 ", sep = ""), rownames(seq))] = paste(sample, "-A2", sep = "");
      rownames(seq)[grep(paste(sample, "-h3 ", sep = ""), rownames(seq))] = paste(sample, "-B1", sep = "");
      rownames(seq)[grep(paste(sample, "-h4 ", sep = ""), rownames(seq))] = paste(sample, "-B2", sep = "")
      
      
      # vyhodim
      # tree = drop.tip(tree, tip = grep(tree$tip.label, pattern = sample, value = T))
      #seq = seq[- grep(sample, rownames(seq)),]
    }
    
    
  }
  
  #tree = drop.tip(tree, tip = c(
  # "barC007_103-h1", "barC007_103-h2", "barC007_103-h3", "barC007_103-h4",
  #"barC010_102-h1", "barC010_102-h2", "barC010_102-h3", "barC010_102-h4",
  #"barC011_101-h1", "barC011_101-h2", "barC011_101-h3", "barC011_101-h4"
  #"ampC041_101-h1", "ampC041_101-h2", "ampC041_101-h3", "ampC041_101-h4"
  #))
  #seq = seq[- grep(paste(toGrep, collapse = "|") , rownames(seq)),]
  
  # write.tree(tree, file = paste(resdir, "/", pattern, ".AB_homeologs.support", sep = ""))
  write.dna(seq, file = paste(resdir, "/", pattern, ".AB_homeologs.fasta", sep = ""), format = "fasta", nbcol = -1, colsep = "")
  
  lck = lock(logfile, exclusive = TRUE, timeout = Inf)
  cat(tocat, file = logfile, append = T)
  unlock(lck)
}

# pouzijes uz na rozdelene stromy na A  a  B
identitiHomeologs_by_1parent_make_dist_table <- function(pattern) {
  distCat = character();  distCat = c("\n", pattern, "\t")
  
  tree = read.tree(file = list.files(treeDir, pattern = pattern, full.names = TRUE))
  tree$tip.label = sub("_$", "", tree$tip.label, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  
  
  
  #############  GRRRRRRRRRRRRRRRRRRRR tu je to strasne skaredo spravene
  barbareoidesInTree = unique(sub("-.*", "", grep("bar", tree$tip.label,value = T)))
  
  
  for (sample in samples$V1) {
    if (! (sample %in% barbareoidesInTree)) {
      distCat = c(distCat, paste(rep("NA", length(parents)), collapse = "\t" ), "\t", paste(rep("NA", length(parents)), collapse = "\t" ), "\t")
      next
    }
      # JE OK? = splna parametre?   - ano, to si riesil v skorsich krokoch
      
      homeolog_A.1 = paste(sample, "-A1", sep = ""); homeolog_A.2 = paste(sample, "-A2", sep = ""); homeolog_B.1 = paste(sample, "-B1", sep = "");  homeolog_B.2 = paste(sample, "-B2", sep = "")
      
      # vyrataj vzdialenosti aliel ku rodicom
      dt_2 = distancesToParents(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, parents)
        
      distCat = c(distCat, paste(dt_2$dist_A_to_parent, collapse = "\t" ), "\t", paste(dt_2$dist_B_to_parent, collapse = "\t" ), "\t")
      
      
    }
  
  lck2 = lock(distTableFile, exclusive = TRUE, timeout = Inf)
  cat(distCat, file = distTableFile, append = T)
  unlock(lck2)
}

identitiHomeologs_by_1parent_REMOVE_if_bellow_treshold <- function(pattern) {
  tocat = character();  tocat = c("\n", pattern, "\t")
  
  tree = read.tree(file = list.files(treeDir, pattern = pattern, full.names = TRUE))
  seq = read.dna(file = list.files(seqDir, pattern = pattern, full.names = TRUE), format = "fasta", as.character = T)
  
  distanceTable1 = calculate_h1..H4_Distances(tree, samples)
  
  for (sample in samples$V1) {
    
    # JE OK? = splna parametre?
    if (  (distanceTable1[[sample]]$dist[2] / distanceTable1[[sample]]$dist[1] ) >  treshold_distance & 
          distanceTable1[[sample]]$twoOrThree[1] ==  "2pairs" &
          distanceTable1[[sample]]$twoOrThree[2] ==  "trio"
    ) {  # je OK
      homeolog_A.1 = paste(sample, distanceTable1[[sample]]$homeolog_A.1[1], sep = ""); homeolog_A.2 = paste(sample, distanceTable1[[sample]]$homeolog_A.2[1], sep = ""); homeolog_B.1 = paste(sample, distanceTable1[[sample]]$homeolog_B.1[1], sep = "");  homeolog_B.2 = paste(sample, distanceTable1[[sample]]$homeolog_B.2[1], sep = "")
      
      # vyrataj vzdialenosti aliel ku rodicom
      dt_2 = distancesToParents(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, parents)
      
      if (
        dt_2[parents[1],]$dist_A_to_parent < dt_2[parents[1],]$dist_B_to_parent #    & dt_2[parents[1],]$ratio > treshold_distance
      ) {
        colnames(dt_2) = c("amara_homeolog", "acris_homeolog")
        
        # zmenim v strome ; presisem seq
        # tree$tip.label[which(tree$tip.label == homeolog_A.1)] = paste(sample, "-A1", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_A.2)] = paste(sample, "-A2", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_B.1)] = paste(sample, "-B1", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_B.2)] = paste(sample, "-B2", sep = "")
        rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
        
      } else if (
        dt_2[parents[1],]$dist_A_to_parent > dt_2[parents[1],]$dist_B_to_parent  # & dt_2[parents[1],]$ratio > treshold_distance
      ){
        colnames(dt_2) = c("acris_homeolog", "amara_homeolog") 
        
        # tree$tip.label[which(tree$tip.label == homeolog_A.1)] = paste(sample, "-B1", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_A.2)] = paste(sample, "-B2", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_B.1)] = paste(sample, "-A1", sep = ""); tree$tip.label[which(tree$tip.label == homeolog_B.2)] = paste(sample, "-A2", sep = "")
        rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-B2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-A2", sep = "")
        
      } else  { # vyhodim
        
        tocat = c(tocat, "ratio_ku_rodicovi_je_rovnake_pre_AaB_Maskujem_N", "\t")
        
        seq = maskujN(sample, seq)
        # musis ich premenovat, je jedno ako, lebo sekvencie su rovnake
        rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = ""); rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = ""); rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = ""); rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
        next
      }
      
      tocat = c(tocat, "** zapisem distances", "\t")
      #sample_distance_table[[sample]] = cbind(sample_distance_table[[sample]], dt_2)
      
    } else {  # nie je OK
      tocat = c(tocat, "___mesplna treshold distance VYHADZUJEM zo seq", "\t")
      
      # vyhodim
      seq = seq[- grep(sample, rownames(seq)),]
      
      
      # tree = drop.tip(tree, tip = grep(tree$tip.label, pattern = sample, value = T))
      
    }
    
    
  }
  
  #tree = drop.tip(tree, tip = c(
  # "barC007_103-h1", "barC007_103-h2", "barC007_103-h3", "barC007_103-h4",
  #"barC010_102-h1", "barC010_102-h2", "barC010_102-h3", "barC010_102-h4",
  #"barC011_101-h1", "barC011_101-h2", "barC011_101-h3", "barC011_101-h4"
  #"ampC041_101-h1", "ampC041_101-h2", "ampC041_101-h3", "ampC041_101-h4"
  #))
  #seq = seq[- grep(paste(toGrep, collapse = "|") , rownames(seq)),]
  
  # write.tree(tree, file = paste(resdir, "/", pattern, ".AB_homeologs.support", sep = ""))
  write.dna(seq, file = paste(resdir, "/", pattern, ".AB_homeologs.fasta", sep = ""), format = "fasta", nbcol = -1, colsep = "")
  
  lck = lock(logfile, exclusive = TRUE, timeout = Inf)
  cat(tocat, file = logfile, append = T)
  unlock(lck)
  
}


chceckHomeologDistancesTo_1_Parent <- function(treeID) {
    
    pattern = strsplit(treeID, split = ".", fixed = T)[[1]][1]
    tocat = character();  tocat = c("\n", pattern, "\t")
    tree = read.tree(file = paste(treeDir, treeID, sep = "/"))
    seq = read.dna(file = list.files(seqDir, pattern = pattern, full.names = TRUE), format = "fasta", as.character = T)
    
      distanceTable1 = calculate_AB_Distances(tree, samples)
    
    for (sample in samples$V1) {
      
      # JE OK? = splna parametre?
      if (  (distanceTable1[[sample]]$dist[2] / distanceTable1[[sample]]$dist[1] ) >  treshold_distance & 
            distanceTable1[[sample]]$twoOrThree[1] ==  "2pairs" &
            distanceTable1[[sample]]$twoOrThree[2] ==  "trio"
      ) {  # je OK
        homeolog_A.1 = paste(sample, distanceTable1[[sample]]$homeolog_A.1[1], sep = ""); homeolog_A.2 = paste(sample, distanceTable1[[sample]]$homeolog_A.2[1], sep = ""); homeolog_B.1 = paste(sample, distanceTable1[[sample]]$homeolog_B.1[1], sep = "");  homeolog_B.2 = paste(sample, distanceTable1[[sample]]$homeolog_B.2[1], sep = "")
        
        # vyrataj vzdialenosti aliel ku rodicom
        dt_2 = distancesToParents(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, parents)
        
        if (
          dt_2[parents[1],]$ratio < treshold_distance
        ) {
          tocat = c(tocat, "ratio_ku_rodicom_je_mensi_ako_tresholdDistance_VYHADZUJEM_SEQ", "\t")
          
          seq = seq[- grep(sample, rownames(seq)),]
          
          # tree = drop.tip(tree, tip = grep(tree$tip.label, pattern = sample, value = T))
          
        } else  { # vyhodim
          tocat = c(tocat, "distance>ratio_OK", "\t")
        }
        
        
      } else {  # nie je OK
        tocat = c(tocat, "___nie_je_OK_VYHADZUJEM_SEQ", "\t")
        
        seq = seq[- grep(sample, rownames(seq)),]
        
        # tree = drop.tip(tree, tip = grep(tree$tip.label, pattern = sample, value = T))
        
      }
      
      
      
    }
    
      # write.tree(tree, file = paste(resdir, "/", pattern, ".AB_homeologs.support", sep = ""))
    write.dna(seq, file = paste(resdir, "/", pattern, ".", prefix, ".fasta", sep = ""), format = "fasta", nbcol = -1, colsep = "")
    
    lck = lock(logfile, exclusive = TRUE, timeout = Inf)
    cat(tocat, file = logfile, append = T)
    unlock(lck)
  }









identitiHomeologs_2parents <- function(pattern) {
  
  tocat = character()
  tocat = c("\n", pattern, "\t")
  
  tree = read.tree(file = list.files(treeDir, pattern = pattern, full.names = TRUE))
  seq = read.dna(file = list.files(seqDir, pattern = pattern, full.names = TRUE), format = "fasta", as.character = T)
  
  distanceTable1 = calculateDistances(tree, samples)
  
  for (sample in samples$V1) {
    
    # JE OK? = splna parametre?
    if (  (distanceTable1[[sample]]$dist[2] / distanceTable1[[sample]]$dist[1] ) >  treshold_distance & 
          distanceTable1[[sample]]$twoOrThree[1] ==  "2pairs" &
          distanceTable1[[sample]]$twoOrThree[2] ==  "trio"
    ) {  # je OK
      
      # vyrataj vzdialenosti aliel ku rodicom
      homeolog_A.1 = paste(sample, distanceTable1[[sample]]$homeolog_A.1[1], sep = "")
      homeolog_A.2 = paste(sample, distanceTable1[[sample]]$homeolog_A.2[1], sep = "")
      homeolog_B.1 = paste(sample, distanceTable1[[sample]]$homeolog_B.1[1], sep = "")
      homeolog_B.2 = paste(sample, distanceTable1[[sample]]$homeolog_B.2[1], sep = "")
      
      dt_2 = distancesToParents(tree, homeolog_A.1, homeolog_A.2, homeolog_B.1, homeolog_B.2, parents)
      
      if (
        dt_2[parents[1],]$dist_A_to_parent < dt_2[parents[1],]$dist_B_to_parent 
        &
        dt_2[parents[1],]$ratio > treshold_distance
        &
        dt_2[parents[2],]$dist_A_to_parent > dt_2[parents[2],]$dist_B_to_parent
        &
        dt_2[parents[2],]$ratio > treshold_distance  
      ) {
        colnames(dt_2) = c("amara_homeolog", "acris_homeolog")
        
        # zmenim v strome
        tree$tip.label[which(tree$tip.label == homeolog_A.1)] = paste(sample, "-A1", sep = "")
        tree$tip.label[which(tree$tip.label == homeolog_A.2)] = paste(sample, "-A2", sep = "")
        tree$tip.label[which(tree$tip.label == homeolog_B.1)] = paste(sample, "-B1", sep = "")
        tree$tip.label[which(tree$tip.label == homeolog_B.2)] = paste(sample, "-B2", sep = "")
        # presisem seq
        rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-A1", sep = "")
        rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-A2", sep = "")
        rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-B1", sep = "")
        rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-B2", sep = "")
        
      } else if (
        dt_2[parents[1],]$dist_A_to_parent > dt_2[parents[1],]$dist_B_to_parent 
        &
        dt_2[parents[1],]$ratio > treshold_distance
        &
        dt_2[parents[2],]$dist_A_to_parent < dt_2[parents[2],]$dist_B_to_parent
        &
        dt_2[parents[2],]$ratio > treshold_distance
      ){
        colnames(dt_2) = c("acris_homeolog", "amara_homeolog") 
        
        # zmenim v strome
        tree$tip.label[which(tree$tip.label == homeolog_A.1)] = paste(sample, "-B1", sep = "")
        tree$tip.label[which(tree$tip.label == homeolog_A.2)] = paste(sample, "-B2", sep = "")
        tree$tip.label[which(tree$tip.label == homeolog_B.1)] = paste(sample, "-A1", sep = "")
        tree$tip.label[which(tree$tip.label == homeolog_B.2)] = paste(sample, "-A2", sep = "")
        # presisem seq
        rownames(seq)[grep(homeolog_A.1, rownames(seq))] = paste(sample, "-B1", sep = "")
        rownames(seq)[grep(homeolog_A.2, rownames(seq))] = paste(sample, "-B2", sep = "")
        rownames(seq)[grep(homeolog_B.1, rownames(seq))] = paste(sample, "-A1", sep = "")
        rownames(seq)[grep(homeolog_B.2, rownames(seq))] = paste(sample, "-A2", sep = "")
        
        
      } else  {
        tocat = c(tocat, "'--nesplna_amara>acris", "\t")
        
        # vyhodim zo stromu
        tree = drop.tip(tree, tip = grep(tree$tip.label, pattern = sample, value = T))
        
        # vyhodim zo seq
        seq = seq[- grep(sample, rownames(seq)),]
        
        next
      }
      
      ######
      #dt_2$ratio = numeric(length = length(parents))
      #for (parent in parents) {
      #  dt_2[parent,]$ratio = max(dt_2[parent,]$amara_homeolog, dt_2[parent,]$acris_homeolog) /min(dt_2[parent,]$amara_homeolog, dt_2[parent,]$acris_homeolog)
      #}
      ######
      tocat = c(tocat, "** zapisem distances", "\t")
      #sample_distance_table[[sample]] = cbind(sample_distance_table[[sample]], dt_2)
      
    } else {
      tocat = c(tocat, "___mesplna treshold_distance", "\t")
      
      #vyhodim zo stromu
      tree = drop.tip(tree, tip = grep(tree$tip.label, pattern = sample, value = T))
      # vyhodim zo seq
      seq = seq[- grep(sample, rownames(seq)),]
    }
    
    
  }
  
  #tree = drop.tip(tree, tip = c(
  # "barC007_103-h1", "barC007_103-h2", "barC007_103-h3", "barC007_103-h4",
  #"barC010_102-h1", "barC010_102-h2", "barC010_102-h3", "barC010_102-h4",
  #"barC011_101-h1", "barC011_101-h2", "barC011_101-h3", "barC011_101-h4"
  #"ampC041_101-h1", "ampC041_101-h2", "ampC041_101-h3", "ampC041_101-h4"
  #))
  
  #seq = seq[- grep(paste(toGrep, collapse = "|") , rownames(seq)),]
  
  write.tree(tree, file = paste(resdir, "/", pattern, ".AB_homeologs.support", sep = ""))
  write.dna(seq, file = paste(resdir, "/", pattern, ".AB_homeologs.fasta", sep = ""), format = "fasta", nbcol = -1, colsep = "")
  
  lck = lock(logfile, exclusive = TRUE, timeout = Inf)
  cat(tocat, file = logfile, append = T)
  unlock(lck)
}
