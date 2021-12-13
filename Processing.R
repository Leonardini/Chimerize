library(Biostrings) ### for the lcprefix and lcsuffix functions
library(tools)      ### for the toTitleCase function
library(stringi)    ### for the stri_reverse function
library(binhf)      ### for the shift function
library(dada2)      ### for the nwalign function
library(igraph)     ### for the topological sort, other functions
library(tidyverse)  ### for the standard tidyverse functions
library(magrittr)   ### for the assignment pipe and aliases
options(warn = 0)   ### avoid hard fails due to warnings

CHARS_PER_LINE = 60L### parameter for the GBF file format
NUM_CHIMERAS = 10L  ### number of chimeras to construct
ALPHABET = c("a","c","g","t")
CLEAN_ONLY = TRUE   ### determines whether only genomes containing letters from the alphabet are admissible
ANI_THRESHOLD = 0.9       ### threshold ANI similarity applied to genes, above which they are considered similar
ANI_THRESHOLD_INTER = 0.8 ### threshold ANI similarity applied to intergenic regions, above which they are considered similar
MAX_DIFF_LEN = 0.1        ### threshold applied to genes, below which they are considered to have a similar length
MAX_DIFF_LEN_INTER = 0.2  ### threshold applied to intergenic regions, below which they are considered to have a similar length
READ_LENGTH = 150L                ### average length of a single read
HALF_READ_LENGTH = 75L            ### used to leave a bit on each end (hard-coded to avoid non-integer results of the division)
INSERT_LENGTH = 3L * READ_LENGTH  ### minimum length of a segment (junction or inter-junctional region)

ProjectDir = "~"                                                    ### the folder with the bacterial sequences where we start
subDir1 = "Salmonella enterica"                                     ### the source of the first part of each chimera
subDir2 = c("Escherichia coli", "Shigella sonnei")                  ### the source of the second part of each chimera
curDir = getwd()
setwd(ProjectDir)

### This function left-rotates a string using reversals (https://www.geeksforgeeks.org/left-rotation-right-rotation-string-2/)
### If you want to right-rotate a string, multiply numChar by -1. The default value of numChar, 0, does not change the string.
shiftString = function(inputString, numChar = 0) {
  numChar %<>%
    mod(str_length(inputString))
  if (numChar > 0) {
    str_sub(inputString, end = numChar) %<>%
      stri_reverse
    str_sub(inputString, start = numChar + 1) %<>%
      stri_reverse
    inputString %<>%
      stri_reverse
  }
  inputString
}

### Modification of the function from https://github.com/WGS-TB/Borrelia/blob/master/Utilities.R to work directly on strings
reverseComplement = function(inputString) {
  if (str_length(inputString) == 0) { ### important corner case: the reverse complement of the empty string is the empty string
    return(inputString)
  }
  revCompInds <- inputString %<>%
    str_split("") %>%
    extract2(1) %>%
    match(ALPHABET) %>%
    rev %>%
    multiply_by(-1) %>%
    add(length(ALPHABET) + 1)
  revCompString <- ALPHABET[revCompInds] %>%
    str_c(collapse = "")
  revCompString
}

### This function extends the previous one to vectors or lists
reverseComplementAll = function(listOfStrings) {
  L <- length(listOfStrings)
  for (ind in seq_len(L)) {
    listOfStrings[[ind]] %<>% 
      reverseComplement
  }
  listOfStrings
}

### Function selects numPairs non-overlapping entries in a given distMatrix with maximum (minimum if minimize = TRUE) total.
### If greedy = TRUE, uses a greedy algorithm; otherwise, uses the Hungarian algorithm (currently not implemented!)
selectDistantPairs = function(distMatrix, numPairs, greedy = TRUE, minimize = FALSE) {
  stopifnot(numPairs <= min(dim(distMatrix)))
  usableRows <- rownames(distMatrix)
  usableCols <- colnames(distMatrix)
  selectedPairs <- matrix(NA, numPairs, 3, dimnames = list(1:numPairs, c("First", "Second", "Distance")))
  optFunction <- ifelse(minimize, min, max)
  if (greedy) {
    for (index in 1:numPairs) {
      redMatrix <- distMatrix[usableRows, usableCols, drop = TRUE]
      curOpt <- optFunction(redMatrix)
      curInds <- which(redMatrix == curOpt, arr.ind = TRUE)[1, ]
      selectedPairs[index, 1] <- usableRows[curInds[1]]
      selectedPairs[index, 2] <- usableCols[curInds[2]]
      selectedPairs[index, 3] <- curOpt
      usableRows %<>% 
        setdiff(selectedPairs[index, 1])
      usableCols %<>% 
        setdiff(selectedPairs[index, 2])
    }
  }
  selectedPairs %<>%
    as_tibble()
  selectedPairs
}

### This function removes any duplicate genes by name from a genome's annotation
removeDuplicateGenes = function(genome) {
  goodRows <- !duplicated(genome[[2]]$geneName)
  genome[[2]] %<>%
    extract(goodRows, , drop = TRUE)
  genome
}

### This function extracts information from the GBFF files in an input directory, possibly saving the results as an RData file
extractAllFiles = function(myDir = subDir1, saveFile = paste0(str_remove_all(toTitleCase(myDir), " "), ".RData")) {
  if (!is.null(saveFile) && saveFile %in% list.files()) {
    load(saveFile)
    return(fullList)
  }
  initDir <- getwd()
  setwd(myDir)
  allDirs <- list.dirs()
  allDirs %<>%
    extract(str_length(.) > 1)
  numDirs <- length(allDirs)
  print(paste("There are", numDirs, "files to extract"))
  fullList <- vector("list", numDirs)
  names(fullList) <- allDirs
  for (curDir in allDirs) {
    print(curDir)
    setwd(curDir)
    curFiles <- list.files()
    goodFile <- curFiles[str_ends(curFiles, '.gbff')]
    curGenomeAndGenes <- extractGenomeAndGenes(goodFile)
    fullList[[curDir]] <- curGenomeAndGenes
    setwd('..')
  }
  setwd(initDir)
  if (!is.null(saveFile) && !(saveFile %in% list.files())) {
    save(fullList, file = saveFile)
  }
  fullList
}

### This function computes all pairwise similarities (numbers of common genes) between two groups of gene-annotated genomes
computeAllDists = function(Results1, Results2) {
  M <- length(Results1)
  N <- length(Results2)
  distMatrix <- matrix(0, M, N, dimnames = list(names(Results1), names(Results2)))
  for (ind1 in 1:M) {
    curRes1 <- Results1[[ind1]]
    curNames1 <- curRes1[[2]]$geneName
    for (ind2 in 1:N) {
      curRes2 <- Results2[[ind2]]
      curNames2 <- curRes2[[2]]$geneName
      curSim <- length(intersect(curNames1, curNames2))
      distMatrix[ind1, ind2] <- curSim
    }
  }
  distMatrix
}

### This function extracts the genome and a representation of its genes from a gbff (GenBank Flat File) format
### It only takes the chromosome, not the plasmids, into account. If keepDuplicates = FALSE, it removes any duplicated genes.
extractGenomeAndGenes = function(fname, keepDuplicates = FALSE) {
  Lines <- fname %>%
    readLines %>%
    str_trim(side = "left")
  startLine <- which(str_starts(Lines, 'CONTIG'))[1]
  contigLength <- Lines %>%
    extract(startLine) %>%
    str_extract("[0-9]+\\)") %>%
    str_remove('\\)') %>%
    parse_integer
  expectNumLines <- ceiling(contigLength/CHARS_PER_LINE)
  stopifnot(str_length(Lines[startLine + expectNumLines + 2]) < 3)
  subLines <- Lines[startLine + 1 + (1:expectNumLines)]
  fullGenome <- subLines %>%
    str_remove("^[0-9]+") %>%
    str_remove_all(" ") %>%
    str_c(collapse = "")
  stopifnot(str_length(fullGenome) == contigLength) ### sanity check
  preContig <- Lines %>%
    extract(1:(startLine - 1))
  genePos <- preContig %>%
    str_detect('^gene ')
  geneLines <- preContig %>%
    extract(genePos) %>%
    str_remove('^gene ') %>%
    str_trim(side = "left")
  reverseComp <- geneLines %>%
    str_detect('^complement')
  geneLines %<>%
    str_remove('^complement\\(') %>%
    str_remove('\\)')
  startPos <- geneLines %>%
    str_extract("^[0-9]+") %>%
    parse_integer                          
  endPos <- geneLines %>%
    str_extract("[0-9]+$") %>%
    parse_integer
  badPos <- (is.na(startPos) | is.na(endPos))
  geneNameLines <- preContig %>%
    extract(which(genePos) + 1)
  geneName <- geneNameLines %>%
    str_remove('^/gene=\"') %>%
    str_remove('^/locus_tag=\"') %>%
    str_remove('\"$')
  allGenes <- tibble(startPos, endPos, reverseComp, geneName) %>% 
    filter(!badPos)
  output <- list(genome = fullGenome, genes = allGenes)
  if (!keepDuplicates) {
    output %<>%
      removeDuplicateGenes
  }
  output
}

### This function is a wrapper around nwalign which computes an approximate value of the average nucleotide identity (ANI)
### If the ratio of the lengths exceeds 1 + lengthMargin (an input parameter), the ANI is not computed and is assumed to be 0
### Note: we assume that two empty strings have an ANI of 1 and an empty string has an ANI of 0 with any non-empty string
computeANI <- function(Seq1, Seq2, lengthMargin) {
  Range <- sort(c(str_length(Seq1), str_length(Seq2)))
  minLen <- Range[1]
  maxLen <- Range[2]
  if (minLen == 0) {
    answer = ifelse(maxLen == 0, 1, 0)
    return(answer)
  }
  if (maxLen/minLen > 1 + lengthMargin) {
    return(0)
  }
  myAlignment <- nwalign(toupper(Seq1), toupper(Seq2), vec = (maxLen < 1000))
  numMatches <- sum(str_split(myAlignment[1], "")[[1]] == str_split(myAlignment[2], "")[[1]])
  output <- numMatches/minLen
  output
}

### This function applies the previous function to two lists of sequences
computePairsANI = function(Seqs1, Seqs2, lengthMargin = MAX_DIFF_LEN) {
  N <- length(Seqs1)
  stopifnot(length(Seqs2) == N)
  output <- rep(0, N)
  for (ind in seq_len(N)) {
    output[ind] <- computeANI(Seqs1[[ind]], Seqs2[[ind]], lengthMargin = lengthMargin)
  }
  output
}

### This function reverse complements the circular chromosome represented by the input genome
reverseComplementGenome = function(genome) {
  genomeLength <- str_length(genome[[1]])
  genome[[1]] %<>%
    reverseComplement
  genome[[2]] %<>%
    mutate(initStartPos = startPos, initEndPos = endPos) %>%
    mutate(startPos =  genomeLength - initEndPos + 1L, endPos =  genomeLength - initStartPos + 1L) %>%
    mutate_at("reverseComp", `n'est pas`) %>%
    dplyr::select(-starts_with('init')) %>%
    arrange(startPos)
  genome
}

### This function shifts the circular chromosome represented by the genome input by the specified amount.
circularShift = function(genome, shiftBy = 0L) {
  genomeLength <- str_length(genome[[1]])
  if (shiftBy != 0) {
    genome[[1]] %<>%
      shiftString(numChar = shiftBy)
    genome[[2]] %<>% ### Note that occasionally, the interval containing a gene may "wrap around" the genome after a rotation!
      mutate(startPos =  dplyr::recode(mod(startPos - shiftBy, genomeLength), `0` = genomeLength)) %>%
      mutate(endPos   =  dplyr::recode(mod(endPos   - shiftBy, genomeLength), `0` = genomeLength)) %>%
      arrange(startPos)
  }
  genome
}

### A version of str_sub that works on circular strings; it is vectorized over start and end, but is NOT vectorized over string!
str_sub_circular = function(string, start, end) {
  L0 <- length(string)
  L1 <- length(start)
  L2 <- length(end)
  L <- max(L0, L1, L2)
  stopifnot(mod(L, L0) == 0 && mod(L, L1) == 0 && mod(L, L2) == 0)
  start <- rep(start, L1/L)
  end   <- rep(end,   L2/L)
  output <- rep("", L)
  regPos <- (start <= end)
  output[regPos] <- str_sub(string, start[regPos], end[regPos])
  invPos <- (start > end)
  output[invPos] <- str_c(str_sub(string, start = start[invPos]), str_sub(string, end = end[invPos]), collapse = "")
  output
}

### This auxiliary function extracts the gene sequence for each gene in an annatation tab, using the startPos, endPos columns
completeWithSequence = function(annotationTab, genomeSeq) {
  annotationTab %<>%
    mutate(seq = str_sub_circular(genomeSeq, startPos, endPos)) %>% 
    mutate(seq = replace(seq, reverseComp, reverseComplementAll(seq[reverseComp])))
  annotationTab
}

### This function peforms an all-pair alingment of the genes between specified boundaries (gene indices) for two input genomes 
pairAlignBlock = function(genome1, from1, to1, genome2, from2, to2) {
  genes1 <- genome1[[2]] %>%
    slice(from1:to1) %>%
    completeWithSequence(genome1[[1]])
  size1 <- nrow(genes1)
  genes2 <- genome2[[2]] %>%
    slice(from2:to2) %>%
    completeWithSequence(genome2[[1]])
  size2 <- nrow(genes2)
  geneSeqs1 <- genes1 %>%
    pull(seq)
  geneSeqs2 <- genes2 %>%
    pull(seq)
  myMat <- matrix(0, size1, size2, dimnames = list(genes1$geneName, genes2$geneName))
  print(paste("There are", size1, "genes from genome 1 and", size2, "genes from genome 2 to process"))
  for (ind1 in 1:size1) {
    print(ind1)
    geneSeq1 <- rep(geneSeqs1[[ind1]], size2)
    myMat[ind1, ] <- computePairsANI(geneSeq1, geneSeqs2, lengthMargin = MAX_DIFF_LEN)
  }
  myMat
}

### This function finds the maximum number of junctions that are non-intersecting (where (i,j) <> (k,l) iff (i-k) * (j-l) > 0)
### It returns a 2-column matrix with gene names (coming from the dimnames of alignmentMatrix) corresponding to the junctions
findJunctions = function(alignmentMatrix, threshold = ANI_THRESHOLD) {
  allGoodPos <- which(alignmentMatrix >= threshold, arr.ind = TRUE) %>%
    as_tibble %>%
    arrange(row, col)
  if (nrow(allGoodPos) > 0) {
    allJunctions <- allGoodPos %>%
      slice(findMaximumCompatiblePairs(allGoodPos))
  } else {
    allJunctions <- allGoodPos
  }
  allJunctions %<>%
    mutate(gene1 = rownames(alignmentMatrix)[row], gene2 = colnames(alignmentMatrix)[col]) %>%
    dplyr::select(-row, -col)
  allJunctions
}

### This function returns the length of the longest affix (prefix if pre = TRUE, suffix otherwise) for string pairs in 2 lists 
lcaffixAll = function(strList1, strList2, pre = TRUE) {
  L <- length(strList1)
  stopifnot(length(strList2) == L)
  output <- rep(NA, L)
  for (ind in seq_len(L)) {
    output[ind] <- ifelse(pre, lcprefix(strList1[[ind]], strList2[[ind]]), lcsuffix(strList1[[ind]], strList2[[ind]]))
  }
  output
}

### This function extends the junctions as much as possible to the left and to the right
finishChimera = function(genome1, genome2, allJunctions) {
  L <- nrow(allJunctions)
  seq1 <- genome1[[1]]
  L1 <- str_length(seq1)
  seq2 <- genome2[[1]]
  L2 <- str_length(seq2)
  genes1 <- genome1[[2]]
  genes2 <- genome2[[2]]
  junctionIntervals1 <- genes1 %>%
    filter(geneName %in% allJunctions$gene1) %>%
    rowid_to_column
  junctionIntervals2 <- genes2 %>%
    filter(geneName %in% allJunctions$gene2) %>%
    rowid_to_column
  junctionIntervals <- inner_join(junctionIntervals1, junctionIntervals2, by = 'rowid') %>%
    dplyr::select(contains('Pos'))
  junctionIntervals %<>%
    mutate(addPre  = lcaffixAll(str_sub(seq1, 1, startPos.x - 1), str_sub(seq2, 1, startPos.y - 1), pre = FALSE )) %>%
    mutate(addPost = lcaffixAll(str_sub(seq1, endPos.x + 1, L1 ), str_sub(seq2, endPos.y + 1, L2 ), pre = TRUE)) %>%
    mutate(startPos.x = startPos.x - addPre , startPos.y = startPos.y - addPre) %>%
    mutate(endPos.x   = endPos.x   + addPost, endPos.y   = endPos.y   + addPost)
  output = list(seq1, seq2, junctionIntervals)
  output
}

### This function finds additional junctions in genes located in the intervals between previously matched genes, initJunctions
completeChimera = function(genome1, genome2, initJunctions) {
  genes1 <- genome1[[2]]
  genes2 <- genome2[[2]]
  numGenes1 <- nrow(genes1)
  numGenes2 <- nrow(genes2)
  L <- length(initJunctions)
  allJunctions <- c()
  bounds1 <- c(which(genes1$geneName %in% initJunctions), numGenes1 + 1)
  bounds2 <- c(which(genes2$geneName %in% initJunctions), numGenes2 + 1)
  stopifnot(all(diff(bounds1) > 0)) ### sanity check: junction genes must be consistently ordered in both genomes 
  stopifnot(all(diff(bounds2) > 0))
  for (id in 1:L) {
    curJunction <- initJunctions[id]
    curAlignment <- pairAlignBlock(genome1, bounds1[id] + 1, bounds1[id + 1] - 1, genome2, bounds2[id] + 1, bounds2[id + 1] - 1)
    curJunctions <- findJunctions(curAlignment, threshold = ANI_THRESHOLD)
    allJunctions %<>%
      bind_rows(tibble(gene1 = curJunction, gene2 = curJunction)) %>%
      bind_rows(curJunctions)
  }
  fullChimera <- finishChimera(genome1, genome2, allJunctions)
  fullChimera
}

### This function from https://lists.nongnu.org/archive/html/igraph-help/2014-11/msg00043.html finds the longest path in a DAG
findLongestPath = function(g) {
  sortOrder <- topological.sort(g)
  for (node in sortOrder) {
    curPredecessors <- V(g)[nei(node, mode = "in")]
    if (length(curPredecessors) > 0) {
      dists <- curPredecessors$dist + 1
      V(g)[node]$dist <- max(dists)
      V(g)[node]$parent <- curPredecessors[which.max(dists)]
    } else {
      V(g)[node]$dist <- 0
      V(g)[node]$parent <- 0
    }
  }
  longestPath <- c()
  curNode <- which.max(V(g)$dist)
  while (curNode != 0) {
    longestPath %<>%
      c(curNode)
    curNode <- V(g)[curNode]$parent
  }
  longestPath %<>%
    rev
  longestPath
}

### Uses a graph algorithm to find the maximum chain of rows in a 2-column matrix with respect to strict lexicographic order.
### In other words, if (i,j) and (k,l) are both chosen, then it must either be that i < k and j < l, or that i > k and j > l.
### Assumes that myMatrix is a tibble with 2 columns, "row" and "col", which is sorted by row first, and ties broken by col.
findMaximumCompatiblePairs = function(myMatrix) {
  M <- nrow(myMatrix)
  if (M == 0) {
    return()
  }
  allEdges <- matrix(NA, 0, 2)
  for (ind in 1:M) {
    curRow <- myMatrix[ind,]
    goodRows <- which(myMatrix$row > curRow$row & myMatrix$col > curRow$col)
    if (length(goodRows) > 0) {
      allEdges %<>%
        rbind(cbind(ind, goodRows))
    }
  }
  if (nrow(allEdges) == 0) {
    return(1)
  }
  G <- graph(edges = t(allEdges))
  output <- findLongestPath(G)
  output
}

### This function constructs a chimera out of a pair of annotated genomes, using common genes as anchors.
### The value of ANI_THRESHOLD is used as a cutoff to determine which stretches of the genomes to use as junctions.
### Note that both genomes get rotated and the second one gets reverse-complemented if necessary for the alignment.
### If skipCompletion = TRUE, only returns the adjusted genomes and the initial junctions, but does not complete the chimera.
constructChimera = function(genome1, genome2, skipCompletion = FALSE) {
  seq1 <- genome1[[1]]
  seq2 <- genome2[[1]]
  genes1 <- genome1[[2]] %>%
    completeWithSequence(seq1)
  genes2 <- genome2[[2]] %>%
    completeWithSequence(seq2)
  commonGenes <- inner_join(genes1, genes2, by = "geneName")
  agreeDirection <- (commonGenes$reverseComp.x == commonGenes$reverseComp.y)
  sameDirection <- ifelse(sum(agreeDirection) >= sum(!agreeDirection), TRUE, FALSE)
  extraFactor <- ifelse(sameDirection, 1, -1)
  commonGenesReduced <- commonGenes %>%
    filter(agreeDirection == sameDirection)
  commonGenesReduced %<>%
    mutate(ANI = computePairsANI(seq.x, seq.y)) %>%
    filter(ANI >= ANI_THRESHOLD)
  stopifnot(nrow(commonGenesReduced) > 0)
  scaledRankY <- rank(extraFactor * commonGenesReduced$startPos.y)
  goodStretches <- findMaximumCompatiblePairs(tibble(row = 1:nrow(commonGenesReduced), col = scaledRankY))
  commonGenesReduced %<>%
    slice(goodStretches)
  row1 <- commonGenesReduced %>%
    slice(1)
  genome1 %<>%
    circularShift(shiftBy = row1$startPos.x - 1L)
  if (sameDirection) {
    genome2 %<>%
      circularShift(shiftBy = row1$startPos.y - 1L)
  } else {
    genome2 %<>%
      circularShift(shiftBy = row1$endPos.y) %>%
      reverseComplementGenome
  }
  initJunctions <- commonGenesReduced %>%
    pull(geneName)
  if (skipCompletion) {
    myChimera <- list(genome1 = genome1, genome2 = genome2, initJunctions = initJunctions)
  } else {
    myChimera <- completeChimera(genome1, genome2, initJunctions)
  }
  myChimera
}

### This function is used to "fix" chimeras built by the original method, which is to say, by merging adjacent genes if possible
### If checkPreMerge = TRUE, the function first tests that the intergenic regions align "well" before merging any adjacent genes
fixChimera = function(preChimera, curChimera, checkPreMerge = TRUE) {
  seq1 <- curChimera[[1]]
  seq2 <- curChimera[[2]]
  stopifnot(preChimera[[1]][[1]] == seq1)
  stopifnot(preChimera[[2]][[1]] == seq2)
  summaryTab1 <- preChimera[[1]][[2]]
  summaryTab2 <- preChimera[[2]][[2]]
  mainTab <- curChimera[[3]] %>%
    mutate_if(is.numeric, as.integer)
  mainTab %<>%  ### adjust to remove possibly incorrect extension information
    mutate(startPos.x = startPos.x + addPre , startPos.y = startPos.y + addPre ) %>%
    mutate(endPos.x   = endPos.x   - addPost, endPos.y   = endPos.y   - addPost) %>%
    select(-addPre, -addPost)
  mainTab %<>%  ### match to the correct gene numbers based on known names
    mutate(gene.x = match(startPos.x, summaryTab1$startPos), gene.y = match(startPos.y, summaryTab2$startPos)) %>%
    rowid_to_column
  stopifnot(all(!is.na(mainTab$gene.x)))
  stopifnot(all(!is.na(mainTab$gene.y)))
  L1 <- nchar(seq1)
  L2 <- nchar(seq2)
  N1 <- nrow(summaryTab1)
  N2 <- nrow(summaryTab2)
  mainTab %<>%  ### annotate with information about genes and check for circular consecutiveness
    mutate(name.x = summaryTab1$geneName[gene.x], reverseComp.x = summaryTab1$reverseComp[gene.x]) %>%
    mutate(name.y = summaryTab2$geneName[gene.y], reverseComp.y = summaryTab2$reverseComp[gene.y]) %>%
    mutate(cons.x = mod(shift(gene.x, -1) - gene.x, N1) == 1, cons.y = mod(shift(gene.y, -1) - gene.y, N2) == 1) %>%
    mutate(cons.x = cons.x & cons.y, cons.y = cons.x)   ### only merge when the genes are consecutive in both genomes
  ### special case: if the last junction extends all the way around the circle, we shift the genomes and the positions
  if (mainTab %>% 
      pull(cons.x) %>% 
      tail(1)) {
    lastRow <- mainTab %>% 
      pull(cons.x) %>%
      `n'est pas` %>%
      which %>% 
      max %>%
      add(1)
    lastRow <- mainTab %>%
      slice(lastRow)
    lastPos.x <- lastRow %>%
      pull(startPos.x)
    lastPos.y <- lastRow %>%
      pull(startPos.y)
    seq1 %<>%
      shiftString(lastPos.x - 1L)
    seq2 %<>%
      shiftString(lastPos.y - 1L)
    mainTab %<>%
      mutate(startPos.x = mod(startPos.x - lastPos.x + 1L, L1), endPos.x = mod(endPos.x - lastPos.x + 1L, L1)) %>%
      mutate(startPos.y = mod(startPos.y - lastPos.y + 1L, L2), endPos.y = mod(endPos.y - lastPos.y + 1L, L2)) %>%
      arrange(startPos.x) %>%
      select(-rowid) %>%
      rowid_to_column
  }
  ### if needed, test for mergeability, then merge
  if (checkPreMerge) {
    xTab <- mainTab %>%
      select('rowid', ends_with('.x')) %>%
      rename_all(function(string) {str_remove(string, '.x')}) %>%
      mutate(startPos = mainTab$endPos.x + 1L, endPos = dplyr::recode(shift(mainTab$startPos.x, -1) - 1L, `0` = L1)) %>%
      filter(cons == TRUE) %>%
      completeWithSequence(seq1)
    yTab <- mainTab %>%
      select('rowid', ends_with('.y')) %>%
      rename_all(function(string) {str_remove(string, '.y')}) %>%
      mutate(startPos = mainTab$endPos.y + 1L, endPos = dplyr::recode(shift(mainTab$startPos.y, -1) - 1L, `0` = L2)) %>%
      filter(cons == TRUE) %>%
      completeWithSequence(seq2)
    xTab %<>%
      mutate(ANI = computePairsANI(seq, yTab$seq, lengthMargin = MAX_DIFF_LEN_INTER), cons.x = (ANI >= ANI_THRESHOLD_INTER))
  } else {
    xTab <- mainTab
  }
  mergeIDs <- xTab %>%
    filter(cons.x == TRUE) %>%
    pull(rowid)
  mergeLogical  <- (mainTab$rowid %in% mergeIDs)
  startBlocks   <- sort(unique(c(1, which(!lag(mergeLogical, 1)))))
  endBlocks     <- which(!mergeLogical)
  newEndPos.x   <- mainTab$endPos.x[endBlocks]
  newEndPos.y   <- mainTab$endPos.y[endBlocks]
  mainTab %<>%
    slice(startBlocks) %>%
    mutate(endPos.x = newEndPos.x, endPos.y = newEndPos.y) %>%
    select(contains('Pos'))
  mainTab %<>%
    mutate(addPre  = lcaffixAll(str_sub(seq1, end = startPos.x - 1), str_sub(seq2, end = startPos.y - 1), pre = FALSE)) %>%
    mutate(addPost = lcaffixAll(str_sub(seq1, start = endPos.x + 1), str_sub(seq2, start = endPos.y + 1), pre = TRUE)) %>%
    mutate(startPos.x = startPos.x - addPre , startPos.y = startPos.y - addPre) %>%
    mutate(endPos.x   = endPos.x   + addPost, endPos.y   = endPos.y   + addPost)
  ### special case: the 1st junction normally starts at position 1, but we want to prepend the common prefix of the two genomes
  L0 <- lcsuffix(seq1, seq2)
  if (L0 > 0) {
    mainTab$addPre[1] <- L0
    seq1 %<>%
      shiftString(-L0)
    seq2 %<>%
      shiftString(-L0)
    mainTab %<>%
      mutate(startPos.x = startPos.x + L0, endPos.x = endPos.x + L0, startPos.y = startPos.y + L0, endPos.y = endPos.y + L0)
    mainTab$startPos.x[1] <- 1
    mainTab$startPos.y[1] <- 1
  }
  output <- list(seq1, seq2, mainTab)
  output
}

### This function merges the adjacent intervals in an input interval table, using L1 and L2 as the length of the two sequences
### The columns that are being merged are specified in the last argument and assumed to always be present in the interval table
mergeAdjacentIntervals = function(intervalTable, L1, L2, mergeCols = c("endPos.x", "endPos.y", "addPost")) {
  intervalTable %<>%
    mutate(gap.x = mod(startPos.x - shift(endPos.x, 1) - 1L, L1), gap.y = mod(startPos.y - shift(endPos.y, 1) - 1L, L2)) %>%
    mutate(adjacent = (gap.x == 0 & gap.y == 0))
  mergeRows <- which(intervalTable$adjacent)
  if (length(mergeRows) > 0) { ### merge from the end towards the start
    for (ind in (length(mergeRows)):1) {
      curPos <- mergeRows[ind]
      intervalTable[curPos - 1, mergeCols] <- intervalTable[curPos, mergeCols]
    }
  }
  intervalTable %<>%
    slice(-mergeRows)
  intervalTable
}

### This function adjusts an input chimera to ensure that no segment or junction is shorter than the specified minimum length
### If forceEven is TRUE and the number of junctions is odd, the shortest one is eliminated
adjustChimera = function(curChimera, minLength = INSERT_LENGTH, forceEven = TRUE) {
  seq1 <- curChimera[[1]]
  seq2 <- curChimera[[2]]
  L1 <- nchar(seq1)
  L2 <- nchar(seq2)
  curTab <- curChimera[[3]]
  curTab %<>%
    mutate(length.x = endPos.x - startPos.x, length.y = endPos.y - startPos.y) %>%
    mutate(gap.x = mod(startPos.x - shift(endPos.x, 1) - 1L, L1), gap.y = mod(startPos.y - shift(endPos.y, 1) - 1L, L2)) %>%
    mutate(badJunction = pmin(length.x, length.y) < minLength, badGap.x = gap.x < minLength, badGap.y = gap.y < minLength)
  redTab <- curTab %>%
    filter(badJunction == FALSE)
  ### if only one of the gaps is small, remove preceding junction
  redTab %<>%
    filter(badGap.x == badGap.y) %>%
    rowid_to_column
  ### if both gaps are small, check the ANI of the resulting extended junction; if it is above threshold, merge; otherwise, cut
  miniTab <- redTab %>%
    filter(badGap.x == TRUE) %>%
    mutate(startPos.x = startPos.x - gap.x, startPos.y = startPos.y - gap.y) %>%
    mutate(eseq1 = str_sub(seq1, startPos.x, endPos.x), eseq2 = str_sub(seq2, startPos.y, endPos.y)) %>%
    mutate(ANI = computePairsANI(eseq1, eseq2), extendable = (ANI > ANI_THRESHOLD))
  killRows <- miniTab %>%
    filter(extendable == FALSE) %>%
    pull(rowid)
  redTab %<>%
    filter(!(rowid %in% killRows))
  saveRows <- miniTab %>%
    filter(extendable == TRUE) %>%
    pull(rowid)
  redTab[redTab$rowid %in% saveRows, ] <- miniTab[miniTab$rowid %in% saveRows, 1:ncol(redTab)]
  redTab %<>%
    mergeAdjacentIntervals(L1, L2)
  if (forceEven && mod(nrow(redTab), 2) == 1) {
    killRow <- redTab %>%
      mutate(totalLength = (endPos.x - startPos.x) + (endPos.y - startPos.y)) %>%
      pull(totalLength) %>%
      which.min
    redTab %<>%
      slice(-killRow)
  }
  redTab %<>%
    select(contains('Pos'), contains('add'))
  output <- list(seq1, seq2, redTab)
  output
}

### This is the driver function for building chimeric genomes; if cleanOnly is TRUE, only uses clean genomes (as defined at top)
mainChimera = function(cleanOnly = CLEAN_ONLY) {
  first <- extractAllFiles(subDir1)
  second <- extractAllFiles(subDir2[1])
  third <- extractAllFiles(subDir2[2])
  second %<>%
    c(third)
  if (cleanOnly) {
    firstGenomes <- sapply(first, function(x) {x[[1]]})
    firstChars <- sapply(ALPHABET, function(y) {str_count(firstGenomes, y)})
    firstClean <- which(rowSums(firstChars) == str_length(firstGenomes))
    first %<>%
      extract(firstClean)
    secondGenomes <- sapply(second, function(x) {x[[1]]})
    secondChars <- sapply(ALPHABET, function(y) {str_count(secondGenomes, y)})
    secondClean <- which(rowSums(secondChars) == str_length(secondGenomes))
    second %<>%
      extract(secondClean)
  }
  myDists <- computeAllDists(first, second)
  bestPairs <- selectDistantPairs(myDists, NUM_CHIMERAS, greedy = TRUE, minimize = TRUE)
  write_csv(bestPairs, path =  "GeneratingPairs.csv")
  chimeraList <- vector("list", NUM_CHIMERAS)
  for (index in 1:NUM_CHIMERAS) {
    print(index)
    curFilename <- paste0("Chimera", index, ".RData")
    updatedFilename <- str_replace(curFilename, ".RData", "Updated.RData")
    finalFilename <- str_replace(curFilename, ".RData", "Final.RData")
    curPair <- bestPairs[index,]
    curFirst <- first[[curPair[[1]]]]
    curSecond <- second[[curPair[[2]]]]
    if (!curFilename %in% list.files()) {
      curChimera <- constructChimera(curFirst, curSecond, skipCompletion = FALSE)
      save(curChimera, file = curFilename)
    } else {
      if (!updatedFilename %in% list.files()) {
        load(curFilename)
        preChimera <- constructChimera(curFirst, curSecond, skipCompletion = TRUE)
        curChimera <- fixChimera(preChimera, curChimera, checkPreMerge = TRUE)
        save(curChimera, file = updatedFilename)
      } else {
        load(updatedFilename)
        curChimera <- adjustChimera(curChimera)
        save(curChimera, file = finalFilename)
      }
    }
    chimeraList[[index]] <- curChimera
  }
  chimeraList
}