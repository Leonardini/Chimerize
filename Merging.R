source("Processing.R")
### codes for the type of sequence being used
FIRST = 1
SECOND = 2
JUNCTION = 3
MACHINE = "NextSeq500v2L150" ### the sequencing machine
MYSEED  = parse_integer(str_c(1:9, collapse = ""))
MFL     = 712L                    ### mean fragment length
SDEV    = 466L                    ### its standard deviation
REMOVE_FROM_ENDS = 1000L ### number of reads to remove from each end of the files to avoid detectability via head or tail
MAX_L = 3E7L ### maximum number of lines to keep from a large file
projectDir  = "~"
startDir    = projectDir
simGenomeDir = projectDir
artDir      = "~/Downloads/Downloaded software/art_src_MountRainier_MacOS"### the folder containing the ART simulator source code
extDir      = "ART_profiler_illumina"                                     ### the folder containing the additional profiles 

### This function extracts the number of lines and reads in each file matching the specified pattern in the specified directory
### If postprocessFilename = TRUE, also splits the name of each file into meaningful parts
getNumReads = function(useDir, pattern, postprocessFilename = FALSE) {
  curDir <- getwd()
  setwd(useDir)
  output <- system2("wc", args = c("-l", pattern), stdout = TRUE)
  outTab <- output %>%
    head(-1) %>%
    str_trim %>%
    str_split_fixed(" ", 2) %>%
    as_tibble(.name_repair = "minimal")
  colnames(outTab) <- c("nReads", "File")
  outTab %<>%
    mutate(nReads = as.integer(parse_integer(nReads)/4)) %>%
    mutate(Filename = File)
  if (postprocessFilename) {
    outTab %<>%
      separate(File, sep = "\\_", into = c("Patient1", "Patient2", "Sample", "Label", "Read", "Extension")) %>%
      unite("Patient", c(Patient1, Patient2)) %>%
      mutate(Sample = parse_integer(str_sub(Sample, start = 2)), Read = parse_integer(str_sub(Read, start = 2))) %>%
      select(-Label, -Extension)
  } else {
    N <- nrow(outTab)
    outTab %<>%
      arrange(File) %>%
      mutate(Sample = rep(1:(N/2), each = 2))
  }
  setwd(curDir)
  outTab
}

### This function changes the extensions of files in an input directory by replacing the inPattern with the outPattern
### It is designed to work with fastq/fq files (or compressed versions thereof), specifically assuming a paired-end read format
changeFileExtensions = function(inputDir, inPattern = "_1.fq.gz", outPattern = ".fq.gz", inds = 1:2, pre = "_") {
  curDir <- getwd()
  setwd(inputDir)
  allFiles <- list.files()
  goodFiles <- allFiles[str_detect(allFiles, inPattern)]
  transGoodFiles <- goodFiles
  for (ind in inds) {
    curInPattern  <- paste0(ind, inPattern )
    curOutPattern <- paste0(pre, ind, outPattern)
    transGoodFiles %<>%
      str_replace(curInPattern, curOutPattern)
  }
  stopifnot(all(goodFiles != transGoodFiles))
  file.rename(goodFiles, transGoodFiles)
  setwd(curDir)
}

### This function uses ART's Illumina profiler to create a noise model in the output directory for reads in an input directory
prepareErrorModel = function(inputDir, outputDir = artDir, ext = extDir, machine = MACHINE) {
  curDir <- getwd()
  setwd(outputDir)
  setwd(ext)
  output <- system2(command = "./art_profiler_illumina", args = c(machine, inputDir, "fastq"), stdout = TRUE)
  setwd(curDir)
  output
}

### From a chimera object defined by its index, prepare a table of positions, segments, and specified coverages for each one
### If fromFirst = TRUE, the first segment comes from the first genome; otherwise, it comes from the second genome
processChimera = function(index, ProjectDir, fromFirst = (index %% 2 == 1), coverages = 2 + (2 * fromFirst - 1) * c(-1, 1, 0), 
                          factor = 1) {
  coverages %<>%
    as.integer %>%
    multiply_by(factor)
  curDir <- getwd()
  setwd(ProjectDir)
  fname <- paste0("Chimera", index, "Final.RData")
  load(fname)
  seq1 <- curChimera[[1]]
  seq2 <- curChimera[[2]]
  curTab <- curChimera[[3]]
  L1 <- nchar(seq1)
  L2 <- nchar(seq2)
  N <- nrow(curTab)
  stopifnot(mod(N, 2) == 0)
  curTab %<>%
    mutate(useSeq = rep(c(!fromFirst, fromFirst) + 1L, N/2)) %>%
    mutate_all(as.integer) %>%
    mutate(junction = TRUE) %>%
    rowid_to_column()
  seqTab <- curTab
  seqTab %<>%
    mutate(startPos.x = mod(shift(endPos.x + 1L, 1), L1), endPos.x = mod(curTab$startPos.x - 1L, L1)) %>%
    mutate(startPos.y = mod(shift(endPos.y + 1L, 1), L2), endPos.y = mod(curTab$startPos.y - 1L, L2)) %>%
    mutate(junction = FALSE)
  fullTab <- bind_rows(curTab, seqTab) %>%
    arrange(rowid, junction)
  firstTab <- fullTab %>%
    filter(useSeq == 1) %>%
    mutate(seq = str_sub_circular(seq1, startPos.x, endPos.x))
  secondTab <- fullTab %>%
    filter(useSeq == 2) %>%
    mutate(seq = str_sub_circular(seq2, startPos.y, endPos.y))
  allSegs <- bind_rows(firstTab, secondTab) %>%
    arrange(rowid, junction) %>%
    mutate(label = junction * 3L + useSeq * (!junction)) %>%
    mutate(coverage = coverages[label]) %>%
    select(rowid, useSeq, junction, seq, label, coverage)
  setwd(curDir)
  allSegs
}

### This function separates a string into pieces of the specified size, with only the last piece being possibly smaller
### This is useful for fixed-width formatting situations such as fasta files
str_chunk = function(sequence, maxSize = 80) {
  L <- str_length(sequence)
  numFullChunks <- floor(L / maxSize)
  numLeft <- mod(L, maxSize)
  useExtra <- (numLeft > 0)
  size <- numFullChunks + useExtra
  output <- vector("character", size)
  curPos <- 0
  for (ind in seq_len(size)) {
    output[ind] <- str_sub(sequence, curPos + 1, curPos + maxSize)
    curPos %<>%
      add(maxSize)
  }
  output
}


createReference = function(seq, index, extraIndex = 0, dir = ProjectDir) {
  fname <- paste0(dir, 'Reference', index, 'I', extraIndex, '.fa')
  f <- file(fname, "w")
  firstLine <- paste("> REFERENCE SEQUENCE", index, "PART", extraIndex)
  Lines <- c(firstLine, str_chunk(toupper(seq)))
  writeLines(Lines, f)
  close(f)
  fname
}

### This function takes in a vector of 2 filenames and returns a table of read names, sequences and corresponding qualities
### If nameSuffixes = TRUE, the read names are just the last numbers preceding /1 and /2, respectively; if not, full names used
### If maxLines is specified, only that many lines are read from the input file; if clean = TRUE, any reads with N's are removed
getReadsAndQualities = function(Filenames, nameSuffixes = TRUE, maxLines = -1L, clean = TRUE) {
  File1 <- Filenames[1]
  File2 <- Filenames[2]
  U <- readLines(File1, n = maxLines)
  V <- readLines(File2, n = maxLines)
  L <- length(U)
  stopifnot(mod(L, 4) == 0)
  stopifnot(length(V) == L)
  L0 <- as.integer(L/4)
  allReads1 <- U[4 * (1:L0) - 2]
  allReads2 <- V[4 * (1:L0) - 2]
  allQuals1 <- U[4 * (1:L0)]
  allQuals2 <- V[4 * (1:L0)]
  allNames1 <- U[4 * (1:L0) - 3]
  allNames2 <- V[4 * (1:L0) - 3]
  if (nameSuffixes) {
    allNames1 %<>% str_extract('[0-9]+\\/1') %>%
      str_remove('\\/1')
    allNames2 %<>% str_extract('[0-9]+\\/2') %>%
      str_remove('\\/2')
  }
  stopifnot(all(allNames1 == allNames2))
  if (nameSuffixes) {
    allNames1 %<>%
      parse_integer
  }
  output <- tibble(name = allNames1, read1 = allReads1, qual1 = allQuals1, read2 = allReads2, qual2 = allQuals2)
  if (clean) {
    output %<>%
      filter(str_detect(read1, 'N', negate = TRUE) && str_detect(read2, 'N', negate = TRUE)) %>%
      arrange(name)
  }
  output
}

### This function extracts the reads in a numbered sequence of paired-end read files of length maxIndex starting from firstFile
### If cleanUp = TRUE, the files (including the corresponding reference file) get removed immediately after getting processed
prepareExtraReads = function(firstFile, maxIndex, extensions = paste0(1:2, ".fq"), cleanUp = FALSE) {
  if (cleanUp) {
    ind <- firstFile %>%
      str_extract("[0-9]+") %>%
      parse_integer
    myDir <- paste0("Chimera", ind)
    dir.create(myDir)
  }
  myList <- vector("list", maxIndex + 1)
  basePos <- firstFile %>%
    str_locate_all('0') %>%
    extract2(1) %>%
    tail(1) %>%
    extract(1)
  baseFileParts <- str_sub(firstFile, start = c(1, basePos + 1), end = c(basePos - 1, -1))
  for (index in 0:maxIndex) {
    curFile <- paste0(baseFileParts[1], index, baseFileParts[2])
    curFiles <- paste0(curFile, extensions)
    myList[[index + 1]] <- getReadsAndQualities(curFiles)
    if (cleanUp) {
      file.rename(curFiles, paste0(myDir, "/", curFiles))
      refFile <- paste0("Reference", ind, "I", index, ".fa")
      file.rename(refFile,  paste0(myDir, "/", refFile))
    }
  }
  output <- do.call(bind_rows, myList)
  output
}

### This is the driver function for simulating reads from an input chimera; the resulting reads are returned as an output table
simulateReads = function(curChimera, index) {
  curDir <- getwd()
  setwd(artDir)
  profs <- paste0(extDir, "/", MACHINE, "R", 1:2, ".txt")
  allCovs <- curChimera$coverage
  ### first part: simulate the full sequence at the lowest coverage
  cov1 <- min(allCovs)
  fullSeq <- str_c(curChimera$seq, collapse = "")
  fullL <- str_length(fullSeq)
  extraIndex <- 0
  ref1 <- createReference(seq = fullSeq, index, extraIndex, dir = ProjectDir)
  out1 <- paste0(ProjectDir, "Chimera", index, "Part", extraIndex, "R")
  result <- system2('./art_illumina', args = c('--paired', '--qprof1', profs[1], '--qprof2', profs[2], '--len', READ_LENGTH,
                                               '--mflen', MFL, '--sdev', SDEV, '--fcov', cov1, '--in', ref1, '--out', out1, '--rndSeed', MYSEED, '--noALN'), stdout = TRUE)
  ### second part: simulate the sequence with the segments at lowest coverage contracted to half a read length on each side
  curChimera %<>%
    mutate(breakpoint = (coverage == cov1), start = c(0L, head(cumsum(nchar(seq)), -1)) + 1L, end = cumsum(nchar(seq)))
  redChimera <- curChimera %>%
    filter(breakpoint == FALSE) %>%
    mutate(start = mod(start - HALF_READ_LENGTH, fullL), end = mod(end + HALF_READ_LENGTH, fullL), coverage = coverage - cov1)
  redSeq1 <- str_sub_circular(fullSeq, redChimera$start, redChimera$end)
  cov2 <- min(redChimera$coverage)
  print(paste("There are", length(redSeq1), "medium fragments to simulate from"))
  for (index1 in 1:length(redSeq1)) {
    print(index1)
    ref2 <- createReference(seq = redSeq1[index1], index = index, extraIndex = index1, dir = ProjectDir)
    out2 <- paste0(ProjectDir, "Chimera", index, "Part", index1, "R")
    res <- system2('./art_illumina', args = c('--paired', '--qprof1', profs[1], '--qprof2', profs[2], '--len', READ_LENGTH,
                                              '--mflen', MFL, '--sdev', SDEV, '--fcov', cov2, '--in', ref2, '--out', out2, '--rndSeed', MYSEED, '--noALN'), stdout = TRUE)
  }
  ### third part: simulate the sequence with the segments at second lowest coverage contracted to half a read length on each side
  cov12 <- cov1 + cov2
  curChimera %<>%
    mutate(breakpoint = (breakpoint | (coverage == cov12)))
  miniChimera <- curChimera %>%
    filter(breakpoint == FALSE) %>%
    mutate(start = mod(start - HALF_READ_LENGTH, fullL), end = mod(end + HALF_READ_LENGTH, fullL), coverage = coverage - cov12)
  redSeq2 <- str_sub_circular(fullSeq, miniChimera$start, miniChimera$end)
  cov3 <- min(miniChimera$coverage)
  stopifnot(all(miniChimera$coverage == cov3))
  print(paste("There are", length(redSeq2), "short fragments to simulate from"))
  for (index2 in 1:length(redSeq2)) {
    print(index2)
    extraIndex <- length(redSeq1) + index2
    ref3 <- createReference(seq = redSeq2[index2], index = index, extraIndex = extraIndex, dir = ProjectDir)
    out3 <- paste0(ProjectDir, "Chimera", index, "Part", extraIndex, "R")
    res <- system2('./art_illumina', args = c('--paired', '--qprof1', profs[1], '--qprof2', profs[2], '--len', READ_LENGTH,
                                              '--mflen', MFL, '--sdev', SDEV, '--fcov', cov3, '--in', ref3, '--out', out3, '--rndSeed', MYSEED, '--noALN'), stdout = TRUE)
  }
  maxIndex <- length(redSeq1) + length(redSeq2)
  initialFilename <- str_remove(out1, ProjectDir)
  setwd(ProjectDir)
  output <- prepareExtraReads(initialFilename, maxIndex, cleanUp = TRUE)
  setwd(curDir)
  output
}

### This auxiliary function estimates the coverage adjustment factor needed to achieve fraction of totalReads with the chimeras
computeFactor = function(totalReads, fraction = 0.1) {
  meanChimeraLength <- 5.1e6
  extraFactor <- 2
  desiredChars <- totalReads * fraction * READ_LENGTH
  myFactor <- as.integer(round(desiredChars / (meanChimeraLength * extraFactor * NUM_CHIMERAS)))
  myFactor
}

### This function operates on the specified pair of filenames (assumed, without checking, to be paired-end reads)
### It first replaces the reads indicated by replaceInds with the first |replaceInds| fakeReads (also assumed to be paired)
### It then removes the reads indicated by removeInds, and finally randomizes the order of the reads in the resulting file
### The results get written into files with "Randomized" included in the name, and the return value is the remaining fakeReads
removeAndReplaceReads = function(filenames, replaceInds, removeInds, fakeReads) {
  stopifnot(length(filenames) == 2)
  outFilenames <- str_replace(filenames, '.fastq', '_Randomized.fastq')
  numRemove <- length(removeInds)
  numReplace <- length(replaceInds)
  usedReads <- fakeReads %>%
    slice(1:numReplace)
  for (fnum in 1:2) {
    print(fnum)
    myLines <- readLines(filenames[fnum])
    numLines <- length(myLines)
    stopifnot(mod(numLines, 4) == 0)
    numReads <- as.integer(numLines/4)
    numRemain <- numReads - numRemove
    if (fnum == 1) {
      finalOrder <- sample(1:numRemain)
    }
    readLines <- 4 * ((1:numReads) - 1) + 2
    qualLines <- 4 * ((1:numReads) - 1) + 4
    ### replacement first
    curReadName <- paste0("read", fnum)
    myLines[readLines[replaceInds]] <- usedReads[[curReadName]]
    curQualName <- paste0("qual", fnum)
    myLines[qualLines[replaceInds]] <- usedReads[[curQualName]]
    ### removal second
    removeLines <- as.vector(outer(1:4, removeInds, function(x, y) {4 * (y - 1) + x}))
    myLines <- myLines[-removeLines]
    ### randomization last; note that the odd lines do not change at all!
    finalLines <- as.vector(outer(c(2, 4), 1:numRemain, function(x, y) {4 * (y - 1) + x}))
    orderLines <- as.vector(outer(c(2, 4), finalOrder, function(x, y) {4 * (y - 1) + x}))
    myLines[finalLines] <- myLines[orderLines]
    ### writing the new file
    f = file(outFilenames[fnum], 'w')
    writeLines(myLines, f)
    close(f)
  }
  fakeReads %<>%
    slice(-(1:numReplace))
  fakeReads
}

### This function injects fake reads into the files specified in fileSummaries, in randomly generated locations
### It also reduces the size of the resulting file by a percentage drawn from a normal distribution with specified mean and sd
### If removeFromEnds is positive, then the removed reads always include that many reads from each end of the file (disguise!)
injectFakeReads = function(fileSummaries, fakeReads, meanReduction = 1, sdReduction = 1, removeFromEnds = REMOVE_FROM_ENDS) {
  curDir <- getwd()
  setwd(startDir)
  N <- nrow(fileSummaries)/2
  set.seed(MYSEED)
  pcReduction <- rnorm(N, meanReduction, sdReduction)
  pcReduction %<>%
    sort
  numFakeReads <- 2 * sum(sapply(fakeReads, nrow))
  numRealReads <- sum(fileSummaries$nReads)
  fracReplace <- numFakeReads/numRealReads
  fileSummaries %<>%
    arrange(nReads) %>%
    mutate(nReplace = as.integer(round(fracReplace * nReads))) %>%
    mutate(nRemove = as.integer(round(nReads * rep(pcReduction, each = 2) / 100))) %>%
    mutate(nRetain = nReads - nRemove)
  delta <- (numFakeReads - sum(fileSummaries$nReplace))/2
  mySign <- ifelse(delta >= 0, 1, -1)
  adjust <- sample(1:N, abs(delta), replace = FALSE)
  fileSummaries %<>%
    mutate(nReplace = nReplace + mySign * (Sample %in% adjust))
  fileSummaries %<>%
    arrange(Sample, Read)
  fakeReads %<>%
    do.call(bind_rows, .)
  fakeReads %<>%
    select(-name) %>%
    mutate(number = sample(1:(numFakeReads/2))) %>%
    arrange(number)
  for (index in 1:N) {
    print(paste("Currently processing sample", index))
    curRows <- fileSummaries %>%
      filter(Sample == index)
    curSize <- curRows$nReads[1]
    curReplaceInds <- sample(1:curSize, curRows$nReplace[1]) ### generate the indices of reads to replace
    partRemoveInds <- sample(REMOVE_FROM_ENDS + (1:(curSize - 2 * REMOVE_FROM_ENDS)), curRows$nRemove[1] - 2 * REMOVE_FROM_ENDS)
    curRemoveInds <- c(seq_len(REMOVE_FROM_ENDS), partRemoveInds, curSize + 1 - rev(seq_len(REMOVE_FROM_ENDS)))
    ### note that a replaced read may later be removed
    curFiles <- curRows %>%
      pull(Filename)
    fakeReads <- removeAndReplaceReads(curFiles, curReplaceInds, curRemoveInds, fakeReads)
  }
  stopifnot(nrow(fakeReads) == 0)
  setwd(curDir)
  return(0)
}

### This function is the main driver for integrating the reads originating from a collection of chimeras into an existing sample
mainMerging = function() {
  numReads <- getNumReads(postprocessFilename = TRUE)
  totalReads <- sum(numReads$nReads)
  useFactor <- computeFactor(totalReads)
  curDir <- getwd()
  setwd(artDir)
  setwd(extDir)
  if (!all(paste0(MACHINE, "R", 1:2, ".txt") %in% list.files())) { ### check if error model is present; if not, compute it first
    prepareErrorModel()
  }
  curList <- vector("list", NUM_CHIMERAS)
  for (index in 1:NUM_CHIMERAS) {
    print(index)
    curChimera <- processChimera(index, factor = useFactor)
    curList[[index]] <- simulateReads(curChimera, index = index)
  }
  save(curList, file = "ChimericReads.RData")
  stopifnot(injectFakeReads(numReads, curList) == 0)
  curList
}

### This function extracts a named list containing the read identifiers found in the read files in useDir matching the pattern
extractIDs = function(useDir = startDir, pattern = "1.fastq$") {
  curDir <- getwd()
  setwd(useDir)
  LF <- list.files(pattern = pattern)
  numFiles <- length(LF)
  outList <- vector("list", numFiles)
  redNames <- str_split_fixed(LF, "_", n = 6) %>%
    extract(, c(3, 6), drop = FALSE) %>%
    apply(1, function(x) {str_c(x, collapse = "_")})
  names(outList) <- redNames
  for (ind in 1:length(LF)) {
    fname <- LF[ind]
    print(fname)
    curLines <- readLines(fname)
    nReads <- length(curLines)/4
    curIDs <- curLines[(1:nReads) * 4 - 3]
    outList[[redNames[ind]]] <- curIDs
  }
  setwd(curDir)
  outList
}

### This function explores the ranges of allowed or observed variation in the identifier parts extracted from a specific sample
getRanges = function(ListOfIDs, numParts = 7) {
  ListOfIDs %<>%
    lapply(function(x) {str_split_fixed(x, ":", n = numParts)})
  ListOfIDs %<>%
    lapply(function(x) {colnames(x) <- paste0("Part", 1:numParts); x %<>% as_tibble})
  ListOfIDs <- lapply(1:length(ListOfIDs), function(y) {x <- ListOfIDs[[y]]; x %<>% mutate(index = y)})
  ListOfIDs %<>%
    do.call(bind_rows, .)
  numUnique <- summarize_all(ListOfIDs, ~n_distinct(.))
  ListOfIDs %<>%
    select(which(numUnique > 1))
  numUnique %<>%
    select(which(numUnique > 1))
  M <- max(ListOfIDs[['index']])
  for (ind in 1:(ncol(numUnique) - 1)) {
    curCol <- ListOfIDs[[ind]]
    if (numUnique[[ind]] == M) {
      myTab <- table(curCol, ListOfIDs[['index']])
      stopifnot(sum(diag(myTab)) == sum(myTab))
    }
  }
  badCols <- head(which(numUnique[1,] == max(ListOfIDs[['index']])), -1)
  ListOfIDs %<>%
    select(-badCols)
  numUnique %<>%
    select(-badCols)
  ListOfIDs %<>%
    mutate(Part7 = str_remove(Part7, '/1')) %>%
    mutate_if(is.character, parse_integer)
  ListOfIDs %<>%
    arrange(index, Part5, Part6, Part7)
  ListOfIDs %<>%
    as.matrix
  testMat <- ListOfIDs[-1,] - ListOfIDs[-nrow(ListOfIDs),]
  numRepeatIndices <- sum(rowSums(testMat[,-1] == 0) == ncol(testMat) - 1)
  stopifnot(numRepeatIndices == 0)
  Part5 <- ListOfIDs[,'Part5']
  Part5 %<>%
    as.character %>%
    str_split_fixed("", n = 4)
  index <- ListOfIDs[,'index']
  Part5Ranges <- sapply(split(Part5, index), 
                        function(x) {
                          apply(matrix(as.numeric(x), ncol = 4), 2, range)
                        }
  )
  Part6 <- ListOfIDs[,'Part6']
  Part6Ranges <- sapply(split(Part6, index), range)
  Part7 <- ListOfIDs[,'Part7']
  Part7Ranges <- sapply(split(Part7, index), range)
  output <- list(Part5Ranges, Part6Ranges, Part7Ranges)
  output
}

### This function creates read files based on an input table; it must have an index column to separate the reads into file pairs
### The other columns should be name (identifier), read1 and read2 (the first & second mate) and qual1 and qual2 (the qualities)
makeReadFiles = function(inputTab, baseFilename) {
  allInds <- inputTab$index
  minIndex <- min(allInds)
  maxIndex <- max(allInds)
  for (curIndex in minIndex:maxIndex) {
    curTab <- inputTab %>%
      filter(index == curIndex)
    curSize <- nrow(curTab)
    curLines1 <- rep("+", 4 * curSize) ### This was incorrect before - the 4 was missing!
    curLines2 <- rep("+", 4 * curSize) ### This was incorrect before - the 4 was missing!
    curIDLines <- 4 * ((1:curSize) - 1) + 1
    curLines1[curIDLines] <- paste0(curTab$name, "/1")
    curLines2[curIDLines] <- paste0(curTab$name, "/2")
    curReadLines <- 4 * ((1:curSize) - 1) + 2
    curLines1[curReadLines] <- curTab$read1
    curLines2[curReadLines] <- curTab$read2
    curQualLines <- 4 * ((1:curSize) - 1) + 4
    curLines1[curQualLines] <- curTab$qual1
    curLines2[curQualLines] <- curTab$qual2
    curFilenames <- paste0(baseFilename, curIndex, 'R', 1:2, '.fastq')
    for (fnum in 1:2) {
      f = file(curFilenames[fnum], 'w')
      if (fnum == 1) {
        writeLines(curLines1, f)
      } else {
        writeLines(curLines2, f)
      }
      close(f)
    }
  }
  return(0)
}

### The driver for creating a sample out of simulated reads (known composition) together with plausible-looking read identifiers
### The RNG seed may be specified; the other arguments denote the names of intermediate and final result files
altProcessing = function(seed = MYSEED, intFile = "RangeData.RData", IDFile = "ReadIDs.RData", sFile = "SimulatedFirstReadsL00") {
  set.seed(seed)
  if (intFile %in% list.files()) {
    load(intFile)
  } else {
    testIDs <- extractIDs()
    allRanges <- getRanges(testIDs)
    save(allRanges, file = intFile)
    rm(testIDs)
  }
  part123 <- "@NB501607:36:HHY2LBGXY"
  numFiles <- 4
  ### strategy: 1) randomly generate the 4 file sizes 2) for each one, randomly generate identifier parts 3) sort and assign
  if (IDFile %in% list.files()) {
    load(IDFile)
    fileSizes <- table(allParts$part4)
  } else {
    numReads  <- getNumReads(paste0(simGenomeDir, "/OriginalFiles"), pattern = "*fullyCleaned*")
    N <- sum(pmin(numReads$nReads, MAX_L/4))/2
    readInds <- sample(1:numFiles, N, replace = TRUE, prob = 0.25 + runif(numFiles)/100)
    fileSizes <- table(readInds)
    part4 <- rep(1:numFiles, fileSizes)
    part5 <- matrix(NA, nrow(allRanges[[1]])/2, N)
    part6 <- rep(NA, N)
    part7 <- rep(NA, N)
    curPos <- 0
    for (index in 1:numFiles) {
      curSpan <- curPos + (1:fileSizes[index])
      curRanges5 <- matrix(allRanges[[1]][, index], nrow = 2)
      for (ind in 1:ncol(curRanges5)) {
        curSubRange5 <- curRanges5[1, ind]:curRanges5[2, ind]
        part5[ind, curSpan] <- sample(curSubRange5, fileSizes[index], replace = TRUE)
      }
      curRanges6 <- allRanges[[2]][, index]
      offset6 <- as.integer(round(runif(1) * (curRanges6[1] - 1000L)))
      curSubRange6 <- (curRanges6[1]:curRanges6[2]) - offset6
      part6[curSpan] <- sample(curSubRange6, fileSizes[index], replace = TRUE)
      curRanges7 <- allRanges[[3]][, index]
      offset7 <- as.integer(round(runif(1) * (curRanges7[1] - 1000L)))
      curSubRange7 <- (curRanges7[1]:curRanges7[2]) - offset7
      part7[curSpan] <- sample(curSubRange7, fileSizes[index], replace = TRUE)
      curPos %<>%
        add(fileSizes[index])
    }
    Part5 <- part5 %>%
      multiply_by(c(10000L, 1000L, 100L, 1L)) %>%
      colSums
    allParts <- tibble(part123 = part123, part4 = part4, part5 = Part5, part6 = as.integer(part6), part7 = as.integer(part7)) %>%
      arrange(part123, part4, part5, part6, part7) %>%
      mutate(allIDs = paste(part123, part4, part5, part6, part7, sep=":"))
    save(allParts, file = IDFile)
  }
  curDir <- getwd()
  setwd(simGenomeDir)
  ### see the way I used to do this previously in LegacyCode.R
  LF1 <- list.files(pattern = "1.fastq$")
  LF2 <- str_replace(LF1, '1.fastq', '2.fastq')
  for (ind in 1:nrow(numReads)) {
    curRow <- numReads %>%
      slice(ind)
    curFile <- curRow %>%
      pull(Filename)
    altFile <- str_c(str_sub(curFile, end = -7), 'Reduced', str_sub(curFile, start = -7))
    if (curRow$nReads > MAX_L/4) {
      system2("split", args = c("-a", 1L, "-l", MAX_L, curFile, altFile))
      file.rename(from = paste0(altFile, "a"), to = curFile)
    }
  }
  system2("cat", args = c(LF1, ">", "allFirstR1.fastq"))
  system2("cat", args = c(LF2, ">", "allSecondR2.fastq"))
  names(fileSizes) <- NULL
  splitPoints <- fileSizes %>%
    cumsum %>%
    multiply_by(4) %>%
    add(1) %>%
    head(-1)
  uBound <- numFiles - 1L
  system2("csplit", args = c("-f", sFile, "-n", 1L, "-s", "allFirstR1.fastq",  splitPoints))
  file.rename(from = paste0(sFile, 0:uBound), to = paste0(sFile, 0:uBound, "R1_1.fastq"))
  system2("csplit", args = c("-f", sFile, "-n", 1L, "-s", "allSecondR2.fastq", splitPoints))
  file.rename(from = paste0(sFile, 0:uBound), to = paste0(sFile, 0:uBound, "R2_2.fastq"))
  for (index in 1:numFiles) {
    fnames <- paste0(sFile, index - 1L, "R", 1:2, "_", 1:2, ".fastq")
    print(fnames)
    curInput <- getReadsAndQualities(fnames, nameSuffixes = FALSE, maxLines = MAX_L, clean = FALSE)
    curIDs <- allParts %>%
      filter(part4 == index) %>%
      pull(allIDs)
    curInput %<>%
      mutate(name = curIDs, index = index)
    makeReadFiles(curInput, "SimulatedFinalReadsL00")
  }
  setwd(curDir)
}

fixReadFiles = function(fnames = list.files()) {
  for (fname in fnames) {
    print(fname)
    curLines <- readLines(fname)
    L <- length(curLines)
    stopifnot(mod(L, 4) == 0)
    L0 <- L/4
    curLines[4 * (1:L0) - 1] <- "+"
    f = file(paste0("Fixed", fname), "w")
    writeLines(curLines, f)
    close(f)
  }
}
