# Chimerize
This repo contains the code used to create chimeric and simulated metagenomic samples

Processing.R contains the main scripts to:
  - process a database of reference genomes
  - find good chimera candidates
  - create draft chimeras
  
Merging.R contains the main scripts to:
  - generate reads from chimeras
  - combine those reads with an existing sample
  - create plausible read IDs for a simulated sample
