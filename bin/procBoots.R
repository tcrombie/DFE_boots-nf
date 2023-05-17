#!/usr/bin/env Rscript
library(magrittr)

#===========================================================================####
# Arguments
#===========================================================================####
# 1 - path to the configBoots.rda file
# 2 - the boot random seed
# 3 - random seed #
# 4 - out dir
args <- c("/projects/b1059/projects/Tim/DFE_boots-nf/test/data/configBoots.rda", 100, 1, "/projects/b1059/projects/Tim/DFE_boots-nf/test")
#args <- commandArgs(trailingOnly = TRUE)
