#!/usr/bin/env Rscript
library(magrittr)
library(dplyr)

#===========================================================================####
# Arguments
#===========================================================================####
# 1 - full path to input data
# 2 - number of boot straps
# 3 - random seed #
# 4 - out dir
#args <- c("/projects/b1059/projects/Tim/DFE_boots-nf/input_data/06_DFE_joined_imputed_geno_10kbMax_pheno.csv", 100, 1, "/projects/b1059/projects/Tim/DFE_boots-nf/test")
args <- commandArgs(trailingOnly = TRUE)

#===========================================================================####
# Setup output dir                                 
#===========================================================================####
dir.create(path = args[4])

#===========================================================================####
# Load data, config boots, output config file                                 
#===========================================================================####
# 1: get the joined phenotype and 10Kb imputed genotype data then filter to RILs and remove lines with all missing genotypes
data <- data.table::fread(args[1]) %>%
  dplyr::filter(grepl(full_id, pattern = "RIL_|RIAIL_")) %>%
  dplyr::filter(rowSums(is.na(.)) != 169) %>% # remove lines with all missing genotypes
  dplyr::select(-p_focal, -p_comp)

# 2: get a line vector
lines <- data %>%
  dplyr::distinct(full_id) %>%
  dplyr::pull(full_id)

# 3: get all allele freqs.
afs <- data %>%
  dplyr::distinct(full_id, .keep_all = T) %>%
  tidyr::gather(variant, genotype, -1:-4) %>% # reshape the variant data to long format
  dplyr::group_by(variant) %>%
  dplyr::mutate(mut.af = sum(genotype, na.rm = T)/sum(!is.na(genotype))) %>%
  dplyr::distinct(variant, mut.af) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(loci = 1:n())

# 4: get loci
loci <- afs$loci

# 5: setup the boot seeds
set.seed(as.numeric(args[3]))
seeds <- sample(1:as.numeric(args[2]), replace = FALSE)

# 6: make config file
config <- tibble::tibble(boots = seeds)

# 7: write the config file
#dir.create(path = paste0(args[4], "/data"))
write.table(config, file = "boots.tsv", quote=FALSE, sep='\t', row.names = F)

# 8: write the configBoots rda file with everything for the boots
configBoots <- list(data, lines, afs)
save(configBoots, file = "configBoots.rda")
