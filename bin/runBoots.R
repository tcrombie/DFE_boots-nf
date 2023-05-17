#!/usr/bin/env Rscript
library(magrittr)
library(dplyr)
library(rio)

#===========================================================================####
# Arguments
#===========================================================================####
# 1 - path to the configBoots.rda file
# 2 - the boot random seed
# 3 - out dir
#args <- c("/projects/b1059/projects/Tim/DFE_boots-nf/test/data/configBoots.rda", 100, "/projects/b1059/projects/Tim/DFE_boots-nf/test")
args <- commandArgs(trailingOnly = TRUE)

#===========================================================================####
# Load data, sample lines, fill in genotypes, clac residual fits, ect.
#===========================================================================####
# load configBoots list
load(args[1])

# pull out elements
data <- configBoots[[1]]
lines <- configBoots[[2]]
afs <- configBoots[[3]]

# get loci
loci <- afs$loci

# 2:sample the lines with replacement and calculate mut_effects
# set the seed
set.seed(args[2])


# sample lines
line.samp <- sample(lines, replace = TRUE)
  
# build full dataset with these lines
dat.samp.list <- NULL
  for(j in 1:length(line.samp)){
    # get the line name
    line <- line.samp[j]
    # get rep data
    dat.samp <- data %>%
      dplyr::filter(full_id == line)
    # add to list
    dat.samp.list[[j]] <-  dat.samp
  }
  
  # bind sampled line data together and assign new id
  dat.samp.df <- data.table::rbindlist(dat.samp.list, idcol = "line.samp.id")
  
  # randomly assign genotypes where missing w/ 50/50 binomial and allele freq binomial sampling
  dat.samp.df2 <- dat.samp.df %>%
    dplyr::distinct(line.samp.id, .keep_all = T) %>%
    tidyr::gather(variant, genotype, -1:-5) %>% # reshape the variant data to long format
    dplyr::left_join(afs) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(geno.b = case_when(is.na(genotype) ~ rbinom(n = 1, size = 1, prob = 0.5),
                                     TRUE ~ genotype),
                  geno.af = case_when(is.na(genotype) ~ rbinom(n = 1, size = 1, prob = mut.af),
                                      TRUE ~ genotype)) %>%
    dplyr::ungroup() 
  
  # get random genos then calculate residuals
  geno.b.df <- dat.samp.df2 %>%
    dplyr::select(line.samp.id, variant, geno.b) %>%
    tidyr::spread(variant, geno.b) %>%
    dplyr::left_join(dplyr::select(dat.samp.df, 1:5))
  geno.b.df.dist <- geno.b.df %>%
    dplyr::distinct(line.samp.id, .keep_all = T) %>%
    dplyr::select(1:(max(loci)+1)) %>%
    tidyr::pivot_longer(cols = 2:(max(loci)+1), names_to = "loci", values_to = "geno.b")
  
  # get af genos
  geno.af.df <- dat.samp.df2 %>%
    dplyr::select(line.samp.id, variant, geno.af) %>%
    tidyr::spread(variant, geno.af) %>%
    dplyr::left_join(dplyr::select(dat.samp.df, 1:5))
  geno.af.df.dist <- geno.af.df %>%
    dplyr::distinct(line.samp.id, .keep_all = T) %>%
    dplyr::select(1:(max(loci)+1)) %>%
    tidyr::pivot_longer(cols = 2:(max(loci)+1), names_to = "loci", values_to = "geno.af")
  
  # make lists to hold data
  geno.b.res.list <- NULL
  geno.af.res.list <- NULL
  for(k in 1:length(loci)){
    # get the proper index offset by 1
    all <- 1+(1:length(loci))
    k.index <- all[-(k)]
    # get b res
    res.b = residuals(lm(data = geno.b.df[, c(k.index, (length(loci)+5))], ln_ci ~ .))
    # get af res
    res.af = residuals(lm(data = geno.af.df[, c(k.index, (length(loci)+5))], ln_ci ~ .))
    
    # assign them to list
    geno.b.res.list[[k]] <- res.b
    geno.af.res.list[[k]] <- res.af
  }
  
  # put them all together, calculate line.sample.means for the residuals
  geno.b.res.df <- as.data.frame(do.call(cbind, geno.b.res.list))
  names(geno.b.res.df) <- paste0(afs$variant)
  
  # make for binomial
  geno.b.res.df2 <- geno.b.res.df %>%
    dplyr::bind_cols(geno.b.df[-2:-(max(loci)+1)]) %>%
    tidyr::pivot_longer(cols = 1:max(loci), names_to = "loci", values_to = "loci_res") %>%
    dplyr::group_by(loci, line.samp.id) %>%
    dplyr::mutate(line.samp.res.mean = mean(loci_res)) %>%
    dplyr::distinct(loci, line.samp.id, .keep_all = T) %>%
    dplyr::left_join(., geno.b.df.dist) %>%
    dplyr::select(-rep) %>% # not useful anymore
    dplyr::group_by(loci, geno.b) %>%
    dplyr::mutate(mean_loci_res = mean(line.samp.res.mean)) %>% # get a mean value for the loci for both genotypes 
    dplyr::ungroup() %>%
    dplyr::group_by(loci) %>%
    dplyr::mutate(mut_mean_loci_res = case_when(geno.b == 1 ~ mean_loci_res,
                                                geno.b == 0 ~ NA_real_,
                                                is.na(geno.b) ~ NA_real_),
                  wt_mean_loci_res = case_when(geno.b == 1 ~ NA_real_,
                                               geno.b == 0 ~ mean_loci_res,
                                               is.na(geno.b) ~ NA_real_)) %>%
    tidyr::fill(mut_mean_loci_res, .direction = "downup") %>%
    tidyr::fill(wt_mean_loci_res, .direction = "downup") %>%
    dplyr::distinct(loci, mut_mean_loci_res, wt_mean_loci_res) %>%
    dplyr::mutate(mut_effect_res = mut_mean_loci_res - wt_mean_loci_res) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(boot.id = as.numeric(args[2]))
  
  # put them all together, calculate line.sample.means for the residuals
  geno.af.res.df <- as.data.frame(do.call(cbind, geno.af.res.list))
  names(geno.af.res.df) <- paste0(afs$variant)
  
  # make for allele freq
  geno.af.res.df2 <- geno.af.res.df %>%
    dplyr::bind_cols(geno.af.df[-2:-(max(loci)+1)]) %>%
    tidyr::pivot_longer(cols = 1:max(loci), names_to = "loci", values_to = "loci_res") %>%
    dplyr::group_by(loci, line.samp.id) %>%
    dplyr::mutate(line.samp.res.mean = mean(loci_res)) %>%
    dplyr::distinct(loci, line.samp.id, .keep_all = T) %>%
    dplyr::left_join(., geno.af.df.dist) %>%
    dplyr::select(-rep) %>% # not useful anymore
    dplyr::group_by(loci, geno.af) %>%
    dplyr::mutate(mean_loci_res = mean(line.samp.res.mean)) %>% # get a mean value for the loci for both genotypes 
    dplyr::ungroup() %>%
    dplyr::group_by(loci) %>%
    dplyr::mutate(mut_mean_loci_res = case_when(geno.af == 1 ~ mean_loci_res,
                                                geno.af == 0 ~ NA_real_,
                                                is.na(geno.af) ~ NA_real_),
                  wt_mean_loci_res = case_when(geno.af == 1 ~ NA_real_,
                                               geno.af == 0 ~ mean_loci_res,
                                               is.na(geno.af) ~ NA_real_)) %>%
    tidyr::fill(mut_mean_loci_res, .direction = "downup") %>%
    tidyr::fill(wt_mean_loci_res, .direction = "downup") %>%
    dplyr::distinct(loci, mut_mean_loci_res, wt_mean_loci_res) %>%
    dplyr::mutate(mut_effect_res = mut_mean_loci_res - wt_mean_loci_res) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(boot.id = as.numeric(args[2]))

#===========================================================================####
# Load data, sample lines, fill in genotypes, clac residual fits, ect.
#===========================================================================####
# save dat
rio::export(geno.af.res.df2, file = paste0("af_mutEffect_", args[2], ".tsv"))
rio::export(geno.b.res.df2, file = paste0("b_mutEffect_", args[2], ".tsv"))

