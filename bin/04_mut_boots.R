library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get the joined phenotype and 10Kb imputed genotype data then filter to RILs and remove lines with all missing genotypes
d1 <- data.table::fread("data/processed/06_DFE_joined_imputed_geno_10kbMax_pheno.csv") %>%
  dplyr::filter(grepl(full_id, pattern = "RIL_|RIAIL_")) %>%
  dplyr::filter(rowSums(is.na(.)) != 169) %>% # remove lines with all missing genotypes
  dplyr::select(-p_focal, -p_comp)

############################################################
# DEFINE ARGUMENTS
data <- d1
i <- "RIL_259"

# testing function
undebug(mutEffBoot2)
debug(mutEffBoot2)
############################################################

test <- mutEffBoot2(data = d1, nboots=3)

# get data for charlie
test.random.df <- test[[1]] %>%
  dplyr::mutate(type = "rand")
test.allele.freq.df <- test[[2]] %>%
  dplyr::mutate(type = "af")

# save for charlie
rio::export(test.random.df, file = "data/processed/test.random.boots.csv")
rio::export(test.allele.freq.df, file = "data/processed/test.allele.freq.boots.csv")

# bind them
all <- rbind(test.random.df, test.allele.freq.df)

# explore distributions
ggplot(all) +
  aes(x = mut_effect_res) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~boot.id*type, ncol = 2)

# look at estimates
mean.additive.mut.eff.rand  <- test.random.df %>%
  dplyr::group_by(boot.id) %>%
  dplyr::mutate(mean.additive.mut.eff = mean(mut_effect_res)) %>%
  dplyr::distinct(mean.additive.mut.eff) %>%
  dplyr::pull(mean.additive.mut.eff)

############################################################
# define the function
mutEffBoot2 <- function(data, nboots=3, seed = 99){
  
  # 1: get a line vector
  lines <- data %>%
    dplyr::distinct(full_id) %>%
    dplyr::pull(full_id)
  
  # get all allele freqs.
  afs <- data %>%
    dplyr::distinct(full_id, .keep_all = T) %>%
    tidyr::gather(variant, genotype, -1:-4) %>% # reshape the variant data to long format
    dplyr::group_by(variant) %>%
    dplyr::mutate(mut.af = sum(genotype, na.rm = T)/sum(!is.na(genotype))) %>%
    dplyr::distinct(variant, mut.af) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(loci = 1:n())
  
  # get loci
  loci <- afs$loci
  
  # 2:sample the lines with replacement and calculate mut_effects
  # set the seed
  set.seed(seed)
  # setup lists
  mut.effs.b <- NULL
  mut.effs.af <- NULL
  # loop over bootstraps
  for (i in 1:nboots) {
    # message
    message(glue::glue("starting boot sample {i} of {nboots}"))
    
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
      dplyr::left_join(select(dat.samp.df, 1:5))
    geno.b.df.dist <- geno.b.df %>%
      dplyr::distinct(line.samp.id, .keep_all = T) %>%
      dplyr::select(1:(max(loci)+1)) %>%
      tidyr::pivot_longer(cols = 2:(max(loci)+1), names_to = "loci", values_to = "geno.b")
    
    # get af genos
    geno.af.df <- dat.samp.df2 %>%
      dplyr::select(line.samp.id, variant, geno.af) %>%
      tidyr::spread(variant, geno.af) %>%
      dplyr::left_join(select(dat.samp.df, 1:5))
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
      dplyr::mutate(mut_effect_res = mut_mean_loci_res - wt_mean_loci_res)
    
    # assign this to the boots
    # put these in var output
    mut.effs.b[[i]] <- geno.b.res.df2
    
    #############################################################
    
    # put them all together, calculate line.sample.means for the residuals
    geno.af.res.df <- as.data.frame(do.call(cbind, geno.af.res.list))
    names(geno.af.res.df) <- paste0(afs$variant)
    
    # make for binomial
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
      dplyr::mutate(mut_effect_res = mut_mean_loci_res - wt_mean_loci_res)
    
    # put this together
    mut.effs.af[[i]] <- geno.af.res.df2
    
  }
  # bind mut effect data together with bootstrap id
  mut.effs.b.df <- data.table::rbindlist(mut.effs.b, idcol = "boot.id")
  mut.effs.af.df <- data.table::rbindlist(mut.effs.af, idcol = "boot.id")
  
  # return dat
  return(list(mut.eff.b = mut.effs.b.df, mut.effect.af = mut.effs.af.df))
  
}