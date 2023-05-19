#!/usr/bin/env Rscript
library(magrittr)

#===========================================================================####
# Arguments
#===========================================================================####
# 1 - merged boot data as .tsv
args <- c("/projects/b1059/projects/Tim/DFE_boots-nf/BOOT_results_20230519/data/mergedBoots.tsv")
#args <- commandArgs(trailingOnly = TRUE)

#===========================================================================####
# Process the bootstrap results
#===========================================================================####
# load parent genos
load("/projects/b1059/projects/Tim/DFE_boots-nf/input_data/parent_genos.rda")

# load the boot strap data
bb <- data.table::fread(args[1]) %>%
  dplyr::filter(geno_meth == "b")

baf <- data.table::fread(args[1]) %>%
  dplyr::filter(geno_meth == "af")

# calculate the bootstrap stats for both methods
bb2 <- bb %>%
  dplyr::group_by(loci) %>%
  dplyr::mutate(bb_mean_mut_effect_res = mean(mut_effect_res),
                lower.conf = unname(quantile(mut_effect_res, probs = 0.05)), #0.00029586
                upper.conf = unname(quantile(mut_effect_res, probs = (1-0.05)))) %>%
  dplyr::distinct(loci, .keep_all = T) %>%
  dplyr::mutate(effect_type = case_when(lower.conf > 0 & upper.conf > 0 ~ "positive",
                                        lower.conf < 0 & upper.conf < 0 ~ "negative",
                                        TRUE ~ "neutral")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(overall_bb_mean_mut_effect_res = mean(bb_mean_mut_effect_res)) %>%
  dplyr::arrange(bb_mean_mut_effect_res) %>%
  dplyr::left_join(., parent_genos) %>%
  dplyr::mutate(color = ifelse(MA530 == 1, "#1B9E77", "#7570B3"),
                shape = ifelse(MA530 == 1, "MA530", "MA563"))

baf2 <- baf %>%
  dplyr::group_by(loci) %>%
  dplyr::mutate(baf_mean_mut_effect_res = mean(mut_effect_res),
                lower.conf = unname(quantile(mut_effect_res, probs = 0.05)), #0.00029586
                upper.conf = unname(quantile(mut_effect_res, probs = (1-0.05)))) %>%
  dplyr::distinct(loci, .keep_all = T) %>%
  dplyr::mutate(effect_type = case_when(lower.conf > 0 & upper.conf > 0 ~ "positive",
                                        lower.conf < 0 & upper.conf < 0 ~ "negative",
                                        TRUE ~ "neutral")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(overall_baf_mean_mut_effect_res = mean(baf_mean_mut_effect_res)) %>%
  dplyr::arrange(baf_mean_mut_effect_res) %>%
  dplyr::left_join(., parent_genos) %>%
  dplyr::mutate(color = ifelse(MA530 == 1, "#1B9E77", "#7570B3"),
                shape = ifelse(MA530 == 1, "MA530", "MA563"))

# plot these
bb_eff_p <- ggplot(bb2) +
  aes(x = factor(loci, levels = loci), y = bb_mean_mut_effect_res,
      ymin = lower.conf, ymax = upper.conf, color = effect_type, shape = shape) +
  geom_pointrange(size = 0.3) +
  scale_color_manual(values = c("positive" = "blue", "negative" = "red", "neutral" = "grey40")) +
  labs(x = "", y = "Mutational effect", color = "effect type", subtitle = "random replacement method") +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 4)) # axis.text.y = element_text(color = bb2$color)
bb_eff_p

baf_eff_p <- ggplot(baf2) +
  aes(x = factor(loci, levels = loci), y = baf_mean_mut_effect_res,
      ymin = lower.conf, ymax = upper.conf, color = effect_type, shape = shape) +
  geom_pointrange(size = 0.3) +
  scale_color_manual(values = c("positive" = "blue", "negative" = "red", "neutral" = "grey40")) +
  labs(x = "", y = "Mutational effect", color = "effect type", subtitle = "AF replacement method") +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 4))
baf_eff_p

b.both.p <- cowplot::plot_grid(bb_eff_p, baf_eff_p, ncol = 2)

# save it
cowplot::ggsave2(b.both.p, filename = "/projects/b1059/projects/Tim/DFE_boots-nf/BOOT_results_20230519/data/mutational.effect.boot.CIs.png", width = 7, height = 10)

# save data
proc.mut.effect.boots <- list(bb2, baf2)
save(proc.mut.effect.boots, file = "/projects/b1059/projects/Tim/DFE_boots-nf/BOOT_results_20230519/data/proc.mut.effect.boots.rda")

#===========================================================================####
# Plot distribution
#===========================================================================####
bb_dist_p <- ggplot(bb2) +
  aes(x = bb_mean_mut_effect_res) +
  geom_histogram(bins = 30) +
  labs(y = "n loci", x = "Mutational effect", subtitle = "random replacement method") +
  geom_vline(xintercept = unique(bb2$overall_bb_mean_mut_effect_res), linetype = 2) +
  annotate(geom = "label", x = 0.01, y = 10, label = glue::glue("mean eff. = {round(unique(bb2$overall_bb_mean_mut_effect_res), digits = 6)}")) +
  theme_bw() +
  labs(x = glue::glue("Mutational effect"))
bb_dist_p

baf_dist_p <- ggplot(baf2) +
  aes(x = baf_mean_mut_effect_res) +
  geom_histogram(bins = 30) +
  labs(y = "n loci", x = "Mutational effect", subtitle = "AF replacement method") +
  geom_vline(xintercept = unique(baf2$overall_baf_mean_mut_effect_res), linetype = 2) +
  annotate(geom = "label", x = 0.01, y = 10, label = glue::glue("mean eff. = {round(unique(baf2$overall_baf_mean_mut_effect_res), digits = 6)}")) +
  theme_bw() +
  labs(x = glue::glue("Mutational effect"))
baf_dist_p

b.both.dist.p <- cowplot::plot_grid(bb_dist_p, baf_dist_p, ncol = 2)

# save it
cowplot::ggsave2(b.both.dist.p, filename = "/projects/b1059/projects/Tim/DFE_boots-nf/BOOT_results_20230519/data/mutational.effect.boot.CIs.dist.png", width = 10, height = 5)

#===========================================================================####
# Plot point estimate correlation between methods
#===========================================================================####
cor.df <- bb2 %>%
  dplyr::select(loci, bb_mean_mut_effect_res) %>%
  dplyr::left_join(dplyr::select(baf2, loci, baf_mean_mut_effect_res))

cor.p <- ggplot(cor.df) +
  aes(x = bb_mean_mut_effect_res, baf_mean_mut_effect_res) +
  geom_point(shape = 21) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
  theme_bw() +
  labs(x = "mutational effect estimate\n(random replacement)", y = "mutational effect estimate\n(AF replacement)")
cor.p

cowplot::ggsave2(cor.p, filename = "/projects/b1059/projects/Tim/DFE_boots-nf/BOOT_results_20230519/data/mutational.effect.boot.method.corr.png", width = 5, height = 5)
