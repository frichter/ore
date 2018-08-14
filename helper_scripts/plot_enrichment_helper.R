# Felix Richter (felix.richter@mssm.edu)
# Plot different combinations of enrichment
#############################################

options(stringsAsFactors=FALSE)
p = c("magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "readr")

## Run once: install packages
lapply(p, install.packages)

## Load packages
lapply(p, require, character.only = TRUE)

## optionally, install interesting color pallettes
library(wesanderson)

## enrichment file directory
enrich_file_dir = "/Users/felixrichter/Dropbox/PhD/alzheimers/enrichment"

## Load data
enrich_file_list = list.files(enrich_file_dir, "_enrichment.txt", full.names = T)
## clean names of enrichment runs
names(enrich_file_list) = gsub(".*/", "", enrich_file_list) %>% gsub(".txt$", "", .) %>% 
  gsub("ad_ore_", "", .) %>% gsub("_10kb", "", .)
enrich_df = map_df(enrich_file_list, read_tsv, .id = "enrichment_run")

## Plot data
p = enrich_df %>% 
  filter(!is.na(gene_ci_lo)) %>% 
  ## Choose expression outlier cut-offs to plot, e.g., (0, 1) for custom:
  filter(expr_cut_off %in% c(2, 2.5, 3)) %>% 
  ## Choose allele frequency cut-off to plot:
  filter(af_cut_off %in% c(1e-4, 1e-3, 1e-2, 0.05)) %>% 
  ## Choose TSS distance cut-off for plots:
  filter(tss_cut_off %in% c(1e4)) %>% 
  ## Choose enrichment run (file name without the txt ending)
  filter(enrichment_run %in% c("utr5", "upstream", "exonic", "splicing", "intergenic", "intronic")) %>% 
  ## convert the value you are coloring by to a factor
  mutate(expr_cut_off = factor(expr_cut_off)) %>% 
  ## pick x-axis: expr_cut_off, af_cut_off, tss_cut_off, enrichment_run:
  ggplot(., aes(x = af_cut_off,
                y = gene_or, ymax = gene_ci_hi, ymin = gene_ci_lo,
  ## pick variable for COLORS: expr_cut_off, af_cut_off, tss_cut_off, enrichment_run
                col = expr_cut_off
                )) +
  ## x-axis log scale only for af_cut_off:
  scale_x_log10() +
  geom_hline(yintercept = 1, col = "grey60") +
  geom_pointrange(fatten = 1, show.legend = T, position = position_dodge(width = 0.3)) + 
  theme_classic() +
  ## Choose colors
  scale_color_manual(values=wes_palette(n=5, name="Zissou1")[c(1,3,5)]) + ##
  facet_wrap(~enrichment_run) +
  xlab("Allele frequency cut-off") + ylab("Odds ratio")
p

## Save the plot
filename = paste0(enrich_file_dir, "/exonic_enrichment_10kb.png")
ggsave(filename, p, width = 3.5, height = 2.25) 


  
  

