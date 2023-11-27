##### Load libraries ##############################

library(gtools)
library(phyloseq)
library(ampvis2)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(ggpubr)
library(MicEco)
library(reshape2)


##### 16S PRE-PROCESSING ##########################
##### RAW DATA -- LOAD AND FORMATTING ############
## OTU IMPORT and TAXONOMY ##

## NOTE: original count/tax tables include another SO248 experiment (DOM addition)
## DADA was run together, since experiment based on the same original seawater 
## This also allows subsequent comparisons, as consistent ASV names were generated
## These can then be tracked over experiments with alginate or DOM additions
## Here, we will subset to alginate-only samples


# Read ASVs (alginate+DOM experiment)
ASV <- read.table(
  "Alg-DOM_v138.1_seqtab.txt",
  h = T,
  row.names = 1,
  check.names = F,
  sep = "\t")

TAX <- read.table(
  "Alg-DOM_v138.1_tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1)

# Import metadata (only alginate experiment)
ENV <- read.table(
  "metadata_edited.txt", 
  h=T, sep = "\t", 
  check.names = F,
  stringsAsFactors=F)
row.names(ENV) = ENV$clip_id

# Subset to alginate samples
ASV <- ASV[names(ASV) %in% ENV$clip_id]

# Mean ASV counts in negative controls 
ASV$negative <- rowMeans(ASV[,c(
  "MQ181221-146_clip",
  "MQ181221-159_clip",
  "MQ181221-136_clip")])
ASV$negative <- round(ASV$negative, 0)

# Subtract negative counts from each column
ASV = ASV - ASV$negative

# If < 0: set to zero
ASV[ASV < 0] <- 0

# Remove more likely contaminants
TAX <- TAX[-grep(
  'Corynebacteriaceae|Bacillaceae|
  Carnobacteriaceae|Burkholderiaceae|
  Streptococcaceae|Propionibacteriaceae|
  Weeksellaceae|Enterococcaceae', TAX$Family),] 

# match TAX after contaminant removal
ASV <- ASV[row.names(TAX),]

# remove NegCtr columns (not needed anymore) and identified outlier (MQ181221-127_clip)
# did not remove other outliers on the basis of read counts- those are removed in rarefaction step
ASV <- ASV %>% dplyr::select(-c(
  "MQ181221-146_clip",
  "MQ181221-159_clip",
  "MQ181221-136_clip",
  "negative",
  "MQ181221-127_clip"))

# Remove NegCtr from metadata (not needed anymore)
ENV <- ENV %>%
  filter(condition!="NegCtr")

# Factorize relevant parameters
ENV$condition <- factor(ENV$condition, 
                        levels = c("CTR","ALG"))
ENV$site <- factor(ENV$site, 
                   levels = c("SouthPacific","EquatorialPacific","NorthPacific"))
ENV$type <- factor(ENV$type, 
                   levels = c("DNA","RNA"))
ENV$fraction <- factor(ENV$fraction, 
                       levels = c("FL","PA"))
ENV$day <- factor(ENV$day, 
                  levels = c("0","1","3","6"))
ENV$exp_cond <- factor(ENV$exp_cond, 
                       levels = c("E1_T0","E1_CTR","E1_ALG",
                                  "E2_T0","E2_CTR","E2_ALG",
                                  "E3_T0","E3_CTR","E3_ALG"))
ENV$fraction_TRUE <- factor(ENV$fraction_TRUE, 
                            levels = c("BULK","FL","PA"))


##### FORMAT TAXONOMY #############################
# Rename BAC-NAs with last known taxrank + "uc"
k <- ncol(TAX)-1
for (i in 2:k) {
  if (sum(is.na(TAX[, i])) >1) {
    temp <- TAX[is.na(TAX[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX[is.na(TAX[, i]), ] <- temp}
  if (sum(is.na(TAX[, i]))==1) {
    temp <- TAX[is.na(TAX[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX[is.na(TAX[, i]),] <- temp
  }
}
TAX[is.na(TAX[, (k+1)]), (k+1)] <- paste(
  TAX[is.na(TAX[, (k+1)]), k], " uc", sep="")

# shorten/modify names
TAX <- TAX %>%
  mutate(across(everything(),~gsub("Clade ","SAR11 Clade ", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("(SAR406 clade)","SAR406", ., fixed=T)))


##### PHYLOSEQ LOAD ###############################

asv = otu_table(ASV, taxa_are_rows=T)
tax = tax_table(TAX)
rownames(tax) <- rownames(asv)

pseq <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
  sample_data(ENV), 
  tax_table(tax))

# Fix tax-IDs
colnames(pseq@tax_table)<- c(
  "Kingdom","Phylum","Class","Order","Family","Genus","Species")

##### FILTRATION AND RAREFACTION ##################

## Filter: ASVs >3 counts in >3% of samples; to rel. ab
pseq.abs = filter_taxa(
  pseq, function(x) sum(x > 3) > (0.03 * length(x)), T)
pseq.rel = transform_sample_counts(
  pseq.abs, function(x) x / sum(x) * 100) 

## Rarefying 
pseq.rar = rarefy_even_depth(
  pseq.abs, sample.size = 10817, 
  rngseed = TRUE, replace = FALSE, 
  trimOTUs = TRUE, verbose = TRUE) 
# 2 samples removedbecause they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 
# MQ180801-_clipMQ180819-45_clip

# Filter: ASVs >3 counts in >3% of samples
pseq.abs.rar = filter_taxa(
  pseq.rar, function(x) sum(x > 0) > (0 * length(x)), T) 

pseq.rel.rar = transform_sample_counts(
  pseq.abs.rar, function(x) x / sum(x) * 100) 


##### AMPVIS LOAD #################################
ampvis <- data.frame(OTU = rownames(
  phyloseq::otu_table(pseq.abs)@.Data),
  phyloseq::otu_table(pseq.abs)@.Data,
  phyloseq::tax_table(pseq.abs)@.Data,
  check.names=F)

# Extract metadata from phyloseq; format and combine
ampvis.env <- data.frame(
  phyloseq::sample_data(pseq.abs), check.names=F)
ampvis.env <- cbind(
  Sample = rownames(ampvis.env), ampvis.env) 
ampvis <- amp_load(ampvis, ampvis.env)

### AMPVIS.RAR
ampvis.rar <- data.frame(OTU = rownames(
  phyloseq::otu_table(pseq.abs.rar)@.Data),
  phyloseq::otu_table(pseq.abs.rar)@.Data,
  phyloseq::tax_table(pseq.abs.rar)@.Data,
  check.names=F)

# Extract metadata from phyloseq; format and combine
ampvis.rar.env <- data.frame(
  phyloseq::sample_data(pseq.abs.rar), check.names=F)
ampvis.rar.env <- cbind(
  Sample = rownames(ampvis.rar.env), ampvis.rar.env) 
ampvis.rar <- amp_load(ampvis.rar, ampvis.rar.env)

##### REMOVAL OF OUTLIERS #########################
# removing "MQ180801-198_clip", "MQ180819-45_clip" from unrarefied dataset
ampvis <- amp_filter_samples(ampvis, !clip_id %in% c("MQ180801-198_clip", "MQ180819-45_clip"), normalise = FALSE)
# 2 samples and 0 OTUs have been filtered 
# Before: 132 samples and 3133 OTUs
# After: 130 samples and 3133 OTUs

ampvis.rar # Just to check how many samples and OTUs this dataset contains; matches unrarefied dataset
# 130 samples and 3133 OTUs



##### MAIN FIGURES ################################
##### Figure 1 ####################################
alginate <- read.csv("SO248_metadata.csv")

mean <- aggregate(alginate[,5:12], by = list(alginate$Site, alginate$Day, alginate$Exp, alginate$Condition), mean, na.rm=TRUE)
sd <- aggregate(alginate[,5:12], by = list(alginate$Site, alginate$Day, alginate$Exp, alginate$Condition), sd, na.rm=TRUE)

namekeymean <-  c("Group.1"="Site", "Group.2"="Day", "Group.3"="Exp", "Group.4"="Condition", "cells_mL"="cells", 
                  "BPP"="BPP", "Leucine_incorp"="Leu_incorp", "Protein"="Protein", "Doubling"="Doubling",
                  "Generation"="Generation", "Glc_turnover_d.1"="Glc_turnover_d.1", "Glc_turnover.rate_d.1"="Glc_turnover.rate_d.1")

names(mean) <- namekeymean[names(mean)]

namekeysd <-  c("Group.1"="Site", "Group.2"="Day", "Group.3"="Exp", "Group.4"="Condition", "cells_mL"="cells_sd", 
                "BPP"="BPP_sd", "Leucine_incorp"="Leu_incorp_sd", "Protein"="Protein_sd", "Doubling"="Doubling_sd",
                "Generation"="Generation_sd", "Glc_turnover_d.1"="Glc_turnover_d.1_sd", "Glc_turnover.rate_d.1"="Glc_turnover.rate_d.1_sd")

names(sd) <- namekeysd[names(sd)]

# Final alginate data frame
sd
sd_long <- gather(sd, factor_sd, sd, cells_sd:Glc_turnover.rate_d.1_sd, factor_key=TRUE)
sd <- sd_long[,6]

mean_long <- gather(mean, factor, mean, cells:Glc_turnover.rate_d.1, factor_key=TRUE)

alg <- cbind(mean_long,sd)

alg$Fraction = "Bulk"
alg$Exp <- gsub("Exp1", "E1", alg$Exp)
alg$Exp <- gsub("Exp2", "E2", alg$Exp)
alg$Exp <- gsub("Exp3", "E3", alg$Exp)
alg$Condition <- gsub("alginate", "ALG", alg$Condition)
alg$Condition <- gsub("control", "CTR", alg$Condition)

alg_sub <- subset(alg, factor != "Leu_incorp")
alg_sub <- subset(alg_sub, factor != "Doubling")
alg_sub <- subset(alg_sub, factor != "Generation")
alg_sub <- subset(alg_sub, factor != "Glc_turnover_d.1")
alg_sub <- subset(alg_sub, factor != "Protein")

alg_sub$factor <- gsub("cells", "Cell counts", alg_sub$factor)
alg_sub$factor <- gsub("BPP", "Bac. prot. prod.", alg_sub$factor)
alg_sub$factor <- gsub("Glc_turnover.rate_d.1", "Glu. turnover", alg_sub$factor)

####Figure 2 - Cell, BPP, Glu turnover, Alpha-div 
alg_meta = alg_sub

alg_meta$frac_cond = paste(alg_meta$Fraction, "_", alg_meta$Condition)
alg_meta$exp_cond =paste(alg_meta$Exp, alg_meta$Condition, sep="_")

alg_meta$factor <- factor(alg_meta$factor, levels = c("Cell counts","Bac. prot. prod.", "Glu. turnover", ordered=TRUE))

alg_meta_line = ggplot(alg_meta, aes(x=Day, y=mean, color=exp_cond, group=frac_cond, fill=exp_cond, shape=exp_cond), na.rm=TRUE) +
  geom_line(linetype="solid", size = 0.5, alpha = 0.4)  + 
  geom_point(size=2.5) + 
  scale_shape_manual(values=c(20, 1, 20, 1, 20, 1))+
  geom_errorbar(alpha = 0.4, aes(ymin=mean-sd, ymax=mean+sd), width = 0.4,  size = 0.5) +
  facet_grid(factor~Exp, scales = "free") +
  scale_color_manual(values=c(
    "E1_CTR"="#117733","E1_ALG"="#117733",
    "E2_CTR"="#882255","E2_ALG"="#882255",
    "E3_CTR"="#332288","E3_ALG"="#332288")) +
  scale_fill_manual(values=c(
    "E1_CTR"="#ffffff","E1_ALG"="#117733",
    "E2_CTR"="#ffffff","E2_ALG"="#882255",
    "E3_CTR"="#ffffff","E3_ALG"="#332288")) +
  theme_classic() + 
  labs(y = "") +
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=11, colour="black"),
        strip.text.y = element_text(size=11, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(alg_meta_line)
dev.copy(pdf,'Fig1B - alg_meta_new_2023.11.11.pdf', width = 4.5, height =4)
dev.off()

##### Figure 2 ####################################
# Figure 2A (contains 78 samples and 3127 ASVs)
PCoA_ALL_DNA = amp_ordinate(amp_subset_samples(
    ampvis, type=="DNA"),
  type = "PCoA",
  filter_species = 0.01,  # have adjusted this; please check no. of filtered ASVs
  distmeasure = "bray",
  transform = "none", 
  sample_color_by = "exp_cond", sample_colorframe = T,
  sample_shape_by = "fraction_TRUE", sample_label_size = 3,
  sample_point_size = 3) +
  scale_color_manual(values=c(
    "E1_T0"="#AAAAAA","E1_CTR"="#999933","E1_ALG"="#117733",
    "E2_T0"="#000000","E2_CTR"="#CC6677","E2_ALG"="#882255",
    "E3_T0"="#DDDDDD","E3_CTR"="#88CCEE","E3_ALG"="#332288")) +
  scale_fill_manual(values=c(
    "E1_T0"="#AAAAAA","E1_CTR"="#999933","E1_ALG"="#117733",
    "E2_T0"="#000000","E2_CTR"="#CC6677","E2_ALG"="#882255",
    "E3_T0"="#DDDDDD","E3_CTR"="#88CCEE","E3_ALG"="#332288")) +
  theme_classic() + 
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(PCoA_ALL_DNA)
dev.copy(pdf,'Fig2A - PCoA_ALL_DNA_2023.11.11.pdf', width = 4, height = 3.6)
dev.off()

# Figure 2B - Panel 1 (E1 - South)
PCoA_E1_DNA = amp_ordinate(
  amp_subset_samples(
    ampvis, site=="SouthPacific" & type=="DNA"),
  type = "PCoA",
  distmeasure = "bray",
  transform = "none",
  sample_color_by = "exp_cond_short", sample_colorframe = T,
  sample_shape_by = "fraction_TRUE", sample_label_size = 3,
  sample_point_size = 3) +
  scale_color_manual(values=c(
    "T0"="#AAAAAA","CTR"="#999933","ALG"="#117733")) +
  scale_fill_manual(values=c(
    "T0"="#AAAAAA","CTR"="#999933","ALG"="#117733")) +
  theme_classic() + 
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(PCoA_E1_DNA)
dev.copy(pdf,'Fig2B - PCoA_E1_DNA_2023.11.11.pdf', width = 2.5, height = 2.2)
dev.off()

# Figure 2B - Panel 2 (E2 - Equator)
PCoA_E2_DNA = amp_ordinate(
  amp_subset_samples(
    ampvis, site=="EquatorialPacific" & type=="DNA"),
  type = "PCoA",
  distmeasure = "bray",
  transform = "none", 
  sample_color_by = "exp_cond_short", sample_colorframe = T,
  sample_shape_by = "fraction_TRUE", sample_label_size = 3,
  sample_point_size = 3) +
  scale_color_manual(values=c(
    "T0"="#000000","CTR"="#CC6677","ALG"="#882255")) +
  scale_fill_manual(values=c(
    "T0"="#000000","CTR"="#CC6677","ALG"="#882255")) +
  theme_classic() + 
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(PCoA_E2_DNA)
dev.copy(pdf,'Fig2B - PCoA_E2_DNA_2023.11.11.pdf', width = 2.5, height = 2.2)
dev.off()


# Figure 2B - Panel 3 (E3 - North)
PCoA_E3_DNA = amp_ordinate(
  amp_subset_samples(
    ampvis, site == "NorthPacific" & type=="DNA"),
  type = "PCoA",
  distmeasure = "bray",
  transform = "none", 
  sample_color_by = "exp_cond_short", sample_colorframe = T,
  sample_shape_by = "fraction_TRUE", sample_label_size = 3,
  sample_point_size = 3) +
  scale_color_manual(values=c(
    "T0"="#DDDDDD","CTR"="#88CCEE","ALG"="#332288")) +
  scale_fill_manual(values=c(
    "T0"="#DDDDDD","CTR"="#88CCEE","ALG"="#332288")) +
  theme_classic() + 
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(PCoA_E3_DNA)
dev.copy(pdf,'Fig2B - PCoA_E3_DNA_2023.11.11.pdf', width = 2.5, height = 2.2)
dev.off()

### Calculating alpha diversity indices for Figure 2C
alpha_div = amp_alphadiv(amp_subset_samples(ampvis.rar, type=="DNA"), 
  measure = NULL, richness = FALSE)
#51 samples and 54 OTUs have been filtered 
#Before: 129 samples and 3146 OTUs
#After: 78 samples and 3092 OTUs

alpha_div_sub = alpha_div %>% select(site,condition,exp,exp_cond, day,fraction_TRUE,
                                     ObservedOTUs,Shannon,Simpson,invSimpson)

alpha_div_long <- gather(alpha_div_sub, richness, value, ObservedOTUs:invSimpson, factor_key=TRUE)

# Aggregate to calculate means and 
alpha_div_long_mean <- aggregate(alpha_div_long[,8], 
                                 by = list(alpha_div_long$site, alpha_div_long$condition, 
                                           alpha_div_long$exp, 
                                           alpha_div_long$exp_cond, 
                                           alpha_div_long$day, 
                                           alpha_div_long$fraction_TRUE,
                                           alpha_div_long$richness), mean, na.rm=TRUE)
namekeymean <-  c("Group.1"="site", "Group.2"="condition", "Group.3"="exp", "Group.4"="exp_cond", 
                  "Group.5" = "day", "Group.6" = "fraction_TRUE", "Group.7" = "richness", "x"="mean")
names(alpha_div_long_mean) <- namekeymean[names(alpha_div_long_mean)]

alpha_div_long_sd <- aggregate(alpha_div_long[,8], 
                               by = list(alpha_div_long$site, alpha_div_long$condition, 
                                         alpha_div_long$exp, 
                                         alpha_div_long$exp_cond, 
                                         alpha_div_long$day, 
                                         alpha_div_long$fraction_TRUE,
                                         alpha_div_long$richness), sd, na.rm=TRUE)
sd <- alpha_div_long_sd[,8]
namekeysd <-  c("x"="sd")
names(sd) <- namekeysd[names(sd)]

alpha_df <- cbind(alpha_div_long_mean, sd)

alpha_df$replace = alpha_df$condition
alpha_df$replace <- gsub("T0", "CTR", alpha_df$replace)
alpha_df$frac_cond = paste(alpha_df$fraction_TRUE,alpha_df$replace, sep="_")

# Shannon Diversity: Figure 2C Panels 1-3
Shannon_div = subset(alpha_df, richness == "Shannon")

#Shannon_div$exp_cond = paste(Shannon_div$exp, Shannon_div$condition, sep="_")

#alpha_sub$Fraction <- gsub("BULK_CTR", "Bulk", alpha_sub$Fraction)

shannon_plot <- ggplot(Shannon_div, aes(x=day, y=mean, color = exp_cond, group = frac_cond))+
  geom_errorbar(alpha = 0.4, aes(ymin=mean-sd, ymax=mean+sd), width = 0.2,  size = 0.5) +
  geom_line(linetype="solid", size = 0.5, alpha = 0.4, aes(group=frac_cond)) + 
  geom_point(size=3, color = "black", aes(shape = frac_cond)) +
  geom_point(size=2.5, aes(shape = frac_cond, color = exp_cond))+ 
  facet_grid(.~exp,scales="free_y") +
  theme_classic() + 
  ylab(expression(Shannon~diversity)) +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values=c(
    "E1_T0"="#AAAAAA","E1_CTR"="#999933","E1_ALG"="#117733",
    "E2_T0"="#000000","E2_CTR"="#CC6677","E2_ALG"="#882255",
    "E3_T0"="#DDDDDD","E3_CTR"="#88CCEE","E3_ALG"="#332288")) +
  scale_fill_manual(values=c(
    "E1_T0"="#AAAAAA","E1_CTR"="#999933","E1_ALG"="#117733",
    "E2_T0"="#000000","E2_CTR"="#CC6677","E2_ALG"="#882255",
    "E3_T0"="#DDDDDD","E3_CTR"="#88CCEE","E3_ALG"="#332288")) +
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(shannon_plot)
dev.copy(pdf,'Fig2C - shannon_line_2023.11.11.pdf', width = 5, height =2.2)
dev.off()
##### Figure 3 ####################################
# contains 78 samples and 3127 ASVs
top20_families = amp_heatmap(amp_subset_samples(ampvis, type=="DNA"),
  tax_aggregate = "Family",
  tax_add = "Class",
  showRemainingTaxa = TRUE,
  tax_show = 20,
  normalise = T,
  plot_values = F,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 15,
  plot_colorscale = "sqrt",
  facet_by = c("site","fraction_TRUE","condition"),
  group_by = "day",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, colour = "black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=9, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(top20_families)
dev.copy(pdf,'Fig3 - top20_families_heatmap_2023.11.11.pdf', width = 9, height = 3.5)
dev.off()


##### Figure 4A ##################################
# (DESeq was conducted on the rarefied dataset)

## Step 1: Run DESeq analyses (PA vs CTR) combined days, experiments separate #
# Total # of upreg OTUs identified: 161 (E1) + 61 (E2) + 82 (E3) = 304 (nrow(rar_PA_upreg))
# Total # of UNIQUE OTUs identified: 284 (nrow(onlyPA_rar_upreg_merge2); onlyPA_rar_upreg_combined_physeq)

# rar E1 d1, d3, d6 combined (161 ASvs upreg, log2fold >=3)
onlyPA_rarE1_DNA <- subset_samples(pseq.abs.rar, exp == "E1" & type=="DNA" & day != "0" & fraction_TRUE != "FL")
onlyPA_rarE1_DNA_DESeq = phyloseq_to_deseq2(onlyPA_rarE1_DNA, ~ fraction_TRUE) # Convert to DESeq2 file
onlyPA_rarE1_DNA_DESeq = DESeq(onlyPA_rarE1_DNA_DESeq, test="Wald", fitType="parametric") # Run test
res = results(onlyPA_rarE1_DNA_DESeq, cooksCutoff = FALSE) # Save Result
alpha = 0.01 #Produce table
sigtab_onlyPA_rarE1 = res[which(res$padj < alpha), ]
sigtab_onlyPA_rarE1 = cbind(as(sigtab_onlyPA_rarE1, "data.frame"), as(tax_table(pseq.abs)[rownames(sigtab_onlyPA_rarE1), ], "matrix"))
sigtab_onlyPA_rarE1_up = sigtab_onlyPA_rarE1 %>% dplyr::filter(log2FoldChange >= 3) # To deconvolute, keep only ASVs with log2FC vals of >3
nrow(sigtab_onlyPA_rarE1_up)

# rar E2 d1,d3,d6 combined (61 ASvs upreg, log2fold >=3)
onlyPA_rarE2_DNA <- subset_samples(pseq.abs.rar, exp == "E2" & type=="DNA" & day != "0" & fraction_TRUE != "FL")
onlyPA_rarE2_DNA_DESeq = phyloseq_to_deseq2(onlyPA_rarE2_DNA, ~ fraction_TRUE) # Convert to DESeq2 file
onlyPA_rarE2_DNA_DESeq = DESeq(onlyPA_rarE2_DNA_DESeq, test="Wald", fitType="parametric") # Run test
res = results(onlyPA_rarE2_DNA_DESeq, cooksCutoff = FALSE) # Save Result
alpha = 0.01 #Produce table
sigtab_onlyPA_rarE2 = res[which(res$padj < alpha), ]
sigtab_onlyPA_rarE2 = cbind(as(sigtab_onlyPA_rarE2, "data.frame"), as(tax_table(pseq.abs)[rownames(sigtab_onlyPA_rarE2), ], "matrix"))
sigtab_onlyPA_rarE2_up = sigtab_onlyPA_rarE2 %>% dplyr::filter(log2FoldChange >= 3) # To deconvolute, keep only ASVs with log2FC vals of >3
nrow(sigtab_onlyPA_rarE2_up)

# rar E3 d1&d3,d6 combined (82 ASvs upreg, log2fold >=3) 
onlyPA_rarE3_DNA <- subset_samples(pseq.abs.rar, exp == "E3" & type=="DNA" & day != "0" & fraction_TRUE != "FL")
onlyPA_rarE3_DNA_DESeq = phyloseq_to_deseq2(onlyPA_rarE3_DNA, ~ fraction_TRUE) # Convert to DESeq2 file
onlyPA_rarE3_DNA_DESeq = DESeq(onlyPA_rarE3_DNA_DESeq, test="Wald", fitType="parametric") # Run test
res = results(onlyPA_rarE3_DNA_DESeq, cooksCutoff = FALSE) # Save Result
alpha = 0.01 #Produce table
sigtab_onlyPA_rarE3 = res[which(res$padj < alpha), ]
sigtab_onlyPA_rarE3 = cbind(as(sigtab_onlyPA_rarE3, "data.frame"), as(tax_table(pseq.abs)[rownames(sigtab_onlyPA_rarE3), ], "matrix"))
sigtab_onlyPA_rarE3_up = sigtab_onlyPA_rarE3 %>% dplyr::filter(log2FoldChange >= 3) # To deconvolute, keep only ASVs with log2FC vals of >3
nrow(sigtab_onlyPA_rarE3_up)

## Step 2: Combine all upreg ASVs into one file and subset phyloseq object
sigtab_onlyPA_rarE1_up$rownames =row.names(sigtab_onlyPA_rarE1_up)
sigtab_onlyPA_rarE2_up$rownames =row.names(sigtab_onlyPA_rarE2_up)
sigtab_onlyPA_rarE3_up$rownames =row.names(sigtab_onlyPA_rarE3_up)

onlyPA_rar_upreg_merge1 <- merge(sigtab_onlyPA_rarE1_up, sigtab_onlyPA_rarE2_up, by="rownames", all=TRUE)
onlyPA_rar_upreg_merge2 <- merge(onlyPA_rar_upreg_merge1, sigtab_onlyPA_rarE3_up, by="rownames", all=TRUE)

onlyPA_rar_upreg_combined_vector = as.vector((onlyPA_rar_upreg_merge2$rownames))

onlyPA_rar_upreg_combined_subset <- subset(otu_table(pseq.abs), rownames(otu_table(pseq.abs)) %in% onlyPA_rar_upreg_combined_vector)

onlyPA_rar_upreg_combined_physeq <- merge_phyloseq(onlyPA_rar_upreg_combined_subset, tax_table(pseq.abs), sample_data(pseq.abs))
# contained 284 taxa (284 taxa no matter the use of pseq.abs.rar or pseq.abs)

## Figure 4A - Classes of PA ASVs thta are differentially enriched in alginate treatments versus controls
sigtab_onlyPA_rarE1_up$Experiment = "E1"
sigtab_onlyPA_rarE2_up$Experiment = "E2"
sigtab_onlyPA_rarE3_up$Experiment = "E3"

rar_PA_upreg = rbind(sigtab_onlyPA_rarE1_up, sigtab_onlyPA_rarE2_up, sigtab_onlyPA_rarE3_up)
# 304

rar_PA_upreg_lim20 = rar_PA_upreg %>% dplyr::filter(log2FoldChange < 20)
# removes 2 ASVs from the visualization; for visual simplicity

ggplot(rar_PA_upreg_lim20, aes(x=log2FoldChange, y=reorder(Class, log2FoldChange), color=Experiment)) +
  geom_point(alpha = 0.2) +  xlim(2,11.5) +
  geom_boxplot(aes(fill=Experiment), alpha = 0.1, size = 0.2) + 
  #geom_point(size=1.5, alpha = 0.3) + 
  facet_grid(.~Experiment) +
  ylab(expression()) +
  #scale_x_continuous(position = 'top') +
  
  scale_color_manual(values=c("E1"="#117733", "E2"="#882255", "E3"="#332288")) +
  scale_fill_manual(values=c("E1"="#117733", "E2"="#882255", "E3"="#332288")) +
  theme_classic() +
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
dev.copy(pdf,'Fig4A - rar_PA_upreg_lim20_box_2023.11.11.pdf', width = 4.5, height = 4)
dev.off()


##### Figure 4B ###################################

## Step 1: Run DESeq analyses (FL vs CTR) combined days, experiment separately 
# Total # of upreg OTUs identified: 25 (E1) + 32 (E2) + 24 (E3) = 81 (nrow(rar_FL_upreg))
# Total # of UNIQUE OTUs identified: 76 (nrow(onlyFL_rar_upreg_merge2); onlyFL_rar_upreg_combined_physeq)

# rar E1 d1&d3 combined (25 ASvs upreg, log2fold >=3)
onlyFL_rarE1_DNA <- subset_samples(pseq.abs.rar, exp == "E1" & type=="DNA" & day != "0" & fraction_TRUE != "PA")
onlyFL_rarE1_DNA_DESeq = phyloseq_to_deseq2(onlyFL_rarE1_DNA, ~ fraction_TRUE) # Convert to DESeq2 file
onlyFL_rarE1_DNA_DESeq = DESeq(onlyFL_rarE1_DNA_DESeq, test="Wald", fitType="parametric") # Run test
res = results(onlyFL_rarE1_DNA_DESeq, cooksCutoff = FALSE) # Save Result
alpha = 0.01 #Produce table
sigtab_onlyFL_rarE1 = res[which(res$padj < alpha), ]
sigtab_onlyFL_rarE1 = cbind(as(sigtab_onlyFL_rarE1, "data.frame"), as(tax_table(pseq.abs)[rownames(sigtab_onlyFL_rarE1), ], "matrix"))
sigtab_onlyFL_rarE1_up = sigtab_onlyFL_rarE1 %>% dplyr::filter(log2FoldChange >= 3) # To deconvolute, keep only ASVs with log2FC vals of >3
nrow(sigtab_onlyFL_rarE1_up)

# rar E2 d1, d3, d6 combined (22 ASvs upreg, log2fold >=3)
onlyFL_rarE2_DNA <- subset_samples(pseq.abs.rar, exp == "E2" & type=="DNA" & day != "0" & fraction_TRUE != "PA")
onlyFL_rarE2_DNA_DESeq = phyloseq_to_deseq2(onlyFL_rarE2_DNA, ~ fraction_TRUE) # Convert to DESeq2 file
onlyFL_rarE2_DNA_DESeq = DESeq(onlyFL_rarE2_DNA_DESeq, test="Wald", fitType="parametric") # Run test
res = results(onlyFL_rarE2_DNA_DESeq, cooksCutoff = FALSE) # Save Result
alpha = 0.01 #Produce table
sigtab_onlyFL_rarE2 = res[which(res$padj < alpha), ]
sigtab_onlyFL_rarE2 = cbind(as(sigtab_onlyFL_rarE2, "data.frame"), as(tax_table(pseq.abs)[rownames(sigtab_onlyFL_rarE2), ], "matrix"))
sigtab_onlyFL_rarE2_up = sigtab_onlyFL_rarE2 %>% dplyr::filter(log2FoldChange >= 3) # To deconvolute, keep only ASVs with log2FC vals of >3
nrow(sigtab_onlyFL_rarE2_up)

# rar E3 d1, d3, d6 combined (10 ASvs upreg, log2fold >=3) 
onlyFL_rarE3_DNA <- subset_samples(pseq.abs.rar, exp == "E3" & type=="DNA" & day != "0" & fraction_TRUE != "PA")
onlyFL_rarE3_DNA_DESeq = phyloseq_to_deseq2(onlyFL_rarE3_DNA, ~ fraction_TRUE) # Convert to DESeq2 file
onlyFL_rarE3_DNA_DESeq = DESeq(onlyFL_rarE3_DNA_DESeq, test="Wald", fitType="parametric") # Run test
res = results(onlyFL_rarE3_DNA_DESeq, cooksCutoff = FALSE) # Save Result
alpha = 0.01 #Produce table
sigtab_onlyFL_rarE3 = res[which(res$padj < alpha), ]
sigtab_onlyFL_rarE3 = cbind(as(sigtab_onlyFL_rarE3, "data.frame"), as(tax_table(pseq.abs)[rownames(sigtab_onlyFL_rarE3), ], "matrix"))
sigtab_onlyFL_rarE3_up = sigtab_onlyFL_rarE3 %>% dplyr::filter(log2FoldChange >= 3) # To deconvolute, keep only ASVs with log2FC vals of >3
nrow(sigtab_onlyFL_rarE3_up)

## Step 2: Combine all upreg ASVs into one file and subset phyloseq object
sigtab_onlyFL_rarE1_up$rownames =row.names(sigtab_onlyFL_rarE1_up)
sigtab_onlyFL_rarE2_up$rownames =row.names(sigtab_onlyFL_rarE2_up)
sigtab_onlyFL_rarE3_up$rownames =row.names(sigtab_onlyFL_rarE3_up)

onlyFL_rar_upreg_merge1 <- merge(sigtab_onlyFL_rarE1_up, sigtab_onlyFL_rarE2_up, by="rownames", all=TRUE)
onlyFL_rar_upreg_merge2 <- merge(onlyFL_rar_upreg_merge1, sigtab_onlyFL_rarE3_up, by="rownames", all=TRUE)

onlyFL_rar_upreg_combined_vector = as.vector((onlyFL_rar_upreg_merge2$rownames))

onlyFL_rar_upreg_combined_subset <- subset(otu_table(pseq.abs), rownames(otu_table(pseq.abs)) %in% onlyFL_rar_upreg_combined_vector)

onlyFL_rar_upreg_combined_physeq <- merge_phyloseq(onlyFL_rar_upreg_combined_subset, tax_table(pseq.abs), sample_data(pseq.abs))
# contained 76 taxa

## Figure 4B - Classes of PA ASVs thta are differentially enriched in alginate treatments versus controls
sigtab_onlyFL_rarE1_up$Experiment = "E1"
sigtab_onlyFL_rarE2_up$Experiment = "E2"
sigtab_onlyFL_rarE3_up$Experiment = "E3"

rar_FL_upreg = rbind(sigtab_onlyFL_rarE1_up, sigtab_onlyFL_rarE2_up, sigtab_onlyFL_rarE3_up)
# nrow(rar_FL_upreg) 79 ASVs

rar_FL_upreg_lim20 = rar_FL_upreg %>% dplyr::filter(log2FoldChange < 20)

ggplot(rar_FL_upreg_lim20, aes(x=log2FoldChange, y=reorder(Class, log2FoldChange), color=Experiment)) +
  geom_point(alpha = 0.2) + xlim(2,11.5) +
  geom_boxplot(aes(fill=Experiment), alpha = 0.1, size = 0.2) + 
  #geom_point(size=1.5, alpha = 0.3) + 
  facet_grid(.~Experiment) +
  ylab(expression()) +
  #scale_x_continuous(position = 'top') +
  scale_color_manual(values=c("E1"="#117733", "E2"="#882255", "E3"="#332288")) +
  scale_fill_manual(values=c("E1"="#117733", "E2"="#882255", "E3"="#332288")) +
  theme_classic() +
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
dev.copy(pdf,'Fig4B - rar_FL_upreg_lim20_box_2023.11.11.pdf', width = 4.5, height = 2.1)
dev.off()

##### Figure 4C ###################################
# Venn diagram of shared PA and FL across sites

# For rar_PA_upreg_combined (n = 284)
pseq.abs.DNA = subset_samples(pseq.abs, type =="DNA" & clip_id != "MQ180801-198_clip" & clip_id != "MQ180819-45_clip")

rar_PA_venn_subset <- subset(otu_table(pseq.abs.DNA), rownames(otu_table(pseq.abs.DNA)) %in% onlyPA_rar_upreg_merge2$rownames)
rar_PA_subset <- merge_phyloseq(rar_PA_venn_subset, tax_table(pseq.abs.DNA), sample_data(pseq.abs.DNA))
rar_PA_venn <- ps_venn(rar_PA_subset, group="exp", quantities = list(type=c("percent","counts"), font = 2), 
                       labels = list(cex = 2), col = "black", fill = c("#DDECE2","#EDDEE6","#E2DFED"))
print(rar_PA_venn)
dev.copy(pdf,'Fig 4C - rar_PA_venn_in_unrar_pseq_2023.11.11.pdf', width = 4, height = 4)
dev.off()

# For rar_FL_upreg_combined
rar_FL_venn_subset <- subset(otu_table(pseq.abs.DNA), rownames(otu_table(pseq.abs.DNA)) %in% onlyFL_rar_upreg_merge2$rownames)
rar_FL_subset <- merge_phyloseq(rar_FL_venn_subset, tax_table(pseq.abs.DNA), sample_data(pseq.abs.DNA))
rar_FL_venn <- ps_venn(rar_FL_subset, group="exp", quantities = list(type=c("percent","counts"), font = 2), 
                       labels = list(cex = 2), col = "black", fill = c("#DDECE2","#EDDEE6","#E2DFED"))
print(rar_FL_venn)
dev.copy(pdf,'Fig4C - rar_FL_venn_in_unrar_pseq_223.11.11.pdf', width = 4, height = 4)
dev.off()



##### Figure 5 ###################################

# Fig. 5A - Particle-associated, top 10
only_PA_rar_upreg_combined_ASVs = amp_subset_taxa(ampvis, tax_vector = onlyPA_rar_upreg_combined_vector, normalise = TRUE, remove = FALSE)
#284 OTUs

only_PA_rar_upreg_combined_heatmap = amp_heatmap(amp_subset_samples(
  only_PA_rar_upreg_combined_ASVs, type=="DNA" & fraction_TRUE =="PA" ),
  tax_aggregate = "OTU",
  tax_add = "Genus",
  showRemainingTaxa = FALSE,
  tax_show = 10,
  normalise = F,
  plot_values = F,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 10,
  plot_colorscale = "sqrt",
  facet_by = c("site"),
  group_by = "day",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, colour="black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=10, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(only_PA_rar_upreg_combined_heatmap)
dev.copy(pdf,'Fig5A - top10_only_PA_rar_upreg_combined_heatmap_2023.11.12.pdf', width = 5, height = 2.5)
dev.off()

# Fig. 5B - Free-living, top 10
only_FL_rar_upreg_combined_ASVs = amp_subset_taxa(ampvis, tax_vector = onlyFL_rar_upreg_combined_vector, normalise = TRUE, remove = FALSE)
#76 OTUs

only_FL_rar_upreg_combined_heatmap = amp_heatmap(amp_subset_samples(
  only_FL_rar_upreg_combined_ASVs, type=="DNA" & fraction_TRUE =="FL" ),
  tax_aggregate = "OTU",
  tax_add = "Genus",
  showRemainingTaxa = FALSE,
  tax_show = 10,
  normalise = F,
  plot_values = F,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 10,
  plot_colorscale = "sqrt",
  facet_by = c("site"),
  group_by = "day",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, colour="black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=10, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(only_FL_rar_upreg_combined_heatmap)
dev.copy(pdf,'Fig5B - top10_only_FL_rar_upreg_combined_heatmap_2023.11.12.pdf', width = 4.5, height = 2.5)
dev.off()

## Fig. 5C - Top 10 diff abundant ASVs shared between PA and FL

# Step 1: identify the top 10 diff abundant ASVs shared between PA and FL
common_PAFL_rar_upreg_combined.df = semi_join(onlyPA_rar_upreg_merge2,onlyFL_rar_upreg_merge2, by = "rownames") 
nrow(common_PAFL_rar_upreg_combined.df)
# 53 ASVs
common_PAFL_rar_upreg_combined_vector = as.vector(common_PAFL_rar_upreg_combined.df$rownames)
common_PAFL_rar_upreg_combined_ASVs = amp_subset_taxa(ampvis, tax_vector = common_PAFL_rar_upreg_combined_vector, normalise = TRUE, remove = FALSE)

common_PAFL_rar_upreg_combined_heatmap = amp_heatmap(amp_subset_samples(
  common_PAFL_rar_upreg_combined_ASVs, type=="DNA"),
  tax_aggregate = "OTU",
  tax_add = "Genus",
  showRemainingTaxa = TRUE,
  tax_show = 10,
  normalise = F,
  plot_values = T,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 10,
  plot_colorscale = "sqrt",
  facet_by = c("site", "fraction_TRUE","condition"),
  group_by = "day",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, colour="black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=10, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "none") 
print(common_PAFL_rar_upreg_combined_heatmap)

# Create vector of top 10 ASVs identified from heatmap above
common_PAFL_rar_upreg_combined_ASVs_top10_vector = c("asv13", "asv33", "asv42", "asv7", "asv53", 
                                                     "asv41", "asv5", "asv10", "asv73", "asv78")
alg_asvs = as.data.frame(otu_table(pseq.rel))
alg_asvs_tr = as.data.frame(t(alg_asvs))
top10_ASVs_df = subset(alg_asvs_tr, select=common_PAFL_rar_upreg_combined_ASVs_top10_vector)
top10_ASVs_df$Sample = row.names(top10_ASVs_df)


env_subset_vector = c("Sample", "site", "exp", "condition", "exp_cond", "day", "fraction_TRUE", "fraction" , "type")
env_subset_df = subset(ampvis.env, select = env_subset_vector)

top10_df = merge(top10_ASVs_df,env_subset_df, by = "Sample") 

top10_ASVs_wide = melt(top10_df, id.vars=c("Sample", "site", "exp", "condition", "exp_cond", "day",
                                           "fraction_TRUE", "fraction", "type"))

top10_mean <- aggregate(top10_ASVs_wide$value, 
                        by = list(top10_ASVs_wide$site, top10_ASVs_wide$exp, 
                                  top10_ASVs_wide$condition, top10_ASVs_wide$exp_cond, top10_ASVs_wide$day,
                                  top10_ASVs_wide$fraction_TRUE, top10_ASVs_wide$fraction, top10_ASVs_wide$type, top10_ASVs_wide$variable), mean, na.rm=TRUE)
namekeymean <-  c("Group.1"="site", "Group.2"="exp", "Group.3"="condition", "Group.4"="exp_cond", 
                  "Group.5" = "day", "Group.6" = "fraction_TRUE", "Group.7" = "fraction", "Group.8" = "type", "Group.9" = "ASV", "x"="Mean")
names(top10_mean) <- namekeymean[names(top10_mean)]

top10_sd <- aggregate(top10_ASVs_wide$value, 
                      by = list(top10_ASVs_wide$site, top10_ASVs_wide$exp, 
                                top10_ASVs_wide$condition, top10_ASVs_wide$exp_cond, top10_ASVs_wide$day,
                                top10_ASVs_wide$fraction_TRUE, top10_ASVs_wide$fraction, top10_ASVs_wide$type, top10_ASVs_wide$variable), sd, na.rm=TRUE)
namekeysd <-  c("Group.1"="site", "Group.2"="exp", "Group.3"="condition", "Group.4"="exp_cond", 
                "Group.5" = "day", "Group.6" = "fraction_TRUE", "Group.7" = "fraction", "Group.8" = "type", "Group.9" = "ASV", "x"="SD")
names(top10_sd) <- namekeysd[names(top10_sd)]

top10_mean$SD = top10_sd$SD

top10_subset = subset(top10_mean, day != "0" & type == "DNA" & fraction_TRUE != "BULK")

# Plot Figure 5C
top10_subset$exp_fraction <- paste(top10_subset$exp, top10_subset$fraction_TRUE)

top10_subset_line = ggplot(top10_subset, aes(x=day, y=Mean, color = exp_fraction, group = exp_fraction), na.rm=TRUE) +
  geom_errorbar(alpha = 0.4, aes(ymin=Mean-SD, ymax=Mean+SD), width = 0.2,  size = 0.3) +
  geom_line(size = 1, alpha = 0.4, aes(linetype = fraction_TRUE)) + 
  geom_point(size=2.5, color = "black", aes(shape = fraction_TRUE)) + geom_point(size=2, aes(color=exp_fraction, shape = fraction_TRUE)) + 
  scale_shape_manual(values = c(17, 15)) +
  facet_wrap(~ASV, ncol = 5, scales = "free") +
  scale_color_manual(values = c("E1 FL" = "#117733", "E1 PA" = "#117733",
                                "E2 FL" = "#882255", "E2 PA" = "#882255",
                                "E3 FL" = "#332288", "E3 PA" = "#332288")) +
  theme_classic() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "bottom")
print(top10_subset_line)
dev.copy(pdf,'Fig5C - top10_subset_line_rar_2023.11.11.pdf', width = 7, height =4.5)
dev.off()

##### Figure 6 ####################################

EEA = read.csv("alg_EEA.csv")
EEA = subset(EEA, Condition != "T0")

# MCAMUF
MCAMUF = subset(EEA, Assay == "MCAMUF")


MCAMUF$Substrate <- factor(MCAMUF$Substrate, levels = c("a-glu","b-glu", "Leu", "AAF-chym", "AAPF-chym",
                                                        "QAR-tryp", "FSR-tryp", ordered=TRUE)) 

MCAMUF$Exp <- gsub("E1", "South", MCAMUF$Exp)
MCAMUF$Exp <- gsub("E2", "Equatorial", MCAMUF$Exp)
MCAMUF$Exp <- gsub("E3", "North", MCAMUF$Exp)

MCAMUF$Condition <- gsub("CTR", "control", MCAMUF$Condition)
MCAMUF$Condition <- gsub("ALG", "alginate", MCAMUF$Condition)
MCAMUF$Condition <- gsub("fold_change", "fold change", MCAMUF$Condition)

MCAMUF$Condition <- factor(MCAMUF$Condition, levels = c("control", "alginate", "fold change", ordered=TRUE))
MCAMUF$Exp <- factor(MCAMUF$Exp, levels = c("South","Equatorial", "North", ordered=TRUE))


MCAMUF_point = ggplot(MCAMUF, aes(x = Substrate, y=Mean, fill=Substrate), na.rm=TRUE) +   
  geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), width = 0.5,  size =0.2) +
  geom_bar(stat = "identity" , position = "dodge", width=1, colour = "#000000", size = 0.2) +
  facet_grid(Condition~Exp, scales = "free") +
  theme_classic() + 
  scale_color_manual(values = c("Leu" = "#999999","AAF-chym" = "#ff8000","AAPF-chym" = "#ff0000",
                                "QAR-tryp" = "#06E00B","FSR-tryp" = "#037106", 
                                "a-glu" = "#09C4DD", "b-glu" = "#751aff")) +
  scale_fill_manual(values = c("Leu" = "#999999","AAF-chym" = "#ff8000","AAPF-chym" = "#ff0000",
                               "QAR-tryp" = "#06E00B","FSR-tryp" = "#037106", 
                               "a-glu" = "#09C4DD", "b-glu" = "#751aff")) +
  ylab(expression(Hydrolysis~rates~(nmol~L^-1~hr^-1))) +
  xlab(expression()) +
  theme(axis.text.x  = element_text(size=9, colour="black", angle=45, hjust = 1),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=11, colour="black"),
        strip.text.y = element_text(size=11, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(MCAMUF_point)

# FLA
FLA = subset(EEA, Assay == "FLA")
FLA = subset(FLA, Substrate != "Fuc")
FLA = subset(FLA, Substrate != "Ara")

FLA$Condition <- gsub("CTR", "control", FLA$Condition)
FLA$Condition <- gsub("ALG", "alginate", FLA$Condition)
FLA$Condition <- gsub("fold_change", "fold change", FLA$Condition)

FLA$Substrate <- gsub("Chn", "chondroitin", FLA$Substrate)
FLA$Substrate <- gsub("Lam", "laminarin", FLA$Substrate)
FLA$Substrate <- gsub("Pul", "pullulan", FLA$Substrate)
FLA$Substrate <- gsub("Xyl", "xylan", FLA$Substrate)

FLA$Exp <- gsub("E1", "South", FLA$Exp)
FLA$Exp <- gsub("E2", "Equatorial", FLA$Exp)
FLA$Exp <- gsub("E3", "North", FLA$Exp)


FLA$Condition <- gsub("CTR", "control", FLA$Condition)
FLA$Condition <- gsub("ALG", "alginate", FLA$Condition)

FLA$Condition <- factor(FLA$Condition, levels = c("control", "alginate", "fold change", ordered=TRUE))
FLA$Exp <- factor(FLA$Exp, levels = c("South","Equatorial", "North", ordered=TRUE))
FLA$Substrate <- factor(FLA$Substrate, levels = c("chondroitin","laminarin", "pullulan", "xylan", ordered=TRUE)) 

FLA_point = ggplot(FLA, aes(x = Substrate, y=Mean, fill=Substrate), na.rm=TRUE) +   
  geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), width = 0.5,  size =0.2) +
  geom_bar(stat = "identity" , position = "dodge", width=1, colour = "#000000", size = 0.2) +
  facet_grid(Condition~Exp) +
  theme_classic() + 
  scale_color_manual(values = c("pullulan" = "#0000FF","laminarin" = "#FFFF30","xylan" = "#E90000",
                                "Fuc" = "#007E00","Ara" = "#cccccc","chondroitin" = "#4DF1F5")) +
  scale_fill_manual(values = c("pullulan" = "#0000FF","laminarin" = "#FFFF30","xylan" = "#E90000",
                               "Fuc" = "#007E00","Ara" = "#cccccc","chondroitin" = "#4DF1F5")) +
  xlab(expression()) +
  ylab(expression(Hydrolysis~rates~(nmol~monomer~L^-1~hr^-1))) +
  theme(axis.text.x  = element_text(size=9, colour="black", angle=45, hjust = 1),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=11, colour="black"),
        strip.text.y = element_text(size=11, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(FLA_point)


# Put it all together!
ggarrange(MCAMUF_point, FLA_point, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          legend = "none",
          align = "hv",
          widths = c(1.4,1))
dev.copy(pdf,'Fig6 - EEA.pdf', width = 7, height = 4)
dev.off()


##### SUPPLEMENTARY FIGURES #######################
##### Figure S1 ##################################

# Figure S1A - DNA and RNA colored by experimental treatment
PCoA_ALL_DNA_v_RNA_exp_cond_colored = amp_ordinate(
  amp_subset_samples(ampvis),
  type = "PCoA",
  distmeasure = "bray",
  transform = "none", 
  sample_color_by = "exp_cond", sample_colorframe = T,
  sample_shape_by = "fraction_TRUE", sample_label_size = 3,
  sample_point_size = 3) +
  scale_color_manual(values=c(
    "E1_T0"="#AAAAAA","E1_CTR"="#999933","E1_ALG"="#117733",
    "E2_T0"="#000000","E2_CTR"="#CC6677","E2_ALG"="#882255",
    "E3_T0"="#DDDDDD","E3_CTR"="#88CCEE","E3_ALG"="#332288")) +
  scale_fill_manual(values=c(
    "E1_T0"="#AAAAAA","E1_CTR"="#999933","E1_ALG"="#117733",
    "E2_T0"="#000000","E2_CTR"="#CC6677","E2_ALG"="#882255",
    "E3_T0"="#DDDDDD","E3_CTR"="#88CCEE","E3_ALG"="#332288")) +
  theme_classic() + 
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(PCoA_ALL_DNA_v_RNA_exp_cond_colored)
dev.copy(pdf,'FigS1A - PCoA_ALL_DNA_v_RNA_exp_cond_colored.pdf', width = 4, height = 3.6)
dev.off()

# Figure S1B - DNA and RNA colored by type (i.e., DNA or RNA)
PCoA_ALL_DNA_v_RNA_type_colored = amp_ordinate(
  amp_subset_samples(ampvis),
  type = "PCoA",
  distmeasure = "bray",
  transform = "none", 
  sample_color_by = "type", sample_colorframe = T,
  sample_shape_by = "fraction_TRUE",
  sample_label_size = 3,
  sample_point_size = 3) +
  theme_classic() + 
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.ticks.x = element_line(colour="black", size=0.1),
        axis.ticks.y = element_line(colour="black", size=0.1),
        strip.text.x = element_text(size=12, colour="black"),
        strip.text.y = element_text(size=12, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "#000000"),
        legend.position = "none")
print(PCoA_ALL_DNA_v_RNA_type_colored)
dev.copy(pdf,'FigS1B - PCoA_ALL_DNA_v_RNA_type_colored.pdf', width = 4, height = 3.6)
dev.off()

##### Figure S2 ##################################
# RNA-based community composition

top20_families = amp_heatmap(amp_subset_samples(
  ampvis, type=="RNA"),
  tax_aggregate = "Family",
  tax_add = "Class",
  showRemainingTaxa = TRUE,
  tax_show = 20,
  normalise = T,
  plot_values = F,
  min_abundance = 0,
  max_abundance = 15,
  plot_colorscale = "sqrt",
  facet_by = c("site","fraction_TRUE","condition"),
  group_by = "day",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, colour = "black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=9, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(top20_families)
dev.copy(pdf,'FigS2 - top20_families_heatmap_01Nov2022_RNA_version_wnumbers.pdf', width = 8.5, height = 3.7)
dev.off()

##### Figure S3 ##################################
# E1
top20_families_E1 = amp_heatmap(amp_subset_samples(
  ampvis, exp == "E1" & type=="DNA" ),
  tax_aggregate = "Family",
  tax_add = "Class",
  showRemainingTaxa = TRUE,
  tax_show = 20,
  normalise = T,
  plot_values = T,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 15,
  plot_colorscale = "sqrt",
  facet_by = c("fraction_TRUE","condition", "day"),
  group_by = "clip_id",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=9, colour="black", angle=45, hjust = 1),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=9, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(top20_families_E1)
dev.copy(pdf,'FigS3A - top20families_E1_replicates_wnumbers_2023.11.11.pdf', width = 10, height = 6)
dev.off()

# E2
top20_families_E2 = amp_heatmap(amp_subset_samples(
  ampvis, exp == "E2" & type=="DNA" ),
  tax_aggregate = "Family",
  tax_add = "Class",
  showRemainingTaxa = TRUE,
  tax_show = 20,
  normalise = T,
  plot_values = T,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 15,
  plot_colorscale = "sqrt",
  facet_by = c("fraction_TRUE","condition", "day"),
  group_by = "clip_id",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=9, colour="black", angle=45, hjust = 1),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=9, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(top20_families_E2)
dev.copy(pdf,'FigS3B - top20families_E2_replicates_wnumbers_2023.11.11.pdf', width = 10, height = 6)
dev.off()

# E3
top20_families_E3 = amp_heatmap(amp_subset_samples(
  ampvis, exp == "E3" & type=="DNA" ),
  tax_aggregate = "Family",
  tax_add = "Class",
  showRemainingTaxa = TRUE,
  tax_show = 20,
  normalise = T,
  plot_values = T,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 15,
  plot_colorscale = "sqrt",
  facet_by = c("fraction_TRUE","condition", "day"),
  group_by = "clip_id",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=9, colour="black", angle=45, hjust = 1),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=9, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(top20_families_E3)
dev.copy(pdf,'FigS3C - top20families_E3_replicates_wnumbers_2023.11.11.pdf', width = 10.3, height = 6)
dev.off()



##### Figure S4 ###################################

# Differentially abundant PA ASVs, expanding and showing relative abundance
PA_rar_upreg_rel_abund_heatmap = amp_heatmap(amp_subset_samples(
  only_PA_rar_upreg_combined_ASVs, type=="DNA"),
  tax_aggregate = "OTU",
  tax_add = "Genus",
  showRemainingTaxa = TRUE,
  tax_show = 25,
  normalise = F,
  plot_values = T,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 10,
  plot_colorscale = "sqrt",
  facet_by = c("site","fraction_TRUE","condition"),
  group_by = "day",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, colour="black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=10, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(PA_rar_upreg_rel_abund_heatmap)
dev.copy(pdf,'FigS4A - PA_rar_upreg_rel_abund_heatmap.pdf', width = 11, height = 8)
dev.off()

# Differentially abundant FL ASVs, expanding and showing relative abundance
FL_rar_upreg_rel_abund_heatmap = amp_heatmap(amp_subset_samples(
  only_FL_rar_upreg_combined_ASVs, type=="DNA"),
  tax_aggregate = "OTU",
  tax_add = "Genus",
  showRemainingTaxa = TRUE,
  tax_show = 25,
  normalise = F,
  plot_values = T,
  plot_values_size = 2.5,
  min_abundance = 0,
  max_abundance = 10,
  plot_colorscale = "sqrt",
  facet_by = c("site","fraction_TRUE","condition"),
  group_by = "day",
  color_vector = c(
    "white", "gray","#3b528b","#21918c","#5ec962", "#FFFF00"))+
  theme_classic() +
  theme(axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, colour="black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=10, colour="black"),
        strip.text.y = element_text(size=7, colour="black"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "#ffffff"),
        legend.position = "right") 
print(FL_rar_upreg_rel_abund_heatmap)
dev.copy(pdf,'FigS4B - FL_rar_upreg_rel_abund_heatmap.pdf', width = 11, height = 8)
dev.off()


##### ADDITIONAL ANALYSES #########################
# count # of bacterial classes and families in PA and FL upreg ASV dataset

unique(rar_PA_upreg$Family) # 63
unique(rar_PA_upreg$Class) #18

FL_upreg_ASVs.df = as.data.frame(tax_table(onlyFL_rar_upreg_combined_physeq))
unique(rar_FL_upreg$Family) # 26
unique(rar_FL_upreg$Class) #9

# Counting Gamma, Alpha, and Bacteroidia ASVs
nrow(subset(rar_PA_upreg, Class == "Gammaproteobacteria")) #110
nrow(subset(rar_FL_upreg, Class == "Gammaproteobacteria")) #51

nrow(subset(rar_PA_upreg, Class == "Alphaproteobacteria")) #81
nrow(subset(rar_FL_upreg, Class == "Alphaproteobacteria")) #13

nrow(subset(rar_PA_upreg, Class == "Bacteroidia")) #69
nrow(subset(rar_FL_upreg, Class == "Bacteroidia")) #5

