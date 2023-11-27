## Merge all runs ##
## Chimera removal, taxonomy with Silva138

require(dada2)
require(ShortRead)
require(ggplot2)
require(gridExtra)

setwd("/isibhv/projects/FRAMdata/mwietz/SO248_meso/Alg-DOM/")

# Load individual seqtabs
sq1 <- readRDS(
  "seqtab_MQ180801.rds")
sq2 <- readRDS(
  "seqtab_MQ180819.rds")
sq3 <- readRDS(
  "seqtab_MQ181221.rds")
sq4 <- readRDS(
  "seqtab_MQ201204.rds")
sq5 <- readRDS(
  "seqtab_MQ210126.rds")

# merge all years/seq-runs
seqtab_full <- mergeSequenceTables(
  sq1,sq2,sq3,sq4,sq5, repeats = "error")

# removing chimeras 
# Identified 155001 bimeras out of 190190 input sequences.
seqtab.nochim <- removeBimeraDenovo(
  seqtab_full, method = "consensus", 
  multithread = T, verbose = T)

# lots of bimeras - but still OK?
sum(seqtab.nochim)/sum(seqtab_full) #0.82 - OK!

# stats
dim(seqtab.nochim)  # 35189 sequences in 302 samples
summary(rowSums(seqtab.nochim)/rowSums(seqtab_full))

# Determine read lengths/size range of amplicons
table(rep(nchar(colnames(seqtab.nochim)), 
          colSums(seqtab.nochim)))

# Remove singletons and junk sequences
# "c" adjusted to size range of amplicons
# ~100k reads are only 265bp
# seems strange - keep in mind
seqtab.nochim2 <- seqtab.nochim[, nchar(
  colnames(seqtab.nochim)) %in% c(382:449) & 
    colSums(seqtab.nochim) > 1]

# stats
dim(seqtab.nochim2) # 14755 sequences
summary(rowSums(seqtab.nochim2))
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))

################################################################

## TAXONOMY ##

# Silva v138
tax1 <- assignTaxonomy(
  seqtab.nochim2, 
  "../../../MolObs/tax_db/silva_nr_v138_train_set.fa.gz", 
  tryRC = TRUE,
  multithread = TRUE)

# outcome: Bacteria 14675 -- Archaea 71 -- Eukaryota 7
table(tax1[, 1])            

# Remove NA on phylum level
# Remove Chloroplasts; mitochondria
sum(is.na(tax1[, 2]))   #132
tax.good <- tax1[!is.na(tax1[, 2]) & tax1[, 1] %in% c(
  "Bacteria","Archaea"),]
tax.good <- tax.good[-c(grep(
  "Chloroplast", tax.good[, 4]), grep(
  "Mitochondria", tax.good[, 5])), ]
seqtab.nochim2.good <- seqtab.nochim2[, rownames(tax.good)]
summary(rowSums(seqtab.nochim2.good))

# Format tables
seqtab.nochim2.print <- t(seqtab.nochim2.good)
tax.print <- tax.good
all.equal(
  rownames(seqtab.nochim2.print), 
  rownames(tax.print)) #TRUE
rownames(seqtab.nochim2.print) <- paste(
  "asv", 1:ncol(seqtab.nochim2.good), sep = "")
rownames(tax.print) <- rownames(seqtab.nochim2.print)

# Export
write.table(
  seqtab.nochim2.print, "Alg-DOM_v138_seqtab.txt", sep="\t", quote=F)
write.table(
  tax.print, "Alg-DOM_v138_tax.txt", sep="\t", quote=F)
uniquesToFasta(
  seqtab.nochim2.good, "Alg-DOM_v138_asv.fasta")

#######################################

# Silva v138.1 with species
tax2 <- assignTaxonomy(
  seqtab.nochim2, 
  "../../../MolObs/tax_db/silva_nr99_v138.1_wSpecies_train_set.fa.gz", 
  tryRC = TRUE,
  multithread = TRUE)

# outcome: Bacteria 14684 -- Archaea 71 
table(tax2[, 1])   

# Remove NA on phylum level
# Remove Chloroplasts; mitochondria
sum(is.na(tax2[, 2]))   #95
tax.good <- tax2[!is.na(tax2[, 2]) & tax2[, 1] %in% c(
  "Bacteria","Archaea"),]
tax.good <- tax.good[-c(grep(
  "Chloroplast", tax.good[, 4]), grep(
    "Mitochondria", tax.good[, 5])), ]
seqtab.nochim2.good <- seqtab.nochim2[, rownames(tax.good)]
summary(rowSums(seqtab.nochim2.good))

# Format tables
seqtab.nochim2.print <- t(seqtab.nochim2.good)
tax.print <- tax.good
all.equal(
  rownames(seqtab.nochim2.print), 
  rownames(tax.print)) #TRUE
rownames(seqtab.nochim2.print) <- paste(
  "asv", 1:ncol(seqtab.nochim2.good), sep = "")
rownames(tax.print) <- rownames(seqtab.nochim2.print)

# Export
write.table(
  seqtab.nochim2.print, "Alg-DOM_v138.1_seqtab.txt", sep="\t", quote=F)
write.table(
  tax.print, "Alg-DOM_v138.1_tax.txt", sep="\t", quote=F)
uniquesToFasta(
  seqtab.nochim2.good, "Alg-DOM_v138.1_asv.fasta")


# save taxdata
save.image("SO248_Alg-DOM_tax.Rdata")

###################################################################################
###################################################################################
