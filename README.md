# Bacterial community dynamics during alginate degradation in the South, Equatorial, and North Pacific Ocean

This repo describes the bioinformatic workflow for the manuscript "Distinct bacterial succession and functional response to alginate in the South, Equatorial, and North Pacific Ocean" by Balmonte and colleagues (doi XXX). Here, we report alginate-induced dynamics of particle-associated and free-living bacterial communities along a Pacific Ocean transect, based on mesocosm experiments onboard expedition SO248 (BacGeoPac) with RV Sonne. 

DNA- and RNA-based reads were processed into amplicon sequence variants, derived from mesocosms with addition of (i) alginate or (ii) algal DOM. Since mesocosms were set up with the same ambient seawater, including shared controls, DADA2 was run on all files together (separated by MiSeq run). This will subsequently allow a coherent interpretation, since the same ASV IDs can be tracked across incubation regimes. For the present paper, the complet ASV table was then subsetted to samples from alginate mesocosms.

The repo contains:

- [ASV generation](./dada.Rmd) and [chimera removal / taxonomy assignment](./MergeChimTax.R)
- Resulting [ASV](./Alg-DOM_v138.1_seqtab.txt) and [taxonomy](./Alg-DOM_v138.1_tax.txt) tables
- [ASV sequences](./Alg-DOM_v138.1_asv.fasta)
- [Metadata](./metadata.txt)
- [Script for analyses and plotting](./SO248_Alginate.R)

