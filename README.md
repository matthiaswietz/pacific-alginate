# Bacterial community dynamics during alginate degradation in the South, Equatorial, and North Pacific Ocean

This repo describes the bioinformatic workflow for the manuscript "Distinct bacterial succession and functional response to alginate in the South, Equatorial, and North Pacific Ocean" by Balmonte et al. (doi XXX). Here, we investigated alginate-induced compositional and functional shifts of particle-associated and free-living bacterial communities along a Pacific Ocean transect Samples originate from mesocosm experiments onboard expedition SO248 [BacGeoPac]() with RV Sonne. DNA- and RNA-based reads were processed into amplicon sequence variants using DADA2. 

The repo contains:

- [DADA2 workflow](./SO248_Alg_DADA2.Rmd)
- Resulting [ASV](./Alg-DOM_v138.1_seqtab.txt) and [taxonomy](./Alg-DOM_v138.1_tax.txt) tables
- [ASV sequences](./Alg-DOM_v138.1_asv.fasta)
- [Metadata](./metadata.txt)
- [Script for analyses and plotting](./metadata.R)

