###################################################################
## GWAS data (06_data_gwas_ukbb_pqtls)                           ##
##                                                               ##
## James Staley                                                  ##
## Email: jrstaley95@ucb.com                                     ##
##                                                               ##
## 14/11/23                                                      ##
###################################################################

###################################################################
##### Set-up #####
###################################################################

##### Clear #####
rm(list = ls())

##### Options #####
options(stringsAsFactors = FALSE)

##### Libraries #####
library(dplyr, warn.conflicts = FALSE)

##### Functions #####
source("./scripts/00_functions.R")

###################################################################
##### Data #####
###################################################################

##### LD SNPs #####
load("./data/01_data_ldmat.Rda")

##### GWAS #####
load("./data/02_data_gwas_sost_region.Rda")
rm(pqtls)

##### pQTLs #####
pqtls <- data.table::fread(
  "./data/05_data_ukbb_pqtls.tsv",
  header = TRUE, data.table = FALSE, sep = "\t"
)
pqtls <- pqtls %>%
  as_tibble() %>%
  mutate(
    chr = as.character(chr)
  )

###################################################################
##### Process pQTLs #####
###################################################################

##### pQTLs #####
pqtls <- pqtls %>%
  select(-c(rsid)) %>%
  inner_join(
    x = ld_snps,
    y = .,
    by = c("chr", "pos")
  ) %>%
  filter((ref == oa & alt == ea) | (ref == ea & alt == oa)) %>%
  mutate(
    af = if_else(ref == ea & alt == oa, 1 - eaf, eaf),
    beta = if_else(ref == ea & alt == oa, -1 * beta, beta)
  ) %>%
  select(rsid, chr, pos, ref, alt, af, beta, se, pvalue)

###################################################################
##### Process GWAS #####
###################################################################

##### GWAS #####
gwas <- gwas %>%
  inner_join(
    x = .,
    y = select(pqtls, rsid, chr, pos, ref, alt),
    by = c("rsid", "chr", "pos", "ref", "alt")
  )

##### Save #####
save(pqtls, gwas, studies, file = "./data/06_data_gwas_ukbb_pqtls.Rda")
