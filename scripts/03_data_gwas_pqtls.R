###################################################################
## GWAS data (03_data_gwas_pqtls)                                ##
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

###################################################################
##### Process GWAS #####
###################################################################

##### Proxy variants #####
# * GCST006867 (T2D) does not contain all of the sclerostin pQTLs (rs66838809 &
#   rs76449013 are not in the dataset)
# * There are no good proxies in the SOST region of these variants (r2 > 0.8)
#   in the GCST006867 dataset
# * The IEU OpenGWAS server also finds no good proxies for GCST006867 (r2 > 0.8)
ld_proxy <- round(ld_mat[, c("rs66838809", "rs76449013"), drop = FALSE]^2, 4)
ld_proxy <- ld_proxy %>%
  as.data.frame() %>%
  mutate(
    proxy = rownames(.)
  ) %>%
  relocate(proxy, .before = 1) %>%
  as_tibble() %>%
  tidyr::pivot_longer(
    cols = c(rs66838809, rs76449013),
    names_to = "rsid",
    values_to = "r2"
  ) %>%
  relocate(rsid, .before = 1) %>%
  filter(r2 >= 0.8) %>%
  arrange(desc(r2)) %>%
  inner_join(
    x = select(ld_snps, rsid),
    y = .,
    by = "rsid"
  ) %>%
  rename(snp = rsid)

gwas %>%
  filter(id == "GCST006867") %>%
  filter(rsid %in% !!ld_proxy$proxy) %>%
  print()

##### Conditional analyses #####
# * CD300LG nearby to SOST is known to be associated with HDL cholesterol &
#   triglycerides
# * Since rs1107747 & rs4793023 are mildly correlated with rs72836567
#   (the top associated variant with CD300LG expression in GTEx, r2 > 0.05),
#   we condition rs1107747 & rs4793023 by rs72836567 for HDL cholesterol &
#   triglycerides
# * This is an approximate analysis conditioning one variant by
#   another for effect sizes only based on the COJO methodology
# * This function gives very similar (often exactly the same) results
#   as the COJO implementation in GCTA

## pQTLs + CD300LG variant
gwas <- gwas %>%
  filter(rsid %in% !!c(unique(pqtls$rsid), "rs72836567"))

## COJO
for (id in c("UKB30760", "UKB30870")) {
  for (snp in c("rs1107747", "rs4793023")) {
    gwas$beta[gwas$id == id & gwas$rsid == snp] <- cojo(
      gwas$beta[gwas$id == id & gwas$rsid == "rs72836567"],
      gwas$beta[gwas$id == id & gwas$rsid == snp],
      gwas$af[gwas$id == id & gwas$rsid == "rs72836567"],
      gwas$af[gwas$id == id & gwas$rsid == snp],
      ld_mat[rownames(ld_mat) == "rs72836567", rownames(ld_mat) == snp],
      studies$n[studies$id == id]
    )
  }
}

## Process
gwas <- gwas %>%
  mutate(
    beta = round(beta, 6),
    pvalue = if_else(
      id %in% c("UKB30760", "UKB30870") &
        rsid %in% c("rs1107747", "rs4793023"),
      signif(2 * pnorm(-abs(beta / se)), 4),
      pvalue
    )
  )

##### GWAS #####
gwas <- gwas %>%
  inner_join(
    x = .,
    y = select(pqtls, rsid, chr, pos, ref, alt),
    by = c("rsid", "chr", "pos", "ref", "alt")
  )

##### Save #####
save(pqtls, gwas, studies, file = "./data/03_data_gwas_pqtls.Rda")
