###################################################################
## GWAS UKBB pQTLs (05_data_ukbb_pqtls)                          ##
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
library(data.table)
library(dplyr, warn.conflicts = FALSE)
library(LDlinkR)

###################################################################
##### Data #####
###################################################################

##### SOST (build 37) #####
chrom <- "17"
start <- 41831099 - 100000
end <- 41836156 + 100000

##### UKBB pQTLs #####
pqtls <- fread(
  "./data/04_data_ukbb_pqtls_region.tsv",
  header = TRUE,
  data.table = FALSE,
  sep = "\t"
)
pqtls <- pqtls %>%
  as_tibble() %>%
  arrange(pvalue) %>%
  filter(
    chr == !!chrom &
      pos >= !!start &
      pos <= !!end
  ) %>%
  mutate(
    snp = if_else(
      grepl("^rs", rsid),
      rsid,
      paste0("chr", chr, ":", pos)
    )
  ) %>%
  filter(pvalue <= 5e-8)

##### LD pruning #####
pqtls_ld <- SNPclip(
  pqtls$snp,
  token = Sys.getenv("ldlink_token"),
  r2_threshold = "0.2"
)
pqtls <- pqtls_ld %>%
  as_tibble() %>%
  filter(grepl("kept", Details)) %>%
  rename(
    snp = RS_Number
  ) %>%
  select(snp) %>%
  inner_join(
    x = pqtls,
    y = .,
    by = "snp"
  ) %>%
  select(-snp)

##### Save #####
fwrite(
  pqtls, "./data/05_data_ukbb_pqtls.tsv",
  quote = FALSE, eol = "\n", sep = "\t"
)
