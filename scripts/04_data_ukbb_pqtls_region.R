###################################################################
## GWAS UKBB pQTLs (04_data_ukbb_pqtls_region)                   ##
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

###################################################################
##### Data #####
###################################################################

##### SOST (build 37) #####
# * An extra 25KB upstream of the gene is added to ensure that all variants
#   that might be needed downstream (i.e., variants associated with CD300LG
#   expression) are extracted
chrom <- "17"
start <- 41831099 - 100000
end <- 41836156 + 125000

##### UKBB pQTL variants #####
snps <- fread(
  paste0(
    Sys.getenv("downloads_folder"),
    "olink_rsid_map_mac5_info03_b0_7_chr17_patched_v2.tsv.gz"
  ),
  header = TRUE,
  data.table = FALSE,
  sep = "\t"
)
snps <- snps %>%
  as_tibble() %>%
  filter(
    POS19 >= !!start &
      POS19 <= !!end
  )

##### UKBB pQTLs #####
pqtls <- fread(
  paste0(
    Sys.getenv("downloads_folder"),
    "SOST_Q9BQB4_OID20204_v1_Cardiometabolic/",
    "discovery_chr17_SOSTQ9BQB4OID20204v1Cardiometabolic.gz"
  ),
  header = TRUE,
  data.table = FALSE,
  sep = " "
)
pqtls <- pqtls %>%
  as_tibble() %>%
  inner_join(
    x = .,
    y = snps,
    by = "ID"
  ) %>%
  rename(
    chr = CHROM,
    pos = POS19,
    ea = ALLELE1,
    oa = ALLELE0,
    eaf = A1FREQ,
    beta = BETA,
    se = SE,
    n = N
  ) %>%
  mutate(
    rsid = gsub("_", ":", rsid),
    eaf = round(eaf, 4),
    beta = round(beta, 6),
    se = round(se, 6),
    pvalue = signif(10^(-LOG10P), 4)
  ) %>%
  select(rsid, chr, pos, ea, oa, eaf, beta, se, pvalue, n) %>%
  arrange(pos)

##### Save #####
fwrite(
  pqtls, "./data/04_data_ukbb_pqtls_region.tsv",
  quote = FALSE, eol = "\n", sep = "\t"
)
