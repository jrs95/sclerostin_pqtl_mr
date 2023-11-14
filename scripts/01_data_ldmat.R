###################################################################
## 1000G genotypes LD matrix (01_data_ldmat)                     ##
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

##### SOST (build 37) #####
# * An extra 25KB upstream of the gene is added to ensure that all variants
#   that might be needed downstream (i.e., variants associated with CD300LG
#   expression) are extracted
chrom <- "17"
start <- 41831099 - 100000
end <- 41836156 + 125000

##### Genotypes #####

## Extract
vcf <- extract_vcf(
  file = paste0(
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    "ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  ),
  chr = chrom,
  start = start,
  end = end,
  remove_index = TRUE
)
vcf <- vcf %>%
  filter(FILTER == "PASS")

## Filter
vcf <- filter_vcf(
  data = vcf,
  samples_file = paste0(
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    "integrated_call_samples_v3.20130502.ALL.panel"
  )
)

## rsID
rsids <- extract_vcf(
  file = paste0(
    "https://ftp.ensembl.org/pub/grch37/release-110/variation/vcf/",
    "homo_sapiens/1000GENOMES-phase_3.vcf.gz"
  ),
  chr = chrom,
  start = start,
  end = end,
  index = "csi",
  remove_index = TRUE
)

## VCF
vcf <- vcf %>%
  left_join(
    select(rsids, `#CHROM`, POS, ID, REF, ALT),
    by = c("#CHROM", "POS", "REF"),
    suffix = c("", "_RSID")
  ) %>%
  mutate(
    ALT_EQUAL = sapply(
      seq_len(nrow(.)),
      function(i) grepl(ALT[i], ALT_RSID[i])
    ),
    ID = if_else(ALT_EQUAL, ID_RSID, ID, ID)
  ) %>%
  select(-c(ID_RSID, ALT_RSID, ALT_EQUAL)) %>%
  mutate(
    ID = if_else(
      (duplicated(ID) | duplicated(ID, fromLast = TRUE)) & ID != ".",
      paste0(ID, ":", REF, ":", ALT),
      ID
    ),
    ID = if_else(
      ID == ".",
      paste0(`#CHROM`, ":", POS, ":", REF, ":", ALT),
      ID
    )
  )

###################################################################
##### Process #####
###################################################################

##### LD matrix #####

## SNPs
ld_snps <- vcf %>%
  select(ID, `#CHROM`, POS, REF, ALT, INFO) %>%
  rename(
    rsid = ID,
    chr = `#CHROM`,
    pos = POS,
    ref = REF,
    alt = ALT,
    af = INFO
  ) %>%
  mutate(
    pos = as.integer(pos),
    af = as.numeric(sub("AF=", "", af))
  )

## Correlation matrix
genos <- vcf %>%
  select(10:ncol(.)) %>%
  as.matrix()
rownames(genos) <- ld_snps$rsid
genos[genos == "0|0"] <- 0
genos[genos %in% c("0|1", "1|0")] <- 1
genos[genos == "1|1"] <- 2
class(genos) <- "numeric"
ld_mat <- round(cor(t(genos)), 4)

##### Save #####
save(ld_snps, ld_mat, file = "./data/01_data_ldmat.Rda")
