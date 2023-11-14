###################################################################
## GWAS data (02_data_gwas_sost_region)                          ##
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
options(timeout = 1200)

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

##### pQTLs #####
pqtls <- data.table::fread(
  "./data/00_data_pqtls.tsv",
  header = TRUE, data.table = FALSE, sep = "\t"
)
pqtls <- pqtls %>%
  as_tibble() %>%
  mutate(
    chr = as.character(chr)
  )

##### Studies #####
studies <- data.table::fread(
  "./data/00_data_studies.tsv",
  header = TRUE, data.table = FALSE, sep = "\t"
)
studies <- studies %>%
  as_tibble()

##### LD SNPs #####
load("./data/01_data_ldmat.Rda")
ld_snps <- ld_snps %>%
  select(-c(af))

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

##### GCST #####
studies_gcst <- studies %>%
  filter(source == "GCST")
gcst <- tibble()
for (i in seq_len(nrow(studies_gcst))) {
  gcst <- gcst %>%
    bind_rows(
      gwas_sumstats(
        file = studies_gcst$link[i],
        id = studies_gcst$id[i],
        chrom = chrom, start = start, end = end,
        chromosome = "chromosome",
        position = "base_pair_location",
        effect_allele = "effect_allele",
        other_allele = "other_allele",
        effect_allele_frequency = "effect_allele_frequency",
        effect = "beta",
        standard_error = "standard_error",
        p_value = "p_value"
      )
    )
  print(paste0("Dataset: ", studies_gcst$id[i], " -- processed"))
}

##### UKB #####

## Neale
# * These datasets were downloaded (the extensions were changed from .bgz
#   to .gz) from the Neale UK Biobank server using the links in the studies
#   manifest.
# * Replace Sys.getenv("downloads_folder") with the file path to the
#   folder where the *.gwas.imputed_v3.both_sexes.varorder.tsv.gz files are
#   located
studies_ukb <- studies %>%
  filter(source == "UKB")
ukb <- tibble()
for (i in seq_len(nrow(studies_ukb))) {
  data_ukb <- data.table::fread(
    paste0(
      Sys.getenv("downloads_folder"),
      sub(".bgz$", ".gz", sub(".*\\/", "", studies_ukb$link[i]))
    ),
    header = TRUE, data.table = FALSE, sep = "\t", showProgress = FALSE
  )
  data_ukb <- data_ukb %>%
    as_tibble() %>%
    filter(
      grepl("17:", variant) &
        low_confidence_variant == FALSE
    ) %>%
    tidyr::separate(
      variant,
      into = c("chr", "pos", "ref", "alt"),
      sep = ":"
    ) %>%
    mutate(
      eaf = if_else(
        minor_allele == alt,
        minor_AF,
        1 - minor_AF
      )
    )
  ukb <- ukb %>%
    bind_rows(
      gwas_sumstats(
        data = data_ukb,
        id = studies_ukb$id[i],
        chrom = chrom, start = start, end = end,
        chromosome = "chr",
        position = "pos",
        effect_allele = "alt",
        other_allele = "ref",
        effect_allele_frequency = "eaf",
        effect = "beta",
        standard_error = "se",
        p_value = "pval"
      )
    )
  print(paste0("Dataset: ", studies_ukb$id[i], " -- processed"))
}

## Pan-UKBB
# * This dataset was downloaded (the extension was changed from .bgz to .gz)
#   from the Pan-UKBB server using the link in the studies manifest
# * Replace Sys.getenv("downloads_folder") with the file path to the
#   folder where the categorical-20002-both_sexes-1065.tsv.gz file is located
studies_pukb <- studies %>%
  filter(source == "PUKB")
ukb <- ukb %>%
  bind_rows(
    gwas_sumstats(
      file = paste0(
        Sys.getenv("downloads_folder"),
        sub(".bgz$", ".gz", sub(".*\\/", "", studies_pukb$link))
      ),
      id = studies_pukb$id,
      chrom = chrom, start = start, end = end,
      chromosome = "chr",
      position = "pos",
      effect_allele = "alt",
      other_allele = "ref",
      effect_allele_frequency = "af_controls_meta_hq",
      effect = "beta_meta_hq",
      standard_error = "se_meta_hq",
      p_value = "neglog10_pval_meta_hq",
      p_value_transform = TRUE
    )
  )
print(paste0("Dataset: ", studies_pukb$id, " -- processed"))

##### MEGASTROKE #####
# * Since chromosome positions are not available in this dataset we add
#   them using the 1000G SNP panel
studies_megastroke <- studies %>%
  filter(source == "MEGASTROKE")
megastroke <- tibble()
for (i in seq_len(nrow(studies_megastroke))) {
  data_megastroke <- data.table::fread(
    studies_megastroke$link[i],
    header = TRUE, data.table = FALSE, sep = "\t", showProgress = FALSE
  )
  data_megastroke <- data_megastroke %>%
    as_tibble() %>%
    inner_join(
      x = .,
      y = ld_snps,
      by = c("variant_id" = "rsid")
    ) %>%
    filter(
      (ref == toupper(other_allele) & alt == toupper(effect_allele)) |
        (ref == toupper(effect_allele) & alt == toupper(other_allele))
    ) %>%
    select(-c(ref, alt))
  megastroke <- megastroke %>%
    bind_rows(
      gwas_sumstats(
        data = data_megastroke,
        id = studies_megastroke$id[i],
        chrom = chrom, start = start, end = end,
        rsid = "variant_id",
        chromosome = "chr",
        position = "pos",
        effect_allele = "effect_allele",
        other_allele = "other_allele",
        effect_allele_frequency = "effect_allele_frequency",
        effect = "beta",
        standard_error = "standard_error",
        p_value = "p_value"
      )
    )
  print(paste0("Dataset: ", studies_megastroke$id[i], " -- processed"))
}

##### DIAGRAM #####
# * This dataset was downloaded from the DIAGRAM downloads page
#   using the links in the studies manifest
# * Replace Sys.getenv("downloads_folder") with the file path to the
#   folder where the DIAMANTE-TA.sumstat.txt.gz file is located
studies_diagram <- studies %>%
  filter(source == "DIAGRAM")
diagram <- gwas_sumstats(
  file = paste0(Sys.getenv("downloads_folder"), "DIAMANTE-TA.sumstat.txt.gz"),
  id = studies_diagram$id,
  chrom = chrom, start = start, end = end,
  rsid = "rsID",
  chromosome = "chromosome(b37)",
  position = "position(b37)",
  effect_allele = "effect_allele",
  other_allele = "other_allele",
  effect = "Fixed-effects_beta",
  standard_error = "Fixed-effects_SE",
  p_value = "Fixed-effects_p-value",
  sep = " "
)
print(paste0("Dataset: ", studies_diagram$id, " -- processed"))

##### CHARGE #####
# * This dataset was downloaded from the CHARGE dbGaP repository
#   (Accession: phs000930) under application #35610
# * Replace Sys.getenv("downloads_folder") with the file path to the
#   folder where the METAL_AAC_Fullresults.txt.gz file is located
# * The effect sizes and standard errors are estimated using the
#   Z-scores and allele frequencies
#   (with SE = 1 / sqrt(2 * af * (1 - af) * (n + z^2)))
file_charge <- paste0(
  Sys.getenv("downloads_folder"), "METAL_AAC_Fullresults.txt.gz"
)
if (file.exists(file_charge)) {
  studies_charge <- studies %>%
    filter(source == "CHARGE")
  data_charge <- data.table::fread(
    file = file_charge,
    header = TRUE, data.table = FALSE, sep = "\t", showProgress = FALSE
  )
  data_charge <- data_charge %>%
    as_tibble() %>%
    mutate(
      SE = 1 / sqrt(2 * AVG_freq1 * (1 - AVG_freq1) * (N + Zscore^2)),
      BETA = Zscore * SE
    )
  charge <- gwas_sumstats(
    data = data_charge,
    id = studies_charge$id,
    chrom = chrom, start = start, end = end,
    rsid = "Markername",
    chromosome = "Chr",
    position = "Position",
    effect_allele = "allele_1",
    other_allele = "allele_2",
    effect_allele_frequency = "AVG_freq1",
    effect = "BETA",
    standard_error = "SE",
    p_value = "P_META"
  )
  print(paste0("Dataset: ", studies_charge$id, " -- processed"))
} else {
  charge <- tibble()
}

##### GEFOS #####
studies_gefos <- studies %>%
  filter(source == "GEFOS")
gefos <- tibble()
for (i in seq_len(nrow(studies_gefos))) {
  gefos <- gefos %>%
    bind_rows(
      gwas_sumstats(
        file = studies_gefos$link[i],
        id = studies_gefos$id[i],
        chrom = chrom, start = start, end = end,
        rsid = c("RSID", "SNP")[i],
        chromosome = "CHR",
        position = "BP",
        effect_allele = c("EA", "ALLELE1")[i],
        other_allele = c("NEA", "ALLELE0")[i],
        effect_allele_frequency = c("EAF", "A1FREQ")[i],
        effect = c("BETA", "logOR")[i],
        standard_error = c("SE", "logOR.SE")[i],
        p_value = "P"
      )
    )
  print(paste0("Dataset: ", studies_gefos$id[i], " -- processed"))
}

##### Overall #####
gwas <- bind_rows(gcst, ukb, megastroke, diagram, charge, gefos)

###################################################################
##### Harmonize #####
###################################################################

##### Min p-value #####
# * All sclerostin pQTLs are available on the 1000G SNP panel
# * However, not all variants in the SOST region are on the 1000G SNP panel,
#   so we compute the minimum p-value in SOST region prior to harmonization to
#   ensure that we are not dropping any important associations
# * The minimum p-value in the SOST region is only different pre and
#   post harmonization for:
#   * HTN (min(pvalue)_preh = 0.00208; min(pvalue)_posth = 0.0592) and
#   * ApoB (min(pvalue)_preh = 0.000692; min(pvalue)_posth = 0.00419)
#   * SVS (min(pvalue)_preh = 0.0172; min(pvalue)_posth = 0.0197)
# * As these associations are not associated at regional-wide level
#   (p-value > 1E-6) we can safely restrict to variants on the 1000G SNP panel

## Min p-value pre-harmonization
minp_preh <- gwas %>%
  mutate(
    variant = paste0(chr, ":", pos, ":", oa, ":", ea)
  )  %>%
  filter(
    pos <= 41836156 + 100000
  ) %>%
  group_by(id) %>%
  summarize(
    variant = variant[which.min(pvalue)],
    pvalue = min(pvalue, na.rm = TRUE),
    .groups = "drop"
  )

## Min p-value post-harmonization
minp_posth <- gwas %>%
  mutate(
    variant = paste0(chr, ":", pos, ":", oa, ":", ea)
  ) %>%
  select(-c(rsid)) %>%
  inner_join(
    x = ld_snps,
    y = .,
    by = c("chr", "pos")
  ) %>%
  filter(
    pos <= 41836156 + 100000 &
      ((ref == oa & alt == ea) | (ref == ea & alt == oa))
  ) %>%
  group_by(id) %>%
  summarize(
    variant = variant[which.min(pvalue)],
    pvalue = min(pvalue, na.rm = TRUE),
    .groups = "drop"
  )

## Min p-value common
minp <- minp_preh %>%
  inner_join(
    x = .,
    y = minp_posth,
    by = "id",
    suffix = c("_preh", "_posth")
  ) %>%
  filter(variant_preh != variant_posth)

print(minp)

##### Harmonize #####
# * This aligns the other and effect alleles to the reference and
#   alternative allele in the 1000G SNP panel, repectively
gwas <- gwas %>%
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
  select(id, rsid, chr, pos, ref, alt, af, beta, se, pvalue) %>%
  inner_join(
    x = select(studies, id),
    y = .,
    by = "id"
  )

##### Studies #####
studies <- studies %>%
  select(-c(link))

##### Save #####
save(pqtls, gwas, studies, file = "./data/02_data_gwas_sost_region.Rda")
