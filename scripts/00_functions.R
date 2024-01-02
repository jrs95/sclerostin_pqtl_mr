#' @title Extract VCF
#'
#' @description Extract regions from tabix-indexed VCF files.
#'
#' @name extract_vcf
#'
#' @import dplyr
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importFrom IRanges IRanges
#'
#' @importFrom Rsamtools scanTabix TabixFile
#'
#' @importFrom stringr str_count str_split_fixed
#'
#' @param file the file to query, this can either be a filepath or a URL
#'
#' @param chr genomic region chromosome
#'
#' @param start genomic region start position (build 37)
#'
#' @param end genomic region end position (build 37)
#'
#' @param index tabix index file type (i.e. tbi or csi); default: `"tbi"`
#'
#' @param remove_index remove index file after the query is performed,
#'   this is useful when extracting data using the URL; default: `FALSE`
#'
#' @param header remove index file after the query is performed,
#'   this is useful when extracting data using the URL; default: `FALSE`
#'
#' @return A `data.frame` in VCF format.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
extract_vcf <- function(file, chr, start, end, index = "tbi",
  remove_index = FALSE, header = "vcf") {

  # Region
  region <- GenomicRanges::GRanges(
    seqnames = as.character(chr),
    ranges = IRanges::IRanges(start = start, end = end),
    strand = "*"
  )

  # Query
  results <- Rsamtools::scanTabix(
    Rsamtools::TabixFile(file, index = paste0(file, ".", index)),
    param = region
  )[[1]]
  ncols <- stringr::str_count(results[1], pattern = "\t") + 1
  results <- results %>%
    stringr::str_split_fixed(pattern = "\t", n = ncols) %>%
    as.data.frame() %>%
    tibble::as_tibble()

  # Header
  if (!is.null(header) && header == "vcf") {
    header <- Rsamtools::headerTabix(
      Rsamtools::TabixFile(file, index = paste0(file, ".", index))
    )
    header <- tail(header$header, n = 1)
    header <- header %>%
      stringr::str_split_fixed(pattern = "\t", n = ncols)
    results <- results %>%
      setNames(header)
  } else if (!is.null(header) && header == "tsv") {
    header <- readr::read_tsv(file, n_max = 0, show_col_types = FALSE)
    results <- results %>%
      setNames(names(header))
  }

  # Remove index
  if (remove_index)
    file.remove(paste0(sub(".*/", "", file), ".", index))

  # Return
  return(results)

}

#' @title Filter VCF
#'
#' @description Filter 1000G VCF files for samples of a particular ancestry and
#'   variants by minor allele frequency and Hardy-Weinburg equilibrium.
#'
#' @name filter_vcf
#'
#' @import dplyr
#'
#' @importFrom data.table fread
#'
#' @importFrom genetics genotype HWE.test
#'
#' @param data `data.frame` with variants in VCF format
#'
#' @param samples_file the filepath or URL for the 1000G samples file
#'
#' @param ancestry the ancestry to keep; default: `"EUR"`
#'
#' @param maf minor allele frequency threshold; default: `0.01`
#'
#' @param hwe hardy-weinburg equilibrium; default: `1e-10`
#'
#' @return A filtered `data.frame` in VCF format.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
filter_vcf <- function(data, samples_file, ancestry = "EUR", maf = 0.01,
  hwe = 1e-10) {

  # Samples
  samples <- data.table::fread(
    input = samples_file, header = TRUE, data.table = FALSE,
    fill = TRUE, sep = "\t", showProgress = FALSE
  )
  samples <- samples %>%
    dplyr::filter(super_pop != !!ancestry) %>%
    dplyr::pull(sample)

  # Results
  results <- data %>%
    dplyr::select(-any_of(samples)) %>%
    dplyr::filter(!grepl(",", REF) & !grepl(",", ALT)) %>%
    dplyr::mutate(
      INFO = sub(paste0(".*", ancestry, "_AF="), "AF=", INFO),
      INFO = sub(";.*", "", INFO),
      MAF = pmin(
        as.numeric(sub("AF=", "", INFO)),
        1 - as.numeric(sub("AF=", "", INFO))
      )
    ) %>%
    dplyr::filter(MAF >= !!maf) %>%
    dplyr::select(-MAF)

  # HWE
  genos <- results %>%
    dplyr::select(-c(1:9)) %>%
    as.data.frame()
  hardy <- function(x) {
    genetics::HWE.test(genetics::genotype(sub("\\|", "/", x)))$test$p.value
  }
  hwe_pvalue <- apply(genos, 1, hardy)
  results <- results %>%
    dplyr::mutate(HWE = !!hwe_pvalue) %>%
    dplyr::filter(HWE >= !!hwe) %>%
    dplyr::select(-HWE)

  # Return
  return(results)

}

#' @title GWAS summary statistics format
#'
#' @description Extract and format GWAS summary statistics for a set of
#'   rsIDs, chromosome positions or a genomic region
#'
#' @name gwas_sumstats
#'
#' @import dplyr
#'
#' @importFROM data.table fread
#'
#' @param file the file containing the GWAS summary statistics, this can either
#'   be a filepath or a URL
#'
#' @param file the file containing the GWAS summary statistics, this can either
#'   be a filepath or a URL
#'
#' @param id study identifier
#'
#' @param snps rsIDs; default: `NULL`
#'
#' @param chrpos chromosome positions (build 37) in the format chr + ':' + pos;
#'   default: `NULL`
#'
#' @param chrom genomic region chromosome; default: `NULL`
#'
#' @param start genomic region start position (build 37); default: `NULL`
#'
#' @param end genomic region end position (build 37); default: `NULL`
#'
#' @param rsid name of the column in the file that depicts the rsID;
#'   default: `""`
#'
#' @param chromosome name of the column in the file that depicts the chromosome
#'
#' @param position name of the column in the file that depicts the position
#'   (build 37)
#'
#' @param effect_allele name of the column in the file that depicts the effect
#'   allele
#'
#' @param other_allele name of the column in the file that depicts the
#'   other allele
#'
#' @param effect_allele_frequency name of the column in the file that depicts
#'   the effect allele frequency; default: `""`
#'
#' @param effect name of the column in the file that depicts the effect size
#'
#' @param standard_error name of the column in the file that depicts the
#'   standard error of the effect size
#'
#' @param p_value name of the column in the file that depicts the p-value of
#'   the association
#'
#' @param p_value_transform transform the p-value from -log10(p-value) to
#'   p-value; default: `FALSE`
#'
#' @param sep file separator; default: `"\t"`
#'
#' @return A `data.frame` with the following columns:
#'
#' \item{id}{study identifier}
#'
#' \item{rsid}{rsID}
#'
#' \item{chr}{chromosome}
#'
#' \item{pos}{position}
#'
#' \item{ea}{effect allele}
#'
#' \item{oa}{other allele}
#'
#' \item{eaf}{effect allele frequency}
#'
#' \item{beta}{association effect size}
#'
#' \item{se}{standard error of beta}
#'
#' \item{pvalue}{association p-value}
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
gwas_sumstats <- function(data = NULL, file = NULL, id,
  snps = NULL, chrpos = NULL, chrom = NULL, start = NULL, end = NULL,
  rsid = "", chromosome, position, effect_allele, other_allele,
  effect_allele_frequency = "", effect, standard_error, p_value,
  p_value_transform = FALSE, sep = "\t") {

  # Data
  if (!is.null(data)) {
    data <- data %>%
      tibble::as_tibble()
  } else if (is.null(data) && !is.null(file)) {
    data <- data.table::fread(
      file, header = TRUE, data.table = FALSE, sep = sep, showProgress = FALSE
    )
    data <- data %>%
      tibble::as_tibble()
  } else {
    stop("both data and file cannot both be null")
  }

  # Process
  if (!is.null(snps)) {
    data <- data %>%
      dplyr::rename(
        rsid = !!rsid,
        chr = !!chromosome,
        pos = !!position
      ) %>%
      dplyr::mutate(
        rsid = as.character(rsid),
        chr = as.character(chr),
        pos = as.integer(pos)
      ) %>%
      dplyr::filter(rsid %in% !!snps)
  } else if (!is.null(chrpos)) {
    data <- data %>%
      dplyr::rename(
        chr = !!chromosome,
        pos = !!position
      ) %>%
      dplyr::mutate(
        chr = as.character(chr),
        pos = as.integer(pos),
        chr_pos = paste0(data$chr, ":", data$pos)
      ) %>%
      dplyr::filter(chr_pos %in% !!chrpos) %>%
      dplyr::select(-c(chr_pos))
    if (rsid != "") {
      data <- data %>%
        dplyr::rename(
          rsid = !!rsid
        ) %>%
        dplyr::mutate(
          rsid = as.character(rsid)
        )
    } else {
      data <- data %>%
        dplyr::mutate(
          rsid = "."
        )
    }
  } else if (!is.null(chrom) && !is.null(start) && !is.null(end)) {
    data <- data %>%
      dplyr::rename(
        chr = !!chromosome,
        pos = !!position
      ) %>%
      dplyr::mutate(
        chr = as.character(chr),
        pos = as.integer(pos)
      ) %>%
      dplyr::filter(
        chr == !!as.character(chrom) &
          pos >= !!start &
          pos <= !!end
      )
    if (rsid != "") {
      data <- data %>%
        dplyr::rename(
          rsid = !!rsid
        ) %>%
        dplyr::mutate(
          rsid = as.character(rsid)
        )
    } else {
      data <- data %>%
        dplyr::mutate(
          rsid = "."
        )
    }
  } else {
    stop("chrpos and chr, start and end should not both be null")
  }
  data <- data %>%
    dplyr::rename(
      ea = !!effect_allele,
      oa = !!other_allele,
      beta = !!effect,
      se = !!standard_error,
      pvalue = !!p_value
    ) %>%
    dplyr::mutate(
      id = !!id,
      rsid = if_else(grepl("^rs", rsid), rsid, ".", "."),
      ea = toupper(ea),
      oa = toupper(oa),
      beta = as.numeric(beta),
      beta = round(beta, 6),
      se = as.numeric(se),
      se = round(se, 6)
    )
  if (p_value_transform) {
    data <- data %>%
      dplyr::mutate(
        pvalue = as.numeric(pvalue),
        pvalue = signif(10^(-pvalue), 4)
      ) %>%
      filter(!is.na(pvalue))
  } else {
    data <- data %>%
      dplyr::mutate(
        pvalue = as.numeric(pvalue),
        pvalue = signif(pvalue, 4)
      ) %>%
      filter(!is.na(pvalue))
  }
  if (effect_allele_frequency != "") {
    data <- data %>%
      dplyr::rename(
        eaf = !!effect_allele_frequency
      ) %>%
      dplyr::mutate(
        eaf = as.numeric(eaf)
      )
  } else {
    data <- data %>%
      dplyr::mutate(eaf = NA)
  }
  data <- data %>%
    dplyr::select(
      id, rsid, chr, pos, ea, oa, eaf, beta, se, pvalue
    )

  # Return
  return(data)

}

#' @title GWAS impute summary statistics
#'
#' @description Impute GWAS summary statistics for a quantitative trait
#'   using the top genetic variant association and LD information
#'
#' @name gwas_impute
#'
#' @import dplyr
#'
#' @param data a `data.frame` of GWAS summary results for the top genetic
#'   variant association with the following columns:
#'   `rsid` (rsID),
#'   `chr` (chromosome),
#'   `pos` (position),
#'   `ea` (effect allele),
#'   `oa` (other allele),
#'   `eaf` (effect allele frequency),
#'   `beta` (association effect size),
#'   `se` (standard error of beta),
#'   `pvalue` (association p-value)
#'
#' @param snpinfo variant information for genetic variants in the region to be
#'   imputed with the following columns:
#'   `rsid` (rsID),
#'   `chr` (chromosome),
#'   `pos` (position),
#'   `ref` (reference allele),
#'   `alt` (alternate allele),
#'   `af` (alternate allele frequency)
#'
#' @param corr correlation matrix for the variants in `snpinfo`
#'
#' @param n number of samples
#'
#' @param n_cases number of cases; default: `NULL`
#'
#' @return A `data.frame` with the following columns:
#'
#' \item{rsid}{rsID}
#'
#' \item{chr}{chromosome}
#'
#' \item{pos}{position}
#'
#' \item{ref}{reference allele}
#'
#' \item{alt}{alternate allele (effect allele)}
#'
#' \item{af}{alternate allele frequency}
#'
#' \item{beta}{association effect size}
#'
#' \item{se}{standard error of beta}
#'
#' \item{pvalue}{association p-value}
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
gwas_impute <- function(data, snpinfo, corr, n, n_cases = NULL) {

  # Data

  ## GWAS
  gwas <- data %>%
    dplyr::select(-rsid) %>%
    dplyr::arrange(as.numeric(pvalue)) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(z = beta / se) %>%
    dplyr::left_join(
      x = snpinfo,
      y = .,
      by = c("chr", "pos")
    ) %>%
    dplyr::mutate(
      af_study = if_else(ref == ea & alt == oa, 1 - eaf, eaf),
      af_study = if_else(
        (ref == oa & alt == ea) | (ref == ea & alt == oa),
        af_study,
        NA
      ),
      z = if_else(ref == ea & alt == oa, -1 * z, z),
      z = if_else(
        (ref == oa & alt == ea) | (ref == ea & alt == oa),
        z,
        NA
      )
    ) %>%
    dplyr::select(rsid, chr, pos, ref, alt, af, af_study, z, se)

  if (any(duplicated(gwas$rsid)))
    stop("rsIDs have to be unique")
  if (all(is.na(gwas$z)))
    stop("all z-statistics are missing")

  ## Correlation
  if (any(rownames(corr) != gwas$rsid))
    stop(
      "the variants in the correlation matrix and the GWAS dataset are not ",
      "the same"
    )
  r <- corr[rownames(corr) == gwas$rsid[!is.na(gwas$z)], , drop = TRUE]
  if (!is.vector(r) || length(r) != nrow(gwas))
    stop("r is not a vector or does not have the same length as gwas has rows")

  # Impute statistics

  ## Z-statistics
  zi <- r * gwas$z[!is.na(gwas$z)]
  gwas <- gwas %>%
    dplyr::mutate(
      zi = zi
    )

  ## Standard errors
  if (is.null(n_cases)) {
    gwas <- gwas %>%
      dplyr::mutate(
        se = if_else(
          is.na(se),
          sqrt(1 / (2 * af * (1 - af) * !!n)),
          se
        )
      )
  } else {
    gwas <- gwas %>%
      dplyr::mutate(
        se = if_else(
          is.na(se),
          sqrt(1 / (2 * af * (1 - af) * (!!n - !!n_cases) * !!n_cases / !!n)),
          se
        )
      )
  }

  ## Effect sizes and p-values
  gwas <- gwas %>%
    dplyr::mutate(
      beta = zi * se,
      pvalue = 2 * pnorm(-abs(beta / se))
    )

  ## Process results
  gwas <- gwas %>%
    dplyr::mutate(
      af_alt = if_else(is.na(af_study), af, af_study),
      beta = round(beta, 6),
      se = round(se, 6),
      pvalue = signif(pvalue, 4)
    )

  # Results
  gwas <- gwas %>%
    dplyr::select(rsid, chr, pos, ref, alt, af_alt, beta, se, pvalue) %>%
    dplyr::rename(af = af_alt)

  # Return
  return(gwas)

}

#' @title Conditional single SNP analysis
#'
#' @description Perform an analysis for SNP 2 conditioned on SNP 1
#'   (effect size only). Assumes that the effect alleles for the two SNPs
#'   are the same across all data inputted.
#'
#' @details This function is based on the COJO methodology developed by
#'   J Yang for GCTA
#'   (PMID: [22426310](https://pubmed.ncbi.nlm.nih.gov/22426310/))
#'
#' @name cojo
#'
#' @param beta1 effect size for SNP 1
#'
#' @param beta2 effect size of SNP 2
#'
#' @param eaf1 effect allele frequency for SNP 1
#'
#' @param eaf2 effect allele frequency for SNP 2
#'
#' @param corr correlation between SNP 1 and SNP 2
#'
#' @param n sample size
#'
#' @return The effect size for SNP 2 adjusted for SNP 1.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
cojo <- function(beta1, beta2, eaf1, eaf2, corr, n) {

  # D (diagonal X'X)
  d1 <- 2 * eaf1 * (1 - eaf1) * n
  d2 <- 2 * eaf2 * (1 - eaf2) * n

  # C
  cc <- sqrt(d1 * d2) * corr

  # Beta2 adjusted for Beta1
  beta2_beta1 <- beta2 - cc * beta1 / d2

  # Return
  return(beta2_beta1)

}

#' @title Regional plot
#'
#' @description Regional plot of associations in a genomic region.
#'
#' @name regional_plot
#'
#' @import dplyr
#'
#' @import geni
#'
#' @import ggiraph
#'
#' @param data `data.frame` with rsID (`rsid`), chromosome (`chr`),
#'   position (`pos`) and p-values (`pvalue`)
#'
#' @param corr correlation `matrix` between variants
#'
#' @param highlights variants to highlight; default: `NULL`
#'
#' @param highlights_title highlighted variants legend title; default: `NULL`
#'
#' @param highlights_colour highlighted variants point colour;
#'   default: `#416F6F`
#'
#' @param highlights_label label highlighted variants points;
#'   default: `TRUE`
#'
#' @param top_marker the variant to plot the correlation statistics of the rest
#'   of the variants against; default: `NULL`
#'
#' @param thresh p-value threshold(s); default: `c(5e-8, 1e-6)`
#'
#' @param thresh_colour p-value threshold(s) colour;
#'   default: `c("darkgreen", "darkred")`
#'
#' @return `regional_plot` returns an interactive regional plot.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
regional_plot <- function(data, corr, highlights = NULL,
  highlights_title = NULL, highlights_colour = "#416F6F",
  highlights_label = TRUE, top_marker = NULL,
  thresh = c(5e-8, 1e-6), thresh_colour = c("darkgreen", "darkred")) {

  # Correlation matrix
  corr <- corr[
    match(data$rsid, rownames(corr)),
    match(data$rsid, colnames(corr)),
    drop = FALSE
  ]

  # Regional plot
  figure <- data %>%
    dplyr::rename(marker = rsid, pvalue = pvalue) %>%
    dplyr::mutate(pos = as.integer(pos)) %>%
    dplyr::select(marker, chr, pos, pvalue) %>%
    fig_region(
      data = .,
      corr = corr,
      build = 37,
      interactive = TRUE,
      top_marker = top_marker,
      highlights = highlights,
      highlights_cat = ifelse(
        !is.null(highlights_title),
        paste0(highlights_title, "   "),
        highlights_title
      ),
      highlights_colours = highlights_colour,
      highlights_label = highlights_label,
      alpha = 0.4,
      highlights_title = "",
      thresh = thresh,
      thresh_colour = thresh_colour,
      axis_text_size = 12,
      axis_title_size = 13,
      genebar_plot_size = 16,
      girafe = FALSE
    )

  # Interactive plot
  figure <- girafe(
    print(figure),
    width_svg = 10, height_svg = 7,
    options = list(
      opts_tooltip(
        css = paste0(
          "background-color:#EEEEEE;",
          "color:black;",
          "padding:10px;",
          "border-radius:5px;"
        )
      ),
      opts_toolbar(position = "topleft"),
      opts_sizing(width = 0.95)
    )
  )

  # Return
  return(figure)

}

#' @title QQ plot
#'
#' @description Quantile-quantile plot of observed vs expected p-values.
#'
#' @name regional_plot
#'
#' @import dplyr
#'
#' @import geni
#'
#' @import ggiraph
#'
#' @param data `data.frame` with rsID (`rsid`) and p-values (`pvalue`)
#'
#' @return `qq_plot` returns an interactive QQ plot.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
qq_plot <- function(data) {

  # Data
  if ("trait" %in% names(data) && "rsid" %in% names(data)) {
    data <- data %>%
      dplyr::mutate(
        text = paste0(
          "Trait: ", trait,
          "<br>rsID: ", rsid,
          "<br>p-value: ", signif(pvalue, 3)
        )
      )
  }  else if ("trait" %in% names(data) && !("rsid" %in% names(data))) {
    data <- data %>%
      dplyr::mutate(
        text = paste0(
          "Trait: ", trait,
          "<br>p-value: ", signif(pvalue, 3)
        )
      )
  } else if (!("trait" %in% names(data)) && "rsid" %in% names(data)) {
    data <- data %>%
      dplyr::mutate(
        text = paste0(
          "rsID: ", rsid,
          "<br>p-value: ", signif(pvalue, 3)
        )
      )
  } else {
    data <- data %>%
      dplyr::mutate(
        text = paste0(
          "p-value: ", signif(pvalue, 3)
        )
      )
  }

  # QQ plot
  figure <- data %>%
    dplyr::select(pvalue, text) %>%
    fig_qq(
      data = .,
      ci = TRUE,
      interactive = TRUE,
      axis_text_size = 10,
      axis_title_size = 11,
      girafe = FALSE
    )

  # Interactive plot
  figure <- girafe(
    print(figure),
    width_svg = 5, height_svg = 5,
    options = list(
      opts_tooltip(
        css = paste0(
          "background-color:#EEEEEE;",
          "color:black;",
          "padding:10px;",
          "border-radius:5px;"
        )
      ),
      opts_toolbar(position = "topleft"),
      opts_sizing(width = 0.6)
    )
  )

  # Return
  return(figure)

}


#' @title Colocalization analysis
#'
#' @description `coloc` performs a colocalization analysis using
#'   the log-ABFs and priors provided.
#'
#' @name coloc
#'
#' @details This function was extracted and adapted from the
#'   \href{https://github.com/chr1swallace/coloc}{coloc} package,
#'   written by Chris Wallace & Claudia Giambartolomei.
#'
#' @param bf1 vector of (log) bayes factors for trait 1
#'
#' @param bf2 vector of (log) bayes factors for trait 2
#'
#' @param p1 prior probability a variant is associated with trait 1;
#'   default: `1e-4`
#'
#' @param p2 prior probability a variant is associated with trait 2;
#'   default: `1e-4`
#'
#' @param p12 prior probability a variant is associated with both traits;
#'   default: `2e-6`
#'
#' @return `coloc` returns a `vector` of colocalization
#'   posterior probabilities.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
coloc <- function(bf1, bf2, p1 = 1e-4, p2 = 1e-4, p12 = 2e-6) {

  # Bayes factors for hypotheses

  ## HO
  h0_bf <- 0

  ## H1 and H2
  h1_bf <- log(p1) + coloc_logsum(bf1)
  h2_bf <- log(p2) + coloc_logsum(bf2)

  ## H3
  h3_bf <- log(p1) + log(p2) +
    coloc_logdiff(
      coloc_logsum(bf1) + coloc_logsum(bf2),
      coloc_logsum(bf1 + bf2)
    )

  ## H4
  h4_bf <- log(p12) + coloc_logsum(bf1 + bf2)

  # Results
  bfs <- c(h0_bf, h1_bf, h2_bf, h3_bf, h4_bf)
  denom <- coloc_logsum(bfs)
  res <- exp(bfs - denom)
  names(res) <- paste0("pph", (seq_along(res)) - 1)

  # Return
  return(res)

}

#' @title Colocalization approximate Bayes factors
#'
#' @description `coloc_bf` computes log approximate Bayes factors.
#'
#' @name coloc_bf
#'
#' @details This function was extracted and adapted from the
#'   \href{https://github.com/chr1swallace/coloc}{coloc} package,
#'   written by Chris Wallace & Claudia Giambartolomei
#'
#' @param z z-statitic of coefficients
#'
#' @param v variance of coefficients
#'
#' @param binary binary phenotype; default: `FALSE`
#'
#' @param pheno_var phenotypic variance; default: `1`
#'
#' @return `coloc_bf` returns a `vector` of Bayes factors.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
coloc_bf <- function(z, v, binary = FALSE, pheno_var = 1) {

  # Variance prior
  if (binary == TRUE) {
    sd_prior <- 0.2
  } else {
    sd_prior <- 0.15 * sqrt(pheno_var)
  }

  # Bayes factors (log)
  r <- sd_prior^2 / (sd_prior^2 + v)
  bf <- 0.5 * (log(1 - r) + (r * z^2))

  # Return
  return(bf)

}

#' @title Colocalization analysis log sum
#'
#' @description `coloc_logsum` computes the log of the sum of the
#'   exponentiated logs taking out the max, i.e. insuring that the sum is
#'   not infinite.
#'
#' @name coloc_logsum
#'
#' @details This function was extracted and adapted from the
#'   \href{https://github.com/chr1swallace/coloc}{coloc} package,
#'   written by Chris Wallace & Claudia Giambartolomei
#'
#' @param x numeric `vector`
#'
#' @return `coloc_logsum` returns the result of:
#'   `max(x) + log(sum(exp(x - max(x))))`.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @noRd
#' @md
coloc_logsum <- function(x) {

  # Compute log-sum
  res <- max(x) + log(sum(exp(x - max(x))))

  # Return
  return(res)

}

#' @title Colocalization analysis log difference
#'
#' @description `coloc_logdiff` computes the log of the
#'   difference of the exponentiated logs taking out the max,
#'   i.e. insuring that the difference is not negative.
#'
#' @name coloc_logdiff
#'
#' @details This function was extracted and adapted from the
#'   \href{https://github.com/chr1swallace/coloc}{coloc} package,
#'   written by Chris Wallace & Claudia Giambartolomei
#'
#' @param x numeric `vector`
#'
#' @param y numeric `vector`
#'
#' @return `coloc_logdiff` returns the results of:
#'   `max(x, y) + log(exp(x - max(x, y)) - exp(y - max(x, y)))`.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @noRd
#' @md
coloc_logdiff <- function(x, y) {

  # Compute log-diff
  res <- max(x, y) + log(exp(x - max(x, y)) - exp(y - max(x, y)))

  # Return
  return(res)

}

#' @title Empty table
#'
#' @description Empty table for captions in document.
#'
#' @name empty_kable
#'
#' @import dplyr
#'
#' @importFrom knitr kable kable_styling
#'
#' @return `empty_kable` returns an empty `kable` table.
#'
#' @author James Staley <jrstaley95@ucb.com>
#'
#' @export
#' @md
empty_kable <- function() {

  # Results
  results <- knitr::kable(tibble(" " = "")[c(), ], col.names = NULL) %>%
    kableExtra::kable_styling(full_width = TRUE)

  # Return
  return(results)

}
