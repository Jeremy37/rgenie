#' rgenie: Analysis of GenIE experiments
#'
#' rgenie is a package to analyze the next-generation sequencing output from a
#' set of GenIE experimental replicates.
#'
#' GenIE (genome-editing interrogation of enhancers) is an experimental method
#' to evaluate the effects of individual SNPs on gene transcription. It is based
#' on targeted CRISPR-Cas9 genome editing in cultured cells to produce indels at
#' a target locus, optionally with a homology-dependent recombination (HDR)
#' construct to create precisely defined genomic changes. Following this, both
#' genomic DNA (gDNA) and RNA are extracted from the same pool of cells, and
#' then cDNA is generated. Both gDNA and cDNA are amplified with locus-specific
#' PCR primers, in multiple replicates, and libraries of the experimental
#' replicates are prepared for next-generation sequencing. The sequencing data
#' from each replicate for each region should be aligned to either the full
#' reference genome or to a locus-specific (amplicon) reference genome for
#' analysis by rgenie.
#'
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @import stringr
#'
#' @section rgenie analyses: Two types of analyses can be run with rgenie: a
#'   grep-based analysis, which uses pattern matching to identify HDR vs.
#'   wild-type (WT) reads, and an alignment-based analysis, which uses all
#'   Cas9-generated alleles, including deletions. The analysis is run for a
#'   given \strong{region}, which must have multiple \strong{replicates} for
#'   both cDNA and gDNA. These analyses take as input a data.frame defining the
#'   region, and a data.frame defining replicates.
#'
#'   The result of an rgenie analysis is a list object with multiple tables,
#'   which depend on whether a grep or deletion analyses was done. Plotting
#'   functions take as input the result of one of these main rgenie analyses.
#'
#' @docType package
#' @name rgenie
#' @seealso \code{\link{grep_analysis}}
#' @seealso \code{\link{deletion_analysis}}
NULL


# Checks that the options are valid, and fills in any needed additional opts
check_del_opts = function(opts) {
  if (is.null(opts$crispr_del_window) || opts$crispr_del_window < 1) {
    stop("Option 'crispr_del_window' must have a positive value.")
  }
  if (opts$del_span_start >= opts$del_span_end) {
    stop("Option 'del_span_start' should be an integer that is less than option 'del_span_end'.")
  }
  return( check_main_opts(opts) )
}

check_main_opts = function(opts) {
  if (is.null(opts$required_match_left) || is.null(opts$required_match_right) ||
      opts$required_match_left < 1 || opts$required_match_right < 1) {
    stop("Options 'required_match_left' abd 'required_match_right' must have positive values.")
  }
  return(opts)
}

# Checks that the region and replicate IDs are distinct
check_regions_and_replicates = function(regions, replicates) {
  if (any(duplicated(regions$name))) {
    duplicated_region = regions$name[ which(duplicated(regions$name))[1] ]
    stop(sprintf("All region names must be unique. Found duplicated region name: %s.", duplicated_region))
  }

  # Check that each replicate has a distinct name; otherwise it will
  # be impossible to distinguish replicates in the output files
  for (region_name in regions$name) {
    df = replicates %>% dplyr::filter(name == region_name)
    if (any(duplicated(df$replicate))) {
      stop(sprintf("Region %s: not all replicates have distinct names. It will be impossible to distinguish results in output files that relate to specific replicates.", region_name))
    }
  }
  if (!all(grepl("cDNA|gDNA", replicates$type))) {
    stop("All replicates must have column 'type' equal to either 'cDNA' or 'gDNA'.")
  }
}


#' Grep-based GenIE analysis
#'
#' For each replicate associated with an input region, \code{grep_analysis}
#' searches for alleles matching the HDR or WT sequences, and returns statistics
#' that indicate whether the HDR:WT ratio differs in cDNA and gDNA.
#'
#' @param regions A data.frame defining GenIE regions.
#' @param replicates A data.frame defining GenIE replicates.
#' @param required_match_left The length of sequence to the left of the HDR site
#'   that must exactly match to identify HDR or WT reads.
#' @param required_match_right The length of sequence to the right of the HDR
#'   site that must exactly match to identify HDR or WT reads.
#' @param min_mapq The minimum mapping quality for reads to be included in the
#'   analysis.
#' @section Details: For a grep analysis, the regions parameter is a data.frame
#'   with a format as follows. All of the column names below must be specified.
#'   \tabular{lccccclll}{ name \tab sequence_name \tab start \tab end \tab
#'   highlight_site \tab cut_site \tab hdr_allele_profile \tab wt_allele_profile
#'   \tab ref_sequence \cr MUL1_rs6700034  \tab MUL1  \tab 1 \tab 21 \tab 11
#'   \tab 9 \tab ----------A---------- \tab ----------C---------- \tab
#'   ACCGCACCCCCCCGGCCTAAC \cr }
#'
#' \tabular{rl}{
#'   \strong{name} \tab A unique identifier for the region. \cr
#'   \strong{sequence_name} \tab The chromosome or amplicon sequence name. \cr
#'   \strong{start} \tab The start coordinate of the amplicon relative to the chromosome or amplicon reference. \cr
#'   \strong{end} \tab The end coordinate of the amplicon relative to the chromosome or amplicon reference (the end coordinate is included in the region). \cr
#'   \strong{highlight_site} \tab The relative position of the site of interest, usually the HDR SNP site. \cr
#'   \strong{cut_site} \tab The relative position of the cut site. \cr
#'   \strong{ref_sequence} \tab The reference sequence for the amplicon region, which must have length (end - start + 1). \cr
#'   \strong{hdr_allele_profile} \tab An allele profile describing the HDR allele. See details below. \cr
#'   \strong{wt_allele_profile} \tab An allele profile describing the WT allele. See details below. \cr
#' }
#'
#' If multiple rows are defined for `regions`, then a separate analysis is run
#' for each region, using matched replicates from `replicates`.
#'
#' The allele_profile columns indicate the positions in the amplicon sequence
#' that must match a given nucleotide for a read to be considered either HDR
#' or WT. This sequence must be the same length as the reference sequence,
#' and all other positions should be "-". The total sequence region that must
#' match is determined by both the positions of specified nucleotides and by
#' the \strong{required_match_left} and \strong{required_match_right}
#' parameters. These parameters give the length of sequence which must match
#' the provided reference sequence to the left of the leftmost specified
#' nucleotide, or to the right of the rightmost specified nucleotide.
#'
#'
#' The replicates parameter is a data.frame with a format as below.
#' \tabular{lccl}{ name \tab replicate \tab type \tab bam \cr MUL1_rs6700034
#' \tab c1.2 \tab cDNA \tab
#' bam_amplicon/MUL1_rs6700034_cDNA_rep1_pcr2.sortedByCoord.bam \cr
#' MUL1_rs6700034 \tab c1.3 \tab cDNA \tab
#' bam_amplicon/MUL1_rs6700034_cDNA_rep1_pcr3.sortedByCoord.bam \cr }
#'
#' \tabular{rl}{
#'   \strong{name} \tab Indicates the region that a given replicate
#'   corresponds with. All replicates matching the name in the regions table
#'   will be used. \cr
#'   \strong{replicate} \tab an ID for the replicate, which must
#'   be unique among replicates for the region. \cr
#'   \strong{type} \tab Must have the value "cDNA" or "gDNA", indicating
#'   whether a given replicate contains data for cDNA or gDNA. \cr
#'   \strong{bam} \tab the path (relative to the working directory) to
#'   a BAM file with sequencing reads for the replicate. \cr
#' }
#'
#' Statistics can only be computed if there are at least 2 replicates of each
#' type (cDNA and gDNA). Replicates are matched to the region based on the
#' \strong{name} column.
#'
#' @section Results: The returned object is a list, where each item is the
#' result for one region. The result for a region (e.g. results[[1]]) is itself
#' a list, with the following items:
#' \tabular{rl}{
#'   \strong{region_stats} \tab Main analysis output, with statistics indicating whether the HDR/WT levels differ in cDNA relative to gDNA. \cr
#'   \strong{replicate_stats} \tab A data.frame with a row for each replicate, which has counts of reads in different categories and some summary values. \cr
#'   \strong{region} \tab 	Details of the input region the result corresponds to. \cr
#'   \strong{replicates} \tab Details of the input replicates the result corresponds to. \cr
#'   \strong{opts} \tab A list containing the options that were given for the analysis. \cr
#'   \strong{type} \tab Has the value “grep_analysis”, and is used by plotting functions that take a full grep_result list as input. \cr
#' }
#'
#' The main output of interest is the `region_stats` field, which is a one-row data.frame with the following values:
#'
#' \tabular{rl}{
#'   \strong{name} \tab Name of the region. \cr
#'   \strong{effect} \tab Estimated effect size - the amount by which the HDR:WT ratio differs in cDNA relative to gDNA. \cr
#'   \strong{effect_sd} \tab Standard deviation of the estimated effect size. \cr
#'   \strong{effect_confint_lo} \tab Lower bound of the 95\% confidence interval for the effect size. \cr
#'   \strong{effect_confint_hi} \tab Upper bound of the 95\% confidence interval for the effect size. \cr
#'   \strong{df_estimated} \tab The degrees of freedom used in the unequal variances t-test, which is estimated from the data. \cr
#'   \strong{pval} \tab A p value from the unequal variants t-test. \cr
#' }
#'
#' @examples
#' grep_results = grep_analysis(regions, replicates)
#' grep_summary_plot(grep_results[[1]])
#' @seealso \code{\link{deletion_analysis}}
#' @seealso \code{\link{grep_summary_plot}}
#' @export
#'
grep_analysis = function(regions,
                         replicates,
                         required_match_left = 10,
                         required_match_right = 10,
                         min_mapq = 0,
                         quiet = FALSE)
{
  check_regions_and_replicates(regions, replicates)
  replicates$replicate = as.character(replicates$replicate)
  region_names = unique(regions$name)

  opts = list(regions,
              replicates,
              required_match_left = required_match_left,
              required_match_right = required_match_right,
              min_mapq = min_mapq,
              quiet = quiet)
  opts = check_main_opts(opts)

  grep_results = list()
  region_index = 1
  for (region_name in region_names) {
    cur_region = regions %>% dplyr::filter(name == region_name)
    cur_replicates = replicates %>% dplyr::filter(name == region_name)
    if (nrow(cur_region) != 1) {
      stop(sprintf("grepAnalysis: Expected a single region row for region %s.", cur_region$name))
    }
    if (nrow(cur_replicates) == 0) {
      warning(sprintf("No replicates found for region %s", region_name))
      next
    }

    grep_result = region_grep_analysis(cur_region, cur_replicates, opts)
    grep_result$region = cur_region
    grep_result$replicates = cur_replicates

    grep_results[[region_index]] = grep_result
    region_index = region_index + 1
  }

  return( grep_results )
}


# Runs a GenIE grep analysis for a single region and its replicates
region_grep_analysis = function(region, replicates, opts)
{
  opts$region = region
  opts$replicates = replicates

  hdr_profile = str_to_upper(region$hdr_allele_profile)
  wt_profile = str_to_upper(region$wt_allele_profile)
  ref_sequence = str_to_upper(region$ref_sequence)
  sequence_name = region$sequence_name

  hdr_profile_chars = strsplit(hdr_profile, "")[[1]]
  wt_profile_chars = strsplit(wt_profile, "")[[1]]
  ref_seq_chars = strsplit(ref_sequence, "")[[1]]

  get_grep_sequence = function(profile_chars, ref_chars) {
    profile_is_letter = is_dna_letter(profile_chars)
    char_positions = which(profile_is_letter)
    minpos = max(1, min(char_positions) - opts$required_match_left)
    maxpos = min(length(ref_chars), max(char_positions) + opts$required_match_right)
    grep_chars = ref_chars
    grep_chars[profile_is_letter] = profile_chars[profile_is_letter]
    paste(grep_chars[minpos:maxpos], collapse = "")
  }
  hdr_seq_grep = get_grep_sequence(hdr_profile_chars, ref_seq_chars)
  wt_seq_grep = get_grep_sequence(wt_profile_chars, ref_seq_chars)
  output(opts, sprintf("HDR grep sequence: %s\n", hdr_seq_grep))
  output(opts, sprintf("WT grep sequence:  %s\n", wt_seq_grep))

  replicate_grep_analyses = list()
  for (i in 1:nrow(replicates)) {
    bam_file = replicates$bam[i]
    replicate = replicates$replicate[i]
    type = replicates$type[i]

    output(opts, sprintf("\n\nAnalysing region %s, replicate %s, file %s\n", region$name, replicate, bam_file, region$sequence_name))
    output(opts, sprintf("HDR allele: %s\n", hdr_profile))
    output(opts, sprintf("WT allele:  %s\n", wt_profile))
    output(opts, sprintf("REF sequence: %s\n", ref_sequence))

    # Get different reference for cDNA if specified
    if (hasName(region, "ref_sequence_cdna") && !is_empty(region$ref_sequence_cdna)) {
      # If a different ref sequence is specified, all the following columns must also be specified.
      if (is_empty(region$sequence_name_cdna) || is_empty(region$hdr_allele_profile_cdna) || is_empty(region$wt_allele_profile_cdna)) {
        stop(sprintf("ref_sequence_cdna specified in regions file for locus %s, but not all of the following columns were specified: sequence_name_cdna, hdr_allele_profile_cdna, wt_allele_profile_cdna", region$name))
      }
      if (type == "cDNA") {
        hdr_profile = str_to_upper(region$hdr_allele_profile_cdna)
        wt_profile = str_to_upper(region$wt_allele_profile_cdna)
        ref_sequence = str_to_upper(region$ref_sequence_cdna)
        sequence_name = region$sequence_name_cdna
      }
    }

    read_data = get_region_reads_from_bam(region$name, replicate, bam_file, sequence_name, region$start, region$end, ref_sequence, opts$min_mapq)
    if (length(read_data$seq) <= 0) {
      stop(sprintf("\nERROR: No reads retrieved at the specified region for replicate %s from file %s\n", replicate, bam_file))
    }

    result = replicate_grep_analysis(read_seqs = read_data$seq,
                                     hdr_seq_grep = hdr_seq_grep,
                                     wt_seq_grep = wt_seq_grep)
    result$name = region$name
    result$replicate = replicate
    result$type = type
    replicate_grep_analyses[[i]] = result
  }

  counts.df = bind_rows(replicate_grep_analyses) %>%
    dplyr::select(name, replicate, type, num_reads, everything())
  counts.df$HDR_WT_ratio = counts.df$num_hdr_reads / counts.df$num_wt_reads
  counts.df$HDR_frac = counts.df$num_hdr_reads / counts.df$num_reads
  counts.df$WT_frac = counts.df$num_wt_reads / counts.df$num_reads

  # Summary of stats from the propagation of errors method
  hdr_ratio_res = NULL
  counts.gDNA = counts.df %>% filter(type == "gDNA")
  counts.cDNA = counts.df %>% filter(type == "cDNA")
  if (nrow(counts.gDNA) < 2 | nrow(counts.cDNA) < 2) {
    warning("Unable to calculate stats with < 2 replicates")
  } else {
    hdr_ratio_res = udp_ratio_estimate(counts.df, replicates.df, numerator = "num_hdr_reads", denominator = "num_wt_reads")
    hdr_ratio_res = hdr_ratio_res %>% mutate(name = region$name) %>% select(name, everything())
  }

  return( list(type = "grep_analysis",
               opts = opts,
               replicate_stats = counts.df,
               region_stats = hdr_ratio_res,
               hdr_seq_grep = hdr_seq_grep,
               wt_seq_grep = wt_seq_grep) )
}


# Gets allele counts with grep for a single replicate
replicate_grep_analysis = function(read_seqs, hdr_seq_grep, wt_seq_grep) {
  read_seqs_char = as.character(read_seqs)
  hdr_read_count = sum(sapply(read_seqs_char, FUN = function(s) grepl(hdr_seq_grep, s, fixed=T)))
  wt_read_count = sum(sapply(read_seqs_char, FUN = function(s) grepl(wt_seq_grep, s, fixed=T)))
  return(list(num_reads = length(read_seqs_char),
              num_hdr_reads = hdr_read_count,
              num_wt_reads = wt_read_count))
}


#' Alignment-based GenIE analysis
#'
#' For each replicate associated with an input region, \code{deletion_analysis}
#' identifies HDR or WT sequences, as well as deletion alleles, and returns
#' statistics that indicate for every allele whether the allele:WT ratio
#' differs in cDNA and gDNA. Statistics are also computed for deletion alleles
#' aggregated together.
#'
#' @param regions A data.frame defining GenIE regions.
#' @param replicates A data.frame defining GenIE replicates.
#' @param required_match_left The length of sequence to the left of the HDR site that must exactly match to identify HDR or WT reads.
#' @param required_match_right The length of sequence to the right of the HDR site that must exactly match to identify HDR or WT reads.
#' @param crispr_del_window The window around the cut site within which any deletion is considered a CRISPR deletion.
#' deletions that do not span the region [cut_site - crispr_del_window, cut_site + crispr_del_window]
#' will be ignored; i.e. such reads can be considered HDR or WT.
#' @param min_mapq The minimum mapping quality for reads to be included in the analysis.
#' @param max_mismatch_frac The maximum fraction of mismatches a read can have and be included in the analysis.
#' @param min_aligned_bases The minimum number of aligned bases within the region of interest for a read to be included in the analysis.
#' @param exclude_multiple_deletions If TRUE, then reads with multiple deletions will be excluded from the analysis.
#' @param allele_profile If TRUE, then the result object will contain data.frames named site_profiles and mismatch_profiles, as detailed in the description below.
#' @param del_span_start An integer that specifies the start of a window, relative to the region's highlight site, within which deletions are counted.
#' @param del_span_end An integer that specifies the end of a window, relative to the region's highlight site, within which deletions are counted.
#' @param quiet If TRUE, then no messages are printing during the analysis.
#'
#' @section Details: For a deletion analysis, the regions parameter is a data.frame
#'   with a format as follows. All of the column names below must be specified.
#'   \tabular{lccccclll}{ name \tab sequence_name \tab start \tab end \tab
#'   highlight_site \tab cut_site \tab hdr_allele_profile \tab wt_allele_profile
#'   \tab ref_sequence \cr MUL1_rs6700034  \tab MUL1  \tab 1 \tab 21 \tab 11
#'   \tab 9 \tab ----------A---------- \tab ----------C---------- \tab
#'   ACCGCACCCCCCCGGCCTAAC \cr }
#'
#' \tabular{rl}{
#'   \strong{name} \tab A unique identifier for the region. \cr
#'   \strong{sequence_name} \tab The chromosome or amplicon sequence name. \cr
#'   \strong{start} \tab The start coordinate of the amplicon relative to the chromosome or amplicon reference. \cr
#'   \strong{end} \tab The end coordinate of the amplicon relative to the chromosome or amplicon reference (the end coordinate is included in the region). \cr
#'   \strong{highlight_site} \tab The relative position of the site of interest, usually the HDR SNP site. \cr
#'   \strong{cut_site} \tab The relative position of the cut site. \cr
#'   \strong{ref_sequence} \tab The reference sequence for the amplicon region, which must have length (end - start + 1). \cr
#'   \strong{hdr_allele_profile} \tab An allele profile describing the HDR allele. See details below. \cr
#'   \strong{wt_allele_profile} \tab An allele profile describing the WT allele. See details below. \cr
#' }
#'
#' If multiple rows are defined for `regions`, then a separate analysis is run
#' for each region, using matched replicates from `replicates`.
#'
#' The allele_profile columns indicate the positions in the amplicon sequence
#' that must match a given nucleotide for a read to be considered either HDR
#' or WT. This sequence must be the same length as the reference sequence,
#' and all other positions should be "-". The total sequence region that must
#' match is determined by both the positions of specified nucleotides and by
#' the \strong{required_match_left} and \strong{required_match_right}
#' parameters. These parameters give the length of sequence which must match
#' the provided reference sequence to the left of the leftmost specified
#' nucleotide, or to the right of the rightmost specified nucleotide.
#'
#'
#' The replicates parameter is a data.frame with a format as below.
#' \tabular{lccl}{ name \tab replicate \tab type \tab bam \cr MUL1_rs6700034
#' \tab c1.2 \tab cDNA \tab
#' bam_amplicon/MUL1_rs6700034_cDNA_rep1_pcr2.sortedByCoord.bam \cr
#' MUL1_rs6700034 \tab c1.3 \tab cDNA \tab
#' bam_amplicon/MUL1_rs6700034_cDNA_rep1_pcr3.sortedByCoord.bam \cr }
#'
#' \tabular{rl}{
#'   \strong{name} \tab Indicates the region that a given replicate
#'   corresponds with. All replicates matching the name in the regions table
#'   will be used. \cr
#'   \strong{replicate} \tab an ID for the replicate, which must
#'   be unique among replicates for the region. \cr
#'   \strong{type} \tab Must have the value "cDNA" or "gDNA", indicating
#'   whether a given replicate contains data for cDNA or gDNA. \cr
#'   \strong{bam} \tab the path (relative to the working directory) to
#'   a BAM file with sequencing reads for the replicate. \cr
#' }
#'
#' Statistics can only be computed if there are at least 2 replicates of each
#' type (cDNA and gDNA). Replicates are matched to the region based on the
#' \strong{name} column.
#'
#' @section Results: The returned object is a list, where each item is the
#' result for one region. The result for a region (e.g. results[[1]]) is itself
#' a list, with the following items:
#' \tabular{rl}{
#'   \strong{region_stats} \tab Main analysis output, with statistics indicating whether the HDR/WT levels differ in cDNA relative to gDNA. \cr
#'   \strong{replicate_stats} \tab A data.frame with a row for each replicate, which has counts of reads in different categories and some summary values. \cr
#'   \strong{region} \tab 	Details of the input region the result corresponds to. \cr
#'   \strong{replicates} \tab Details of the input replicates the result corresponds to. \cr
#'   \strong{opts} \tab A list containing the options that were given for the analysis. \cr
#'   \strong{type} \tab Has the value “deletion_analysis”, and is used by plotting functions that take a full del_result list as input. \cr
#' }
#'
#' The main output of interest is the `region_stats` field, which is a one-row data.frame with the following values:
#'
#' \tabular{rl}{
#'   \strong{name} \tab Name of the region. \cr
#'   \strong{hdr_rate_gDNA} \tab Fraction of reads in gDNA identified as HDR. \cr
#'   \strong{hdr_rate_cDNA} \tab Fraction of reads in cDNA identified as HDR. \cr
#'   \strong{wt_rate_gDNA} \tab Fraction of reads in gDNA identified as WT \cr
#'   \strong{wt_rate_cDNA} \tab Fraction of reads in cDNA identified as WT. \cr
#'   \strong{del_rate_gDNA} \tab Fraction of reads in gDNA identified as having a CRISPR deletion. \cr
#'   \strong{del_rate_cDNA} \tab Fraction of reads in cDNA identified as having a CRISPR deletion. \cr
#'   \strong{hdr_effect} \tab HDR allele: Estimated effect size - the amount by which the HDR:WT ratio differs in cDNA relative to gDNA. \cr
#'   \strong{hdr_effect_sd} \tab HDR allele: Standard deviation of the estimated effect size. \cr
#'   \strong{hdr_effect_confint_lo} \tab HDR allele: Lower bound of the 95\% confidence interval for the effect size. \cr
#'   \strong{hdr_effect_confint_hi} \tab HDR allele: Upper bound of the 95\% confidence interval for the effect size. \cr
#'   \strong{hdr_df_estimated} \tab THDR allele: he degrees of freedom used in the unequal variances t-test, which is estimated from the data. \cr
#'   \strong{hdr_pval} \tab HDR allele: A p value from the unequal variants t-test. \cr
#' }
#'
#' The fields above beginning with `hdr_` give statistics relating to the HDR allele.
#' There are 6 equivalent fields that begin with `del_` which relate to all deletions.
#' There are also 6 equivalent fields that begin with `del_window_`, which relate to
#' deletions that are contained within a window around the `highlight_site` (defined in
#' the region input), the extent of which is determined by `crispr_del_window` parameter.
#'
#' @section Additional fields: The deletion_analysis result object additionally has the following fields:
#' \tabular{rl}{
#'   \strong{replicate_qc} \tab A data.frame of summary information for each replicate,
#'   with counts of reads in different categories, editing rates per replicate, and
#'   quality control summary information. \cr
#'   \strong{replicate_alleles} \tab A data.frame of summary information for each unique
#'   allele in each replicate. \cr
#'   \strong{region_alleles} \tab A data.frame of summary information for each unique
#'   allele averaged across all replicates. \cr
#'   \strong{replicate_allele_fractions} \tab A data.frame of summary information for
#'   the top 20 alleles across all replicates, which is used for replicate quality control. \cr
#'   \strong{allele_effect} \tab A data.frame of GenIE effect size estimates (difference
#'   in allele:WT ratio in cDNA vs. gDNA) for each unique allele with sufficient reads
#'   across replicates. \cr
#'   \strong{site_profiles} \tab A data.frame which indicates read counts for each combination of observed
#'   nucleotides at each specified position in the wt_allele_profile, separately for each replicate. For
#'   example, if 1 position is defined, then there will be up to 5 unique combinations in each replicate,
#'   accounting for A, C, G, T, or * (deletion). If 2 positions are defined, there are up to 25 combinations.
#'   This can be useful to check whether the fraction of WT or edited reads is as expected. In particular,
#'   if a heterozygous site is targeted, and only one of the haplotype alleles actually gets edited, then
#'   you expect depletion of only one of the nucleotides at the SNP site. \cr
#'   \strong{mismatch_profiles} \tab A data.frame which contains every unique read profile,
#'   including mismatches, for reads considered as WT or HDR. This can help to identify if
#'   there are an abundance of "HDR" or "WT" reads with mismatches at a particular position,
#'   or an excess of mismatches in general. \cr
#'   \strong{replicate_list} \tab Internal data for each replicate. Not likely to be of interest to the user. \cr
#'   \strong{type} \tab  \cr
#' }
#'
#' @examples
#' del_results = deletion_analysis(regions, replicates)
#' deletion_plots(del_results[[1]])
#' @seealso \code{\link{grep_analysis}}
#' @seealso \code{\link{deletion_plots}}
#' @seealso \code{\link{deletion_summary_plot}}
#' @seealso \code{\link{experiment_summary_plot}}
#' @seealso \code{\link{deletion_alleles_plot}}
#' @seealso \code{\link{deletion_profile_plot}}
#' @seealso \code{\link{replicate_summary_plot}}
#' @seealso \code{\link{replicate_qc_plot}}
#' @seealso \code{\link{allele_effect_plot}}
#' @seealso \code{\link{get_variance_components}}
#' @seealso \code{\link{variance_components_plot}}
#' @seealso \code{\link{power_plots}}
#' @seealso \code{\link{bind_results}}
#' @export
#'
deletion_analysis = function(regions,
                             replicates,
                             required_match_left = 10,
                             required_match_right = 10,
                             crispr_del_window = 100,
                             min_mapq = 0,
                             max_mismatch_frac = 0.05,
                             min_aligned_bases = 50,
                             exclude_multiple_deletions = F,
                             exclude_nonspanning_reads = T,
                             allele_profile = F,
                             del_span_start = -20,
                             del_span_end = 20,
                             quiet = FALSE)
{
  check_regions_and_replicates(regions, replicates)
  replicates$replicate = as.character(replicates$replicate)
  region_names = unique(regions$name)

  opts = list(regions = regions,
              replicates = replicates,
              required_match_left = required_match_left,
              required_match_right = required_match_right,
              crispr_del_window = crispr_del_window,
              min_mapq = min_mapq,
              max_mismatch_frac = max_mismatch_frac,
              min_aligned_bases = min_aligned_bases,
              exclude_multiple_deletions = exclude_multiple_deletions,
              exclude_nonspanning_reads = exclude_nonspanning_reads,
              allele_profile = allele_profile,
              del_span_start = del_span_start,
              del_span_end = del_span_end,
              quiet = quiet,
              qc_max_alleles = 20,  # Not currently exposed
              qc_min_avg_allele_fraction = 0.005,  # Not currently exposed
              qc_exclude_wt = T)  # Not currently exposed
  opts = check_del_opts(opts)

  #library(profvis)
  #profvis({
  #})
  del_results = list()
  region_index = 1
  for (region_name in region_names) {
    cur_region = regions %>% dplyr::filter(name == region_name)
    cur_replicates = replicates %>% dplyr::filter(name == region_name)
    if (nrow(cur_region) != 1) {
      stop(sprintf("delAnalysis: Expected a single region row for region %s.", cur_region$name))
    }
    cur_region = as.list(cur_region[1,])
    if (nrow(cur_replicates) == 0) {
      warning(sprintf("No replicates found for region %s", region_name))
      next
    }

    delResult = region_del_analysis(region = cur_region, replicates = cur_replicates, opts)
    delResult$region = cur_region
    delResult$replicates = cur_replicates

    del_results[[region_index]] = delResult
    region_index = region_index + 1
  }

  return( del_results )
}


# Runs a GenIE deletion analysis for a single region and its replicates
region_del_analysis = function(region, replicates, opts) {
  opts$region = region
  opts$replicates = replicates

  hdr_profile = str_to_upper(region$hdr_allele_profile)
  wt_profile = str_to_upper(region$wt_allele_profile)
  ref_sequence = str_to_upper(region$ref_sequence)
  sequence_name = region$sequence_name

  if (region$end < 1 || region$start < 1) {
    stop(sprintf("ERROR: region '%s' start and end coords (%d, %d) should be positive integers.", region$name, region$start, region$end))
  }
  region_length = region$end - region$start + 1
  if (region_length < 1) {
    stop(sprintf("ERROR: region '%s' length (start - end: %d) should be a positive integer.", region$name, region_length))
  }
  if (region$highlight_site < 1 || region$highlight_site > region_length) {
    stop(sprintf("ERROR: region '%'s highlight_site (%d) should be an integer between 1 and the length of the region.", region$name, region$highlight_site))
  }
  if (region$cut_site < 1 || region$cut_site > region_length) {
    stop(sprintf("ERROR: region '%'s cut_site (%d) should be an integer between 1 and the length of the region.", region$name, region$cut_site))
  }
  sites = list(start = region$start,
               end = region$end,
               highlight_site = region$highlight_site,
               cut_site = region$cut_site)

  replicate_del_analyses = list()
  for (i in 1:nrow(replicates)) {
    bam_file = replicates$bam[i]
    replicate = replicates$replicate[i]
    type = replicates$type[i]

    output(opts, sprintf("\n\nAnalysing region %s, replicate %s, file %s\n", region$name, replicate, bam_file, region$sequence_name))
    output(opts, sprintf("HDR allele: %s\n", hdr_profile))
    output(opts, sprintf("WT allele:  %s\n", wt_profile))
    output(opts, sprintf("REF sequence: %s\n", ref_sequence))

    read_data = get_region_reads_from_bam(region$name, replicate, bam_file, sequence_name, region$start, region$end, ref_sequence, opts$min_mapq)
    if (length(read_data$seq) <= 0) {
      stop(sprintf("\nERROR: No reads retrieved at the specified region for replicate %s from file %s\n", replicate, bam_file))
    }
    aligned_read_data = get_aligned_reads(read_data, ref_sequence, sites$start, sites$end)

    result = replicate_del_analysis(name = region$name,
                                    replicate = replicate,
                                    type = type,
                                    sites = sites,
                                    hdr_profile = hdr_profile,
                                    wt_profile = wt_profile,
                                    ref_sequence = ref_sequence,
                                    read_data = aligned_read_data,
                                    opts = opts)
    result$name = region$name
    result$replicate = replicate
    result$type = type
    replicate_del_analyses[[i]] = result
  }
  # Merge together results from all replicates
  replicate_counts = bind_rows(lapply(replicate_del_analyses, FUN = function(res) res$counts))
  replicate.udp.df = bind_rows(lapply(replicate_del_analyses, FUN = function(res) res$udp.df))

  # Determine which UDPs are shared between cDNA and gDNA
  gDNACounts = replicate.udp.df %>% filter(type == "gDNA") %>% group_by(udp) %>%
    summarise(udpcount_gDNA = sum(num_reads))
  cDNACounts = replicate.udp.df %>% filter(type == "cDNA") %>% group_by(udp) %>%
    summarise(udpcount_cDNA = sum(num_reads))
  replicate.udp.df = replicate.udp.df %>%
    left_join(gDNACounts, by=c("udp")) %>%
    left_join(cDNACounts, by=c("udp"))
  threshold = 1
  replicate.udp.df$udp_sharing = "unclear"
  replicate.udp.df$udp_sharing[replicate.udp.df$udpcount_gDNA >= threshold & replicate.udp.df$udpcount_cDNA >= threshold] = "both"
  replicate.udp.df$udp_sharing[replicate.udp.df$udpcount_gDNA > 0 & is.na(replicate.udp.df$udpcount_cDNA)] = "gDNA only"
  replicate.udp.df$udp_sharing[replicate.udp.df$udpcount_cDNA > 0 & is.na(replicate.udp.df$udpcount_gDNA)] = "cDNA only"

  stats_res = get_region_del_stats(replicate_counts, replicates)
  stats_res$region_summary$name = region$name

  qc_metrics = replicate_qc_metrics(stats_res$replicate_stats, replicate.udp.df,
                                    opts$qc_max_alleles, opts$qc_min_avg_allele_fraction, opts$qc_exclude_wt)
  replicate_stats = qc_metrics$replicate_qc
  #replicate_stats = stats_res$replicate_stats %>%
  #  left_join(qc_metrics$replicate_qc, by=c("replicate", "type"))

  merged.udp.df = summarise(replicate.udp.df %>% group_by(name, type, udp),
                            num_reads = sum(num_reads),
                            is_hdr_allele = first(is_hdr_allele),
                            is_wt_allele = first(is_wt_allele),
                            has_any_deletion = first(has_any_deletion),
                            has_crispr_deletion = first(has_crispr_deletion),
                            has_deletion_in_window = first(has_deletion_in_window),
                            deletion_start = first(deletion_start),
                            deletion_end = first(deletion_end),
                            deletion2_start = first(deletion2_start),
                            deletion2_end = first(deletion2_end),
                            avg_seq_length = mean(avg_seq_length),
                            avg_mismatch_count = mean(avg_mismatch_count),
                            udp_sharing = first(udp_sharing)) %>%
    arrange(-num_reads) %>% ungroup()

  uns.df = get_uns_data(replicate.udp.df, replicates.df = replicates, region_name = region$name)
  if (!is.null(uns.df)) {
    uns.df = uns.df %>% dplyr::mutate(name = region$name) %>% dplyr::select(name, everything())
  }

  # Make a table of the "site profiles" - combinations of alleles at sites of interest
  site.profiles.df = NULL
  wt_hdr.df = NULL
  if (opts$allele_profile) {
    replicate_wt_hdr = lapply(replicate_del_analyses, FUN = function(res) res$wt_hdr.df)
    wt_hdr.df = bind_rows(replicate_wt_hdr)

    replicate_site_profiles = lapply(replicate_del_analyses, FUN = function(res) res$site.profiles.df)
    site.profiles.df = bind_rows(replicate_site_profiles)
    site.profiles.df$replicate = as.character(site.profiles.df$replicate)
    # Combine replicates and add merged counts as a replicate named "all"
    site.profile.merged = site.profiles.df %>% group_by(sites_profile, type) %>%
      dplyr::summarise(name = first(name),
                       replicate = "all",
                       count = sum(count)) %>%
      dplyr::select(name, replicate, type, sites_profile, count)
    site.profiles.df = bind_rows(site.profiles.df, site.profile.merged) %>%
      dplyr::arrange(replicate, type, sites_profile)
  }

  result_list = list(type = "deletion_analysis",
                     opts = opts,
                     replicate_list = replicate_del_analyses,
                     replicate_qc = qc_metrics$replicate_qc,
                     replicate_allele_fractions = qc_metrics$replicate_allele_fractions,
                     replicate_alleles = replicate.udp.df,
                     region_alleles = merged.udp.df,
                     allele_effect = uns.df,
                     site_profiles = site.profiles.df,
                     mismatch_profiles = wt_hdr.df,
                     replicate_stats = replicate_stats,
                     region_stats = stats_res$region_summary)
  return(result_list)
}

replicate_del_analysis = function(name, replicate, type, sites, hdr_profile , wt_profile, ref_sequence, read_data, opts)
{
  counts = list(name = name, replicate = replicate, type = type)

  region_length = (sites$end - sites$start + 1)
  if (nchar(wt_profile) != region_length) {
    stop(sprintf("The WT allele profile given has length %d, but the region size (end - start + 1) is %d", nchar(wt_profile), region_length))
  }
  if (is.na(hdr_profile)) {
    hdr_profile = ""
  }

  replicate_name = paste(name, replicate, type, sep="_")
  reads.df = read_data$alignedReads

  output(opts, sprintf("%d starting cDNA reads\n", sum(reads.df$count)))
  output(opts, sprintf("%d cDNA reads of %d total (%.1f%%) were soft-clipped\n", read_data$num_softclipped, read_data$num_reads, 100.0 * read_data$num_softclipped / read_data$num_reads))
  output(opts, sprintf("%d cDNA reads of %d total (%.1f%%) were hard-clipped\n", read_data$num_hardclipped, read_data$num_reads, 100.0 * read_data$num_hardclipped / read_data$num_reads))
  output(opts, sprintf("%d cDNA reads of %d total (%.1f%%) had insertions\n", read_data$num_insertion, read_data$num_reads, 100.0 * read_data$num_insertion / read_data$num_reads))
  output(opts, sprintf("%d cDNA reads of %d total (%.1f%%) were likely primer dimers\n", read_data$num_primerdimer, read_data$num_reads, 100.0 * read_data$num_primerdimer / read_data$num_reads))
  counts$num_softclipped = read_data$num_softclipped
  counts$num_hardclipped = read_data$num_hardclipped
  counts$num_insertion = read_data$num_insertion
  counts$num_primerdimer = read_data$num_primerdimer
  counts$num_reads = read_data$num_reads

  hdr_profile_chars = strsplit(hdr_profile, "")[[1]]
  wt_profile_chars = strsplit(wt_profile, "")[[1]]
  ref_seq_chars = strsplit(ref_sequence, "")[[1]]

  reads.df$read_chars = sapply(reads.df$region_read, FUN=function(s) strsplit(s, ""))
  reads.df$seq_length = sapply(reads.df$read_chars, FUN=function(s) sum(is_dna_letter(s)))
  #reads.df$seq_length = sapply(reads.df$read_chars, FUN=function(s) sum(s %in% c("A", "C", "G", "T")))  # SLOWER
  #reads.df$seq_length = sapply(reads.df$region_read, FUN=function(s) str_count(s, "[ACGT]"))  # SLOWER
  #reads.df$seq_length = sapply(reads.df$region_read, FUN=function(s) nchar(gsub("[-*]", "", s)))  # SLOWER

  reads.df$mismatch_count = sapply(reads.df$read_chars, FUN=get_mismatch_chars_count, ref_seq_chars)
  #reads.df$mismatch_count = sapply(reads.df$read_chars, get_mismatch_chars_count, ref_seq_chars)  # SLOWER
  #reads.df$mismatch_count = sapply(reads.df$region_read, get_mismatch_count, ref_sequence)  # SLOWER

  reads.df$spanning_read = T
  span_site = sites$cut_site
  if (!is.na(sites$highlight_site)) {
    span_site = sites$highlight_site
  }
  reads.df$spanning_read = sapply(reads.df$region_read, FUN=function(s) (substring(s, span_site, span_site) != '-'))
  #reads.df$spanning_read = sapply(reads.df$read_chars, FUN=function(s) (s[span_site] != '-'))
  if (opts$exclude_nonspanning_reads) {
    counts$reads_excluded_nonspanning = sum((!reads.df$spanning_read) * reads.df$count)
    output(opts, sprintf("%d of %d reads (%.2f%%) excluded due to not spanning the site of interest (position %d)\n",
                counts$reads_excluded_nonspanning, counts$num_reads, 100.0 * counts$reads_excluded_nonspanning / counts$num_reads,
                span_site))
    reads.df = reads.df[reads.df$spanning_read, ]
  }

  # Exclude reads that don't cover enough of the region of interest
  exclude_for_overlap = (reads.df$seq_length < opts$min_aligned_bases)
  counts$reads_excluded_for_minoverlap = sum(exclude_for_overlap * reads.df$count)
  output(opts, sprintf("%d of %d reads (%.2f%%) excluded due to aligning to less than %d bp in the region of interest\n",
              counts$reads_excluded_for_minoverlap, counts$num_reads, 100.0 * counts$reads_excluded_for_minoverlap / counts$num_reads,
              opts$min_aligned_bases))
  reads.df = reads.df[!exclude_for_overlap, ]

  # Identify unique deletion profile (UDP) for each read
  #reads.df$udp = get_read_udps(reads.df$region_read, wt_profile_chars)
  reads.df$udp = sapply(reads.df$read_chars, FUN=get_read_chars_udp, wt_profile_chars)
  #reads.df$udp = sapply(reads.df$region_read, FUN=get_read_udp, wt_profile_chars) # SLOWER

  reads.df$is_wt_allele = (reads.df$udp == wt_profile)
  # make sure we only count as WT those reads which actually cover the HDR site
  reads.df$is_wt_allele[!reads.df$spanning_read] = NA

  reads.df$has_any_deletion = sapply(reads.df$udp, FUN=function(s) grepl("[*]", s))
  if (!is.na(opts$crispr_del_window)) {
    startPos = max(1, sites$cut_site - opts$crispr_del_window)
    endPos = min(region_length, sites$cut_site + opts$crispr_del_window)
    reads.df$has_crispr_deletion = (!reads.df$is_wt_allele & grepl("[*]", substr(reads.df$udp, startPos, endPos)) )

    # We change the UDPs to "zero out" any deletions that don't overlap
    # with the "editing region" defined by the cut site and crispr_del_window
    updateUDPEdits = function(udp, cut_site, crispr_del_window) {
      dels = str_locate_all(udp, "[*]+")[[1]]
      if (nrow(dels) == 0)
        return(udp)
      for (i in 1:nrow(dels)) {
        start = dels[i,1]
        end = dels[i,2]
        if (start > cut_site + crispr_del_window | end < cut_site - crispr_del_window) {
          substr(udp, start, end) <- strrep("-", end - start + 1)
        }
      }
      return(udp)
    }
    if (sum(reads.df$has_any_deletion) > 0) {
      reads.df$udp[reads.df$has_any_deletion & !reads.df$is_wt_allele] = sapply(reads.df$udp[reads.df$has_any_deletion & !reads.df$is_wt_allele], updateUDPEdits, sites$cut_site, opts$crispr_del_window)
      # If we have zeroed out a deletion, then it's important that we also zero out
      # the flag for whether it has any deletion.
      reads.df$has_any_deletion[reads.df$has_any_deletion] = sapply(reads.df$udp[reads.df$has_any_deletion], FUN=function(s) grepl("[*]", s))
    }
  } else {
    # Accept deletions anywhere in the read
    reads.df$has_crispr_deletion = reads.df$has_any_deletion
  }
  # Include deletions anywhere in the read
  reads.df$has_multiple_deletions = reads.df$has_any_deletion
  if (sum(reads.df$has_multiple_deletions) > 0) {
    reads.df$has_multiple_deletions[reads.df$has_multiple_deletions] =
      sapply(reads.df$udp[reads.df$has_multiple_deletions], FUN=function(s) grepl("[*]+[^*]+[*]+", s))
  }

  # Do this again, to count as WT any reads which had non-CRISPR deletions zeroed out
  reads.df$is_wt_allele = (reads.df$udp == wt_profile)
  # make sure we only count as WT those reads which actually cover the HDR site
  reads.df$is_wt_allele[!reads.df$spanning_read] = NA
  # If the WT profile has a deletion, then these should not be counted as CRISPR deletions
  reads.df$has_crispr_deletion[reads.df$is_wt_allele] = F

  n_wt_vars = sum(is_dna_letter(wt_profile_chars) & wt_profile_chars != ref_seq_chars)
  if (n_wt_vars > 0) {
    reads.df$mismatch_count[reads.df$is_wt_allele] = reads.df$mismatch_count[reads.df$is_wt_allele] - n_wt_vars
  }

  # Identify which reads are HDR.
  reads.df$is_hdr_allele = F
  if (hdr_profile != "") {
    if (nchar(hdr_profile) != region_length) {
      stop(sprintf("The HDR allele profile given has length %d, but the region size (end - start + 1) is %d", nchar(hdr_profile), region_length))
    }
    # Check that the HDR profile and WT profile have DNA characters at
    # the same positions
    #
    #if (any((hdr_profile_chars == "-") != (wt_profile_chars == "-"))) {
    #  stop("Error: HDR profile and WT profile should both indicate the expected sequence letters at the same positions")
    #}
    if (grepl("B|D|H|V", hdr_profile)) {
      hdr_profile_set = get_hdr_profiles(hdr_profile)
      reads.df$is_hdr_allele = sapply(reads.df$udp, function(udp) any(udp == hdr_profile_set))
    } else {
      reads.df$is_hdr_allele = (reads.df$udp == hdr_profile)
    }

    if (!any(hdr_profile_chars == "*")) {
      # If the HDR profile doesn't have a deletion itself, then we don't expect any HDR alleles to have
      # deletions near the CRISPR cut site
      reads.df$is_hdr_allele = reads.df$is_hdr_allele & !reads.df$has_crispr_deletion
    }
    # We don't want to count the HDR site itself as a mismatch
    n_hdr_vars = sum(is_dna_letter(hdr_profile_chars) & hdr_profile_chars != ref_seq_chars)
    if (n_hdr_vars > 0) {
      reads.df$mismatch_count[reads.df$is_hdr_allele] = reads.df$mismatch_count[reads.df$is_hdr_allele] - n_hdr_vars
    }
  }
  if (any(reads.df$mismatch_count < 0)) {
    stop("ERROR: something went wrong - got a negative mismatch count.")
  }

  # Exclude reads with too many mismatches
  exclude_for_mismatches = (reads.df$mismatch_count / reads.df$seq_length) > opts$max_mismatch_frac
  counts$reads_excluded_for_mismatches = sum(exclude_for_mismatches * reads.df$count)
  output(opts, sprintf("%d of %d reads (%.2f%%) excluded due to having more than %.2f%% of the read being mismatches\n",
              counts$reads_excluded_for_mismatches, counts$num_reads, 100.0 * counts$reads_excluded_for_mismatches / counts$num_reads,
              opts$max_mismatch_frac * 100))
  reads.df = reads.df[!exclude_for_mismatches, ]

  n_wt_reads = sum(reads.df$is_wt_allele * reads.df$count, na.rm = T)
  if (n_wt_reads < 1) {
    warning("ERROR: no wild-type reads found. Check the WT allele profile. You may also want to check that your reads span the edit site, given your read length and amplicon coords.")
  } else if (n_wt_reads < 100) {
    warning(sprintf("Warning: only %d wild-type reads found in experiment.", n_wt_reads))
  }

  counts$reads_excluded_for_multiple_deletions = 0
  if (opts$exclude_multiple_deletions) {
    counts$reads_excluded_for_multiple_deletions = sum(reads.df$num_reads[reads.df$has_multiple_deletions])
    counts$udps_excluded_for_multiple_deletions = sum(reads.df$has_multiple_deletions)
    output(opts, sprintf("%d of %d UDPs (%d of %d reads, %.2f%%) excluded due to having more than one separate deletion\n",
                counts$udps_excluded_for_multiple_deletions, counts$num_udps,
                counts$reads_excluded_for_multiple_deletions, counts$num_reads, 100.0 * counts$reads_excluded_for_multiple_deletions / counts$num_reads))
    reads.df = reads.df[!reads.df$has_multiple_deletions, ]
  }

  counts$num_wt_reads = n_wt_reads
  counts$num_hdr_reads = sum(reads.df$is_hdr_allele * reads.df$count, na.rm = T)
  counts$num_deletion_reads = sum((reads.df$has_crispr_deletion & !reads.df$is_hdr_allele) * reads.df$count, na.rm = T)
  counts$num_kept_reads = sum(reads.df$count, na.rm = T)

  # Count deletion reads where the deletion is contained within a specific window
  window_start = min(region_length, max(1, span_site + opts$del_span_start))
  window_end = min(region_length, max(1, span_site + opts$del_span_end))
  reads.df$has_deletion_in_window = ( grepl("[*]", substr(reads.df$udp, window_start, window_end)) &
                                     !grepl("[*]", substr(reads.df$udp, 1, window_start)) &
                                     !grepl("[*]", substr(reads.df$udp, window_end, region_length)) )
  counts$num_deletion_reads_in_window = sum((reads.df$has_deletion_in_window & !reads.df$is_hdr_allele) * reads.df$count, na.rm = T)

  # Aggregate reads according to their UDP
  udp.df = summarise(reads.df %>% group_by(udp),
                     num_reads = sum(count),
                     is_hdr_allele = first(is_hdr_allele),
                     is_wt_allele = first(is_wt_allele),
                     spanning_read = first(spanning_read),
                     has_any_deletion = first(has_any_deletion),
                     has_crispr_deletion = first(has_crispr_deletion),
                     has_deletion_in_window = first(has_deletion_in_window),
                     has_multiple_deletions = first(has_multiple_deletions),
                     avg_seq_length = mean(seq_length),
                     avg_mismatch_count = mean(mismatch_count)) %>%
    mutate(name = name, replicate = replicate, type = type) %>%
    select(name, replicate, type, everything()) %>%
    arrange(-num_reads)
  counts$num_udps = nrow(udp.df)

  udp.df$has_edit = udp.df$has_crispr_deletion | udp.df$is_hdr_allele
  counts$num_edit_reads = sum(udp.df$num_reads[udp.df$has_edit], na.rm = T)

  # Order UDPs by the deletion start position and plot
  get_deletion_loc = function(i) {
    if (udp.df$has_crispr_deletion[i]) {
      return(sapply(udp.df$udp[i], FUN=function(udp) str_locate_all(udp, "[\\*]+")))
    }
    return(NA)
  }
  udp_dels = sapply(1:nrow(udp.df), FUN=get_deletion_loc)
  udp.df$deletion_start = sapply(udp_dels, FUN=function(x) {if (length(x) == 1) return(NA); x[1,1]})
  udp.df$deletion_end = sapply(udp_dels, FUN=function(x) {if (length(x) == 1) return(NA); x[1,2] + 1})
  udp.df$deletion2_start = sapply(udp_dels, FUN=function(x) {if (length(x) == 1) return(NA); if(nrow(x) > 1) {x[2,1]} else {NA}})
  udp.df$deletion2_end = sapply(udp_dels, FUN=function(x) {if (length(x) == 1) return(NA); if(nrow(x) > 1) {x[2,2]+1} else {NA}})
  #udp.df$deletion_length = sapply(udp_dels, FUN=function(x) x[2]-x[1]+1)

  # Make a table which has just the different versions of the WT and HDR alleles
  wt_hdr.df = NULL
  site.profiles.df = NULL
  if (opts$allele_profile) {
    wt_hdr.reads.df = reads.df %>% dplyr::filter(is_wt_allele | is_hdr_allele)
    #wt_hdr.reads.df$mismatch_profile = get_read_mismatch_profile(wt_hdr.reads.df$region_read, ref_seq_chars)
    wt_hdr.reads.df$mismatch_profile = sapply(wt_hdr.reads.df$read_chars, FUN=get_read_chars_mismatch_profile, ref_seq_chars)

    wt_hdr.df = summarise(wt_hdr.reads.df %>% group_by(mismatch_profile),
                          num_reads = sum(count),
                          mismatch_count = first(mismatch_count),
                          spanning_read = first(spanning_read),
                          is_wt_allele = first(is_wt_allele),
                          is_hdr_allele = first(is_hdr_allele)) %>%
      mutate(name = name, replicate = replicate, type = type) %>%
      dplyr::select(name, replicate, type, everything()) %>%
      dplyr::arrange(-num_reads)

    # Make a table which has all variations of read sequences at the "profile" sites,
    # i.e. those sites used to identify WT and HDR alleles
    #reads.df$sites_profile = get_read_site_profile(reads.df$region_read, wt_profile_chars)
    profile_positions = which(is_dna_letter(wt_profile_chars))
    reads.df$sites_profile = sapply(reads.df$region_read, FUN=get_read_site_profile, profile_positions)
    site.profiles.df = reads.df %>% group_by(sites_profile) %>%
      summarise(count = sum(count))  %>%
      mutate(name = name, replicate = replicate, type = type) %>%
      dplyr::select(name, replicate, type, everything())
  }
  return(list(udp.df = udp.df,
              wt_hdr.df = wt_hdr.df,
              site.profiles.df = site.profiles.df,
              counts = counts))
}


get_region_del_stats = function(replicate_data, replicates.df) {
  replicate_data$HDR_WT_ratio = (replicate_data$num_hdr_reads / replicate_data$num_wt_reads)
  replicate_data$DEL_WT_ratio = (replicate_data$num_deletion_reads / replicate_data$num_wt_reads)
  replicate_data$HDR_rate = replicate_data$num_hdr_reads / replicate_data$num_kept_reads
  replicate_data$DEL_rate = replicate_data$num_deletion_reads / replicate_data$num_kept_reads
  replicate_data$editing_rate = replicate_data$num_edit_reads / replicate_data$num_kept_reads
  replicate_data$WT_rate = replicate_data$num_wt_reads / replicate_data$num_kept_reads

  # Summary of stats from the propagation of errors method
  hdr_ratio_res = NULL
  del_ratio_res = NULL
  del_ratio_res_window = NULL
  method = ""
  stats.gDNA = replicate_data %>% filter(type == "gDNA")
  stats.cDNA = replicate_data %>% filter(type == "cDNA")

  if (nrow(stats.gDNA) < 2 | nrow(stats.cDNA) < 2) {
    warning("Unable to calculate stats with < 2 replicates")
    stats.summary = NULL
  } else {
    denom = "num_wt_reads"
    hdr_ratio_res = udp_ratio_estimate(replicate_data, replicates.df, numerator = "num_hdr_reads", denominator = denom)
    del_ratio_res = udp_ratio_estimate(replicate_data, replicates.df, numerator = "num_deletion_reads", denominator = denom)
    del_ratio_res_window = udp_ratio_estimate(replicate_data, replicates.df, numerator = "num_deletion_reads_in_window", denominator = denom)

    hdr_wt_res = as.list(hdr_ratio_res)
    names(hdr_wt_res) = paste0("hdr_", names(hdr_wt_res))
    del_wt_res = as.list(del_ratio_res)
    names(del_wt_res) = paste0("del_", names(del_wt_res))
    del_wt_window_res = as.list(del_ratio_res_window)
    names(del_wt_window_res) = paste0("del_window_", names(del_wt_window_res))
    stats.summary = c(hdr_rate_gDNA = mean(stats.gDNA$HDR_rate),
                      hdr_rate_cDNA = mean(stats.cDNA$HDR_rate),
                      del_rate_gDNA = mean(stats.gDNA$DEL_rate),
                      del_rate_cDNA = mean(stats.cDNA$DEL_rate),
                      wt_rate_gDNA = mean(stats.gDNA$WT_rate),
                      wt_rate_cDNA = mean(stats.cDNA$WT_rate),
                      hdr_wt_res,
                      del_wt_res,
                      del_wt_window_res)
  }

  result_list = list(replicate_stats = replicate_data,
                     region_summary = as_tibble(stats.summary))
  return(result_list)
}


#' Returns quality control statistics for results from a deletion analysis.
#'
#' @param del_result A data.frame defining GenIE regions.
#' @param max_alleles Param description
#' @param min_avg_allele_fraction Param description
#' @param exclude_wt Param description
#' @return Returns quality control statistics for results from a deletion analysis.
#'
#' @examples
#' del_results = deletion_analysis(regions, replicates)
#' metrics = replicate_qc_metrics(del_results[[1]])
#' @seealso \code{\link{deletion_analysis}}
#' @seealso \code{\link{replicate_summary_plot}}
#' @seealso \code{\link{replicate_qc_plot}}
#' @export
#'
replicate_qc_metrics = function(replicates, replicate_alleles, max_alleles = 20, min_avg_allele_fraction = 0.005, exclude_wt = FALSE) {
  # Here we get various metrics for each replicate that enable visual or automatic
  # QC. One of the key inputs is, for each replicate, the UDP fraction for the
  # top N UDPs. When plotted, this should enable visually identifying outlier
  # samples based on the variability of the UDP fractions.

  # We need to have a value for each of the top N UDPs. To do this,
  # spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
  # fill = 0 for when a replicate is missing, and then gather back.
  getType = function(type_rep) { sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][1]) }
  getReplicate = function(type_rep) { sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][2]) }
  replicate.udp.filled.df = replicate_alleles %>%
    dplyr::mutate(type_replicate = paste(type, replicate, sep="^")) %>%
    dplyr::select(udp, type_replicate, is_wt_allele, num_reads) %>%
    tidyr::spread(type_replicate, num_reads, fill = 0) %>%
    tidyr::gather(key="type_replicate", value="num_reads", -udp, -is_wt_allele) %>%
    dplyr::mutate(type = getType(type_replicate),
                  replicate = getReplicate(type_replicate)) %>%
    dplyr::select(udp, type, is_wt_allele, replicate, num_reads)

  replicate.totalreads.df = replicate_alleles %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(replicate_total_reads = sum(num_reads))
  udp.totalreads.df = replicate_alleles %>%
    dplyr::group_by(udp) %>%
    dplyr::summarise(udp_total_reads = sum(num_reads),
                     udp_reads_gDNA = sum(num_reads[type == "gDNA"]),
                     udp_reads_cDNA = sum(num_reads[type == "cDNA"]),
                     is_wt_allele = first(is_wt_allele)) %>%
    dplyr::arrange(-udp_total_reads)

  if (exclude_wt) {
    replicate.udp.filled.df = replicate.udp.filled.df %>% dplyr::filter(!is_wt_allele)
    udp.totalreads.df = udp.totalreads.df %>% dplyr::filter(!is_wt_allele)
  }

  if (max_alleles < nrow(udp.totalreads.df)) {
    udp.totalreads.df = udp.totalreads.df[1:max_alleles,]
  }

  replicate.udp.filled.df = replicate.udp.filled.df %>%
    dplyr::inner_join(udp.totalreads.df, by = "udp") %>%
    dplyr::left_join(replicate.totalreads.df, by="replicate") %>%
    dplyr::mutate(udp_fraction = num_reads / replicate_total_reads) %>%
    dplyr::arrange(-udp_total_reads)
  udp.avg.fractions = replicate.udp.filled.df %>%
    dplyr::group_by(udp) %>%
    dplyr::summarise(avg_udp_fraction = mean(udp_fraction),
                     avg_gdna_udp_fraction = mean(udp_fraction[type == "gDNA"]))

  replicate.udp.filled.df = replicate.udp.filled.df %>%
    dplyr::left_join(udp.avg.fractions, by = "udp") %>%
    dplyr::filter(avg_udp_fraction >= min_avg_allele_fraction)

  num_udps = length(unique(replicate.udp.filled.df$udp))
  replicate.udp.filled.df$udp = factor(replicate.udp.filled.df$udp, levels=unique(replicate.udp.filled.df$udp))
  replicate.udp.filled.df$udp_id = factor(as.integer(replicate.udp.filled.df$udp))

  # Compute a statistic which, for each replicate, is the mean deviation of the
  # replicate's UDP fractions from the mean UDP fractions across replicates.
  replicate.udp_avg_deviation.df = replicate.udp.filled.df %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(avg_udp_deviation = mean(abs(udp_fraction - avg_udp_fraction)),
                     type = first(type))

  # Compute the deviation relative to the mean of gDNA only
  replicate.udp_avg_deviation_from_gDNA.df = replicate.udp.filled.df %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(avg_udp_deviation_from_gDNA = mean(abs(udp_fraction - avg_gdna_udp_fraction)),
                     type = first(type))

  # Compute the average RELATIVE deviation (CV) across UDPs for each replicate.
  replicate.udp_rel_deviation_from_gDNA.df = replicate.udp.filled.df %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(avg_udp_rel_deviation_from_gDNA = mean(abs(udp_fraction / avg_gdna_udp_fraction - 1)),
                     type = first(type))

  replicate_qc = replicates %>%
    dplyr::left_join(replicate.udp_avg_deviation.df, by=c("replicate", "type")) %>%
    dplyr::left_join(replicate.udp_avg_deviation_from_gDNA.df, by=c("replicate", "type")) %>%
    dplyr::left_join(replicate.udp_rel_deviation_from_gDNA.df, by=c("replicate", "type"))

  # Use multiple metrics across replicates to calculate outlier scores using KNN
  # and isolation forest, separately for cDNA and gDNA
  n_cDNA = sum(replicates$type == "cDNA")
  n_gDNA = sum(replicates$type == "gDNA")
  replicate_qc$outlier_score_knn = NA
  if (n_cDNA >= 4 && !any(is.na(replicate_qc$avg_udp_deviation))) {
    stats.cDNA.df = replicate_qc %>%
      filter(type == "cDNA") %>%
      select(num_udps, HDR_WT_ratio, DEL_WT_ratio, avg_udp_deviation, avg_udp_deviation_from_gDNA)
    # Compute KNN outlier metric. First scale the data, and for UDP deviation set low values to
    # zero since these are good and shouldn't be considered outliers.
    stats_scaled <- scale(stats.cDNA.df)
    stats_scaled[is.na(stats_scaled)] = 0
    stats_scaled[, "avg_udp_deviation"] = sapply(stats_scaled[, "avg_udp_deviation"], function(x) max(0, x))
    stats_scaled[, "avg_udp_deviation_from_gDNA"] = sapply(stats_scaled[, "avg_udp_deviation_from_gDNA"], function(x) max(0, x))
    stats_nn <- FNN::get.knn(stats_scaled, k = max(2, floor(n_cDNA / 2)))
    stats_nnd <- rowMeans(stats_nn$nn.dist)
    replicate_qc$outlier_score_knn[replicate_qc$type == "cDNA"] = stats_nnd
  }
  if (n_gDNA >= 4 && !any(is.na(replicate_qc$avg_udp_deviation))) {
    stats.gDNA.df = replicate_qc %>%
      filter(type == "gDNA") %>%
      select(num_udps, HDR_WT_ratio, DEL_WT_ratio, avg_udp_deviation, avg_udp_deviation_from_gDNA)
    stats_scaled <- scale(stats.gDNA.df)
    stats_scaled[is.na(stats_scaled)] = 0
    stats_scaled[, "avg_udp_deviation"] = sapply(stats_scaled[, "avg_udp_deviation"], function(x) max(0, x))
    stats_scaled[, "avg_udp_deviation_from_gDNA"] = sapply(stats_scaled[, "avg_udp_deviation_from_gDNA"], function(x) max(0, x))
    stats_nn <- FNN::get.knn(stats_scaled, k = max(2, floor(n_gDNA / 2)))
    stats_nnd <- rowMeans(stats_nn$nn.dist)
    replicate_qc$outlier_score_knn[replicate_qc$type == "gDNA"] = stats_nnd
  }

  return( list(replicate_qc = replicate_qc,
               replicate_allele_fractions = replicate.udp.filled.df) )
}


get_uns_data = function(replicate.udp.df, replicates.df, region_name) {
  udp_types = unique(replicate.udp.df$type)
  if (!("cDNA" %in% udp_types & "gDNA" %in% udp_types)) {
    warning("Cannot compute allele UNS statistics without both cDNA and gDNA replicates")
    return(NULL)
  }

  # We need to have the same number of replicates for every UDP. We do this
  # by spreading the replicates out into columns (cDNA_1, cDNA_2, etc.), with
  # fill = 0 for when a replicate is missing, and then gather back into separate
  # replicates
  getType = function(type_rep) { as.character(sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][1])) }
  getReplicate = function(type_rep) { as.character(sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][2])) }
  replicate.udp.filled.df = replicate.udp.df %>%
    dplyr::mutate(type_replicate = paste(type, replicate, sep="^")) %>%
    dplyr::select(udp, num_reads, is_hdr_allele, is_wt_allele, has_crispr_deletion, deletion_start, deletion_end, type_replicate) %>%
    tidyr::spread(type_replicate, num_reads, fill = 0) %>%
    tidyr::gather(key="type_replicate", value="num_reads", -(udp:deletion_end)) %>%
    dplyr::mutate(type = getType(type_replicate),
                  replicate = getReplicate(type_replicate)) %>%
    dplyr::select(-type_replicate)

  # We use the same function to calculate stats for all UDPs as we do for the HDR allele.
  # For this we need a column with num_wt_reads.
  replicate.read_counts.df = replicate.udp.filled.df %>%
    dplyr::select(replicate, type, is_wt_allele, num_reads) %>%
    dplyr::group_by(replicate, type) %>%
    dplyr::summarise(replicate_num_reads = sum(num_reads),
                     num_wt_reads = sum(num_reads[is_wt_allele])) %>%
    ungroup()

  # Results not likely to be stable if the WT count is too low
  all_wt_gDNA_count = sum(replicate.read_counts.df %>% dplyr::filter(type == "gDNA") %>% .$num_wt_reads)
  all_wt_cDNA_count = sum(replicate.read_counts.df %>% dplyr::filter(type == "cDNA") %>% .$num_wt_reads)
  if (all_wt_gDNA_count < 100) {
    warning(sprintf("%s wild-type gDNA count (%d) is low and may lead to unstable estimates.", region_name, all_wt_gDNA_count))
    if (all_wt_gDNA_count < 1) {
      warning(sprintf("%s wild-type gDNA count is zero!", region_name))
      return(NULL)
    }
    if (all_wt_gDNA_count < min_gDNA_count) {
      warning(sprintf("%s wild-type gDNA count (%d) is below minimum (%d).", region_name, all_wt_gDNA_count, min_gDNA_count))
      return(NULL)
    }
  }
  if (all_wt_cDNA_count < 100) {
    warning(sprintf("%s wild-type cDNA count (%d) is low and may lead to unstable estimates.", region_name, all_wt_cDNA_count))
    if (all_wt_cDNA_count < 1) {
      warning(sprintf("%s wild-type cDNA count is zero!", region_name))
      return(NULL)
    }
  }

  # Get total read counts by UDP, for gDNA and cDNA
  udp.total_counts.df = replicate.udp.filled.df %>%
    group_by(udp) %>%
    dplyr::summarise(gDNA_total_count = sum(num_reads[type == "gDNA"]),
                     cDNA_total_count = sum(num_reads[type == "cDNA"]),
                     total_count = gDNA_total_count + cDNA_total_count)

  # Following this, we have a row per replicate per UDP, with counts for cDNA
  # and gDNA, which are zero if the UDP wasn't observed in the replicate.
  replicate.dels.df = replicate.udp.filled.df %>%
    dplyr::left_join(udp.total_counts.df, by="udp") %>%
    dplyr::left_join(replicate.read_counts.df %>% dplyr::select(-type), by="replicate")
#  dplyr::filter(gDNA_total_count >= min_gDNA_count & cDNA_total_count >= min_cDNA_count) %>%

  # Summarize details per uDP
  udp.dels.df = replicate.dels.df %>%
    group_by(udp) %>%
    dplyr::filter(!is_wt_allele) %>%
    summarise(gDNA_total_count = first(gDNA_total_count),
              cDNA_total_count = first(cDNA_total_count),
              total_count = first(total_count),
              mean_cDNA_count = mean(num_reads[type == "cDNA"], na.rm = T),
              mean_gDNA_count = mean(num_reads[type == "gDNA"], na.rm = T),
              cDNA_ratio = mean(num_reads[type == "cDNA"] / num_wt_reads[type == "cDNA"], na.rm = T),
              sd_cDNA_ratio = sd(num_reads[type == "cDNA"] / num_wt_reads[type == "cDNA"], na.rm = T),
              se_cDNA_ratio = sd_cDNA_ratio / sqrt(n()),
              gDNA_ratio = mean(num_reads[type == "gDNA"] / num_wt_reads[type == "gDNA"], na.rm = T),
              sd_gDNA_ratio = sd(num_reads[type == "gDNA"] / num_wt_reads[type == "gDNA"], na.rm = T),
              se_gDNA_ratio = sd_gDNA_ratio / sqrt(n()),
              is_hdr_allele = first(is_hdr_allele),
              is_wt_allele = first(is_wt_allele),
              has_crispr_deletion = first(has_crispr_deletion),
              deletion_start = first(deletion_start),
              deletion_end = first(deletion_end)) %>%
    dplyr::arrange(!is_wt_allele, !is_hdr_allele, desc(total_count))

  # Get UNS and associated stats per UDP
  udp.uns_stats.df = replicate.dels.df %>%
    dplyr::filter(!is_wt_allele) %>%
    group_by(udp) %>%
    do( udp_ratio_estimate(., replicates.df, numerator = "num_reads", denominator = "num_wt_reads") )

  udp.dels.df = udp.dels.df %>%
    dplyr::left_join(udp.uns_stats.df %>% rename(uns = effect, uns_se = effect_sd, uns_confint_lo = effect_confint_lo,
                                                 uns_confint_hi = effect_confint_hi, uns_df_est = df_estimated),
                     by="udp") %>%
    dplyr::arrange(!is_wt_allele, !is_hdr_allele, -total_count)

  return(udp.dels.df)
}


fit_variance_components = function(replicate.udp.spread.df, replicates.type.df) {
  replicate.udp.spread.df = replicate.udp.spread.df %>% as.data.frame()
  rownames(replicate.udp.spread.df) = replicate.udp.spread.df$udp

  region_name = replicates.type.df$name[1]
  region_type = replicates.type.df$type[1]
  # Get all columns of the replicate dataframe which begin with "replicate",
  # but ensure that they have more than one factor level
  replicate_cols = colnames(replicates.type.df)[grepl("^vp_", colnames(replicates.type.df))]
  remove_cols = c()
  for (col in replicate_cols) {
    replicates.type.df[, col] = as.factor(replicates.type.df[, col, drop=T])
    if (length(unique(replicates.type.df[, col, drop=T])) <= 1) {
      warning(sprintf("Column %s cannot be used in variance components analysis since it has only one level for region %s.", col, region_name))
      remove_cols = c(col, remove_cols)
    }
  }
  replicate_cols = setdiff(replicate_cols, remove_cols)
  if (length(replicate_cols) == 0) {
    warning(sprintf("Region %s: no replicate variables are available to use in variance components analysis.", region_name))
    return(NULL)
  }
  replicates.type.df = replicates.type.df %>% select(replicate, one_of(replicate_cols))
  colnames(replicates.type.df) = gsub("vp_", "", colnames(replicates.type.df))
  replicate_cols = gsub("vp_", "", replicate_cols)

  formulaStr = sprintf("~ (1|%s)", replicate_cols[1])
  for (col in replicate_cols[-1]) {
    formulaStr = sprintf("%s + (1|%s)", formulaStr, col)
  }

  # IMPORTANT: Make sure that the columns of the expression input to
  # fitExtractVarPartMoel are in the same order of samples as the rows
  # of the "data" input!
  replicate.udp.expr = replicate.udp.spread.df[, replicates.type.df$replicate]

  # Remove any UDPs which have zero variance (or variance components call will fail)
  udpVariance = apply(replicate.udp.expr, MARGIN = 1, var)
  if (any(udpVariance == 0)) {
    warning(sprintf("Region %s, %s: removing %d UDPs with zero variance - cannot fit variance components for these.\n", region_name, region_type, sum(udpVariance == 0)))
    replicate.udp.expr = replicate.udp.expr[udpVariance > 0,]
  }
  if (nrow(replicate.udp.expr) == 0) {
    warning(sprintf("Region %s, %s: no UDPs to fit variance for after filtering\n", region_name, region_type))
    return(NULL)
  }

  message(sprintf("Region %s, %s: fitting variance components with model:\n    UDP_fraction %s\n", region_name, region_type, formulaStr))
  vp <- variancePartition::fitExtractVarPartModel(exprObj = replicate.udp.expr,
                                                  formula = as.formula(formulaStr),
                                                  data = replicates.type.df,
                                                  useWeights = F)
  #View(data.frame(vp))
  vp
}


#' Performs a variance components estimate for each deletion allele based on
#' the replicate metadata provided.
#'
#' @param del_result The result from deletion_analysis
#' @param replicates A data.frame defining GenIE replicate metadata.
#' @param allele_min_reads The minimum number of reads that a deletion allele must have across all replicates to be included.
#' @param allele_min_fraction The minimum fraction of total reads that a deletion allele must have across all replications to be included.
#' @return Returns a list with tables vp_cDNA and vp_gDNA, which partition variance according to the metadata columns that begin with "replicate_" in the 'replicates' parameter.
#' @examples
#' del_results = deletion_analysis(regions, replicates)
#' vc = get_variance_components(del_results[[1]], replicates)
#' variance_components_plot(vc)
#' @seealso \code{\link{deletion_analysis}}
#' @seealso \code{\link{variance_components_plot}}
#' @export
#'
get_variance_components = function(del_result, replicates, allele_min_reads = 100, allele_min_fraction = 0.001) {
  check_is_del_result(del_result)
  # Check that the input replicate information includes all replicates
  # in the del_result
  if (!all(del_result$replicates$name %in% replicates$name)) {
    stop("Not all replicates in the del_result are described in the replicates input dataframe.")
  }
  replicate_cols = colnames(replicates)[grepl("^vp_", colnames(replicates))]
  if (length(replicate_cols) < 1) {
    stop("No columns beginning with 'vp_' found in the replicates input dataframe.")
  }
  # Check that each replicate column has more than 1 distinct value,
  # and fewer distinct values than the number of replicates
  is_valid_col = function(col) {
    n = length(unique(replicates[,col,drop=T]))
    return(n > 1 & n < nrow(replicates))
  }
  validCols = sapply(replicate_cols, is_valid_col)
  if (any(!validCols)) {
    warning(sprintf("Not including the following replicate column(s) in variance analysis: %s",
                    paste(replicate_cols[!validCols], collapse = ", ")))
    replicate_cols = replicate_cols[validCols]
  }
  # Ensure that replicate cols are treated as factors
  for (col in replicate_cols) {
    replicates[, col] = as.factor(replicates[, col, drop=T])
  }

  # Subset replicate info dataframe to those replicates in replicate_alleles
  replicate_alleles = del_result$replicate_alleles
  replicates = replicates %>% dplyr::filter(replicate %in% unique(replicate_alleles$replicate))

  filter_udps = function(udp.df, allele_min_reads, allele_min_fraction) {
    # Spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
    # fill = 0 for when a replicate has no value for a given UDP, then
    # gather again.
    udp.filled.df = udp.df %>%
      dplyr::select(udp, type, replicate, num_reads) %>%
      tidyr::spread(replicate, num_reads, fill = 0) %>%
      tidyr::gather(key="replicate", value="num_reads", -udp, -type)
    totalReads = sum(udp.filled.df$num_reads)
    udp.summary = udp.filled.df %>% group_by(udp, type) %>%
      dplyr::summarise(sum = sum(num_reads), frac = sum / totalReads)

    # Remove UDPs with too few total reads
    udp.filled.df %>% dplyr::inner_join(udp.summary %>% dplyr::filter(sum >= allele_min_reads, frac > allele_min_fraction), by=c("udp", "type"))
  }

  # This function returns a value based on the "num_reads" for each replicate,
  # but which is what we use to determine variance components. Raw number of
  # reads is not really suitable. The recommended value is "read_fraction",
  # where this represents the fraction of reads that a UDP represents for a
  # given replicate.
  get_value_for_variance = function(replicate_alleles, method) {
    replicate.numreads = replicate_alleles %>%
      dplyr::group_by(replicate) %>%
      dplyr::summarise(replicate_num_reads = sum(num_reads))
    replicate_alleles = replicate_alleles %>%
      dplyr::left_join(replicate.numreads, by="replicate") %>%
      dplyr::mutate(udp_fraction = num_reads / replicate_num_reads)
    return(replicate_alleles$udp_fraction)
  }

  # Function to get variance components result for cDNA and gDNA for a single region
  get_region_variance_components = function(replicate_alleles, replicates.type.df, dna_type, allele_min_reads, allele_min_fraction) {
    region = replicate_alleles$name[1]
    #replicates.type.df %>% dplyr::filter(name == region, type == dna_type)
    replicate.udp.type.df = filter_udps(replicate_alleles %>% dplyr::filter(type == dna_type),
                                       allele_min_reads = allele_min_reads, allele_min_fraction = allele_min_fraction)

    replicate.udp.type.df$value = get_value_for_variance(replicate.udp.type.df, method)

    # Spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
    # fill = 0 for when a replicate has no value for a given UDP
    replicate.udp.type.spread.df = replicate.udp.type.df %>%
      dplyr::select(udp, frac, replicate, value) %>%
      tidyr::spread(replicate, value, fill = 0)
    #udpFraction = rowMeans(replicate.udp.type.spread.df %>% dplyr::select(-udp))
    vp = fit_variance_components(replicate.udp.spread.df = replicate.udp.type.spread.df %>% filter(frac > 0.001) %>% select(-frac),
                                 replicates.type.df = replicates.type.df)
    vp.df = NULL
    if (!is.null(vp)) {
      vp.df = data.frame(vp) %>%
        mutate(udp = rownames(vp)) %>%
        left_join(replicate.udp.type.spread.df %>% select(udp, frac), by="udp") %>%
        mutate(name = region) %>%
        select(name, udp, frac, everything())
    }
    return(vp.df)
  }

  #########################################################################
  dna_type = "cDNA"
  regions = unique(replicate_alleles$name)
  dflist = list()
  for (i in 1:length(regions)) {
    dflist[[i]] = get_region_variance_components(replicate_alleles %>% filter(name == regions[i]),
                                                 replicates %>% dplyr::filter(name == regions[i], type == dna_type),
                                                 dna_type, allele_min_reads, allele_min_fraction)
  }
  vp_cDNA = as_tibble(bind_rows(dflist))

  # Do the same for gDNA
  dna_type = "gDNA"
  regions = unique(replicate_alleles$name)
  dflist = list()
  for (i in 1:length(regions)) {
    dflist[[i]] = get_region_variance_components(replicate_alleles %>% filter(name == regions[i]),
                                                 replicates %>% dplyr::filter(name == regions[i], type == dna_type),
                                                 dna_type, allele_min_reads, allele_min_fraction)
  }
  vp_gDNA = as_tibble(bind_rows(dflist))

  opts = list(del_result = del_result,
              replicates = replicates,
              allele_min_reads = allele_min_reads,
              allele_min_fraction = allele_min_fraction)
  return(list(name = del_result$region$name,
              vp_cDNA = vp_cDNA,
              vp_gDNA = vp_gDNA,
              opts = opts))
}


get_udp_stats = function(replicate_udps, dna_type, allele_min_reads) {
  # Spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
  # fill = 0 for when a replicate has no value for a given UDP, then
  # gather again
  replicate_udp_filled = replicate_udps %>%
    dplyr::filter(type == dna_type) %>%
    dplyr::select(udp, replicate, num_reads) %>%
    tidyr::spread(replicate, num_reads, fill = 0)

  replicate_udp_filled$total_count = apply(replicate_udp_filled[,-1], MARGIN = 1, FUN = sum)
  replicate_udp_filled = replicate_udp_filled %>%
    dplyr::filter(total_count > allele_min_reads) %>%
    dplyr::select(-total_count)
  replicate_udp_filled = replicate_udp_filled %>%
    tidyr::gather(key="replicate", value="num_reads", -udp)

  replicate_reads = replicate_udp_filled %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(replicate_reads = sum(num_reads))
  replicate_udp_filled = replicate_udp_filled %>%
    dplyr::left_join(replicate_reads, by="replicate")
  mean_replicate_reads = mean(replicate_reads$replicate_reads)
  stats.df = replicate_udp_filled %>% dplyr::group_by(udp) %>%
    dplyr::summarise(type = dna_type,
                     udp_mean = mean(num_reads),
                     udp_sd = sd(num_reads),
                     udp_cv = udp_sd / udp_mean,
                     udp_frac_mean = mean(num_reads) / mean_replicate_reads,
                     udp_frac_sd = sd(num_reads / replicate_reads),
                     udp_frac_cv = udp_frac_sd / udp_frac_mean)
  stats.df = stats.df %>%
    dplyr::left_join(replicate_udps %>% select(udp, is_hdr_allele, is_wt_allele) %>% unique(), by="udp")
  stats.df
}


#' Performs estimates of power to detect effects of alleles at different
#' read fractions, given the variance observed in the del_result replicates.
#'
#' @param del_result The result from deletion_analysis
#' @param allele_min_reads The minimum number of reads that a deletion allele must have across all replicates to be included.
#' @param WT_fraction If specified, then the model will assume this fraction of WT reads
#' @return Returns...
#' @examples
#' del_results = deletion_analysis(regions, replicates)
#' pwr = power_analysis(del_results[[1]], replicates)
#' power_plots(pwr)
#' @seealso \code{\link{deletion_analysis}}
#' @seealso \code{\link{power_plots}}
#' @export
#'
power_analysis = function(del_result, allele_min_reads = 100, WT_fraction = NA) {
  opts = list(allele_min_reads = allele_min_reads, WT_fraction = WT_fraction)
  check_is_del_result(del_result)
  replicate_alleles = del_result$replicate_alleles
  # We determine power based on %editing, i.e. the % of all reads represented
  # by a given edit. The %editing can be considered equivalent to mean UDP read
  # count across replicates.
  # To determine power, we need the standard deviation expected across cDNA
  # replicates for a given %editing. We can get this by fitting a curve to
  # the empirical curve of coefficient of variation (CV) vs. UDP read count.

  # An error may occur in NLS, so we need to catch it
  res = tryCatch({
    # First fit a curve to a table of CV vs. mean read count.
    ## TODO: consider changing the modeling form to a free exponent, e.g.
    ## udp_frac_cv ~ a + b * udp_frac_mean^c
    ## along with constraints such as a > 0, b > 0, -1 <= c <= 0

    coefs = list()
    # Restrict to UDPs > 0.1% fraction, since noise at very low fractions
    # has a big effect on the fit
    cDNA_stats = get_udp_stats(replicate_alleles, "cDNA", allele_min_reads) %>%
      filter(udp_frac_mean > 0.001)
    fit.cDNA <- nls(udp_frac_cv ~ a + b * (1/sqrt(udp_frac_mean)), data = cDNA_stats,
                    weights = sqrt(cDNA_stats$udp_mean), start = list(a = 0.1, b = 0.1))
    coefs$cDNA <- coefficients(fit.cDNA)

    gDNA_stats = get_udp_stats(replicate_alleles, "gDNA", allele_min_reads) %>%
      filter(udp_frac_mean > 0.001)
    fit.gDNA <- nls(udp_frac_cv ~ a + b * (1/sqrt(udp_frac_mean)), data = gDNA_stats,
                    weights = sqrt(gDNA_stats$udp_mean), start = list(a = 0.1, b = 0.1))
    coefs$gDNA <- coefficients(fit.gDNA)
  }
  , error = function(e) e, warning = function(w) w)
  if (is(res, "error")) {
    stop("Error in NLS fitting to UDP CV vs. UDP fraction.\nPower results cannot be calculated")
  }

  udp_stats = rbind(cDNA_stats, gDNA_stats) %>%
    dplyr::mutate(allele_type = if_else(is_hdr_allele, "HDR", if_else(is_wt_allele, "WT", "Del")))

  power_res = list(opts = opts,
                   coefs = coefs,
                   udp_stats = udp_stats)

  n_cDNA_rep = length(del_result$replicates %>% dplyr::filter(type == "cDNA") %>% .$replicate)
  n_gDNA_rep = length(del_result$replicates %>% dplyr::filter(type == "gDNA") %>% .$replicate)

  # The general formula for propagation of uncertainty is (where f = A/B):
  #     sd(f) = f * sqrt( (sd(A) / A)^2 + (sd(B) / B)^2 - 2*cov(A,B) )
  # In some cases we may be able to assume that cov(A,B) is zero, e.g.
  # where A and B are independent replicates. Ignoring the covariance
  # can only lead to a larger estimate of the uncertainty, so is conservative.

  # Our effect size is the ratio:  (UDP_frac_cDNA / WT_frac_cDNA) / (UDP_frac_gDNA / WT_frac_gDNA)
  # This is a lot of parameters to vary, so we also allow specifying a
  # WT fraction to use, rather than using the actual WT CV.
  if (!is.na(WT_fraction)) {
    if (WT_fraction <= 0 || WT_fraction > 1) {
      stop(sprintf("Invalid WT_fraction (%f). WT_fraction should be greater than 0 and less than 1.", WT_fraction))
    }
    power_res$cDNA_WT_CV = fit_UDP_CV(coefs$cDNA, WT_fraction)
    power_res$gDNA_WT_CV = fit_UDP_CV(coefs$gDNA, WT_fraction)
  } else {
    power_res$cDNA_WT_CV = udp_stats %>% filter(type == "cDNA", is_wt_allele) %>% .$udp_frac_cv
    power_res$gDNA_WT_CV = udp_stats %>% filter(type == "gDNA", is_wt_allele) %>% .$udp_frac_cv
  }

  # Create a table with combinations of UDP fractions and effect sizes for calculating power
  power.combinations = expand.grid(udp_fraction = c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5),
                         effect_size = c(1.05, 1.1, 1.2, 1.5, 2.0))
  # Add to this a combination of either the true numcDNA/numgDNA for each locus, or assume
  # a fixed value of 8/4.
  power.df = power.combinations %>% mutate(num_cDNA_rep = n_cDNA_rep, num_gDNA_rep = n_gDNA_rep)
  if (!(n_cDNA_rep == 8 & n_gDNA_rep == 4)) {
    # This check means we don't duplicate the same power combinations for 8/4 reps. But more
    # importantly, the spread call below fails if there are duplicate combinations!
    power.df = bind_rows(power.df, power.combinations %>% mutate(num_cDNA_rep = 8, num_gDNA_rep = 4))
  }
  power_res$power = power.df %>%
    mutate(power = calc_power(power_res, effect_size, udp_fraction, udp_fraction, num_cDNA_rep, num_gDNA_rep),
           udp_fraction = paste("udp frac", udp_fraction)) %>%
    tidyr::spread(key="udp_fraction", value="power")

  # Determine power as the allocation of replicates to cDNA vs. gDNA changes.
  get_power_df = function(power_res, nRep, udpFraction) {
    df1.1 = tibble(num_cDNA = seq(1:(nRep-1)), effectSize = 1.1, udpFraction = udpFraction)
    df1.1$power = calc_power(power_res, df1.1$effectSize, df1.1$udpFraction, df1.1$udpFraction, df1.1$num_cDNA, nRep - df1.1$num_cDNA)
    df1.2 = tibble(num_cDNA = seq(1:(nRep-1)), effectSize = 1.2, udpFraction = udpFraction)
    df1.2$power = calc_power(power_res, df1.2$effectSize, df1.2$udpFraction, df1.2$udpFraction, df1.2$num_cDNA, nRep - df1.2$num_cDNA)
    df = rbind(df1.1, df1.2) %>%
      mutate(nReps = nRep, cDNA_rep_fraction = num_cDNA / nRep) %>%
      select(nReps, num_cDNA, cDNA_rep_fraction, everything())
    #df$effectSize = sprintf("%.1fx", df$effectSize)
    df
  }

  power_res$replicate_allocation = expand.grid(udp_fraction = c(0.01, 0.1), nReps = c(6, 10, 15, 25)) %>%
    group_by(nReps, udp_fraction) %>%
    group_map(~ get_power_df(power_res, .y$nReps, .y$udp_fraction) %>% select(-nReps, -udpFraction))

  return(power_res)
}

fit_UDP_CV = function(coefs, udp_frac_mean) {
  coefs["a"] + coefs["b"]/sqrt(udp_frac_mean)
}

calc_power = function(power_res, effect, cDNA_frac, gDNA_frac, n_cDNA, n_gDNA) {
  # We use the propagation of uncertainty formula to get the SE of UDP_frac / WT_frac.
  se_cDNA_WT_ratio = function(power_res, frac, n_cDNA) {
    frac * sqrt((fit_UDP_CV(power_res$coefs$cDNA, frac)/sqrt(n_cDNA))^2 + (power_res$cDNA_WT_CV/sqrt(n_cDNA))^2)
  }
  se_gDNA_WT_ratio = function(power_res, frac, n_gDNA) {
    frac * sqrt((fit_UDP_CV(power_res$coefs$gDNA, frac)/sqrt(n_gDNA))^2 + (power_res$gDNA_WT_CV/sqrt(n_gDNA))^2)
  }
  se_ratio_cDNA_gDNA = function(power_res, cDNA_frac, gDNA_frac, n_cDNA, n_gDNA) {
    (cDNA_frac / gDNA_frac) * sqrt( (se_cDNA_WT_ratio(power_res, cDNA_frac, n_cDNA) / cDNA_frac)^2 +
                                    (se_gDNA_WT_ratio(power_res, gDNA_frac, n_gDNA) / gDNA_frac)^2 )
  }

  # Note that sd(A)/A is the coefficient of variation, which is why we use CV below.
  # Here we are computing the SE of the cDNA / gDNA ratio estimate, where this
  # is based on all replicates done for each. Hence, we use the standard error
  # for each of these estimates (cDNA:WT ratio and gDNA:WT ratio). The SE of
  # these estimates is the SD / sqrt(n).
  se_ratio_cDNA_gDNA_old = function(cDNA_frac, gDNA_frac, n_cDNA, n_gDNA) {
    (cDNA_frac / gDNA_frac) * sqrt( (fitcDNAUDP_CV(cDNA_frac) / sqrt(n_cDNA))^2 +
                                      (fitgDNAUDP_CV(gDNA_frac) / sqrt(n_gDNA))^2 )
  }

  # This seems wrong - I don't know why I did this previously. In any case,
  # it only leads to a slightly larger estimated standard error of the ratio.
  # sd_cDNA_udp_frac = function(frac) { frac * sqrt(fitcDNAUDP_CV(frac)^2 + fitcDNAUDP_CV(1-frac)^2) }
  # sd_gDNA_udp_frac = function(frac) { frac * sqrt(fitgDNAUDP_CV(frac)^2 + fitgDNAUDP_CV(1-frac)^2) }
  # sd_ratio_cDNA_gDNA = function(cDNA_frac, gDNA_frac, n_cDNA, n_gDNA) {
  #   (cDNA_frac / gDNA_frac) * sqrt( (sd_cDNA_udp_frac(cDNA_frac) / cDNA_frac / sqrt(n_cDNA))^2 +
  #                                     (sd_gDNA_udp_frac(gDNA_frac) / gDNA_frac / sqrt(n_gDNA))^2 )
  # }

  effect_se = se_ratio_cDNA_gDNA(power_res, cDNA_frac, gDNA_frac, n_cDNA, n_gDNA)
  Zscore = (effect - 1) / effect_se
  pvalue = 2 * pnorm(-abs(Zscore))
  1 - pvalue
}

get_read_udp = function(seq, profile_chars = character(0)) {
  seq_chars = strsplit(seq, "")[[1]]
  get_read_chars_udp(strsplit(seq, "")[[1]], profile_chars)
}

get_read_chars_udp = function(seq_chars, profile_chars = character(0)) {
  if (length(profile_chars) > 0) {
    if (length(profile_chars) != length(seq_chars)) {
      stop(sprintf("Length of sequence to get UDP (%d chars) should be the same as the HDR profile length (%d chars).", length(seq_chars), length(profile_chars)))
    }
    seq_chars[seq_chars != '*' & profile_chars == "-"] = "-"
  } else {
    seq_chars[seq_chars != '*'] = "-"
  }
  return(str_c(seq_chars, collapse = ""))
}

# When there are "ambiguous" DNA letters in a sequence, we can make a set of
# sequences that represent all possible values. E.g. H indicates A or C or T (but not G).
# If two positions have ambiguous letters, then the combinations multiply.
# To handle this case, I use recursion below, without checking whether the combinations
# may become too large.
get_hdr_profiles = function(hdr_profile) {
  hdr_profile_set = c()
  if (grepl("B|D|H|V", hdr_profile)) {
    hits = str_locate_all(hdr_profile, "B|D|H|V")[[1]]
    for (i in 1:nrow(hits)) {
      nt = str_sub(hdr_profile, hits[i,"start"], hits[i,"start"])
      hdr1 <- hdr2 <- hdr3 <- hdr_profile
      if (nt == "B") {
        str_sub(hdr1, hits[i,"start"], hits[i,"start"]) = "C"
        str_sub(hdr2, hits[i,"start"], hits[i,"start"]) = "G"
        str_sub(hdr3, hits[i,"start"], hits[i,"start"]) = "T"
      } else if (nt == "D") {
        str_sub(hdr1, hits[i,"start"], hits[i,"start"]) = "A"
        str_sub(hdr2, hits[i,"start"], hits[i,"start"]) = "G"
        str_sub(hdr3, hits[i,"start"], hits[i,"start"]) = "T"
      } else if (nt == "H") {
        str_sub(hdr1, hits[i,"start"], hits[i,"start"]) = "A"
        str_sub(hdr2, hits[i,"start"], hits[i,"start"]) = "C"
        str_sub(hdr3, hits[i,"start"], hits[i,"start"]) = "T"
      } else if (nt == "V") {
        str_sub(hdr1, hits[i,"start"], hits[i,"start"]) = "A"
        str_sub(hdr2, hits[i,"start"], hits[i,"start"]) = "C"
        str_sub(hdr3, hits[i,"start"], hits[i,"start"]) = "G"
      }
      hdr_profile_set = c(get_hdr_profiles(hdr1), get_hdr_profiles(hdr2), get_hdr_profiles(hdr3))
    }
  } else {
    hdr_profile_set = hdr_profile
  }
  hdr_profile_set
}

# This method was slower unless there are many non-empty "profile chars"
# get_read_site_profile_OLD = function(seq, profile_chars) {
#   seq_chars = strsplit(seq, "")[[1]]
#   if (length(profile_chars) != length(seq_chars)) {
#     stop(sprintf("Length of sequence to get read profile (%d chars) should be the same as the profile length (%d chars).", length(seq_chars), length(profile_chars)))
#   }
#   out_chars = seq_chars[is_dna_letter(profile_chars)]
#   return(paste0(out_chars, collapse = ""))
# }

# profile_positions should have the indices of DNA letters to get in the
# input sequence
get_read_site_profile = function(seq, profile_positions) {
  out_str = strrep(" ", times = length(profile_positions))
  for (i in 1:length(profile_positions)) {
    substr(out_str, i, i) = substr(seq, profile_positions[i], profile_positions[i])
  }
  return(out_str)
}

# get_read_mismatch_profile = function(seq, ref_sequence) {
#   #print(sprintf("Seq:%s\t%s", seq, ref_sequence))
#   seq_chars = strsplit(seq, "")[[1]]
#   ref_chars = strsplit(ref_sequence, "")[[1]]
#   get_read_chars_mismatch_profile(seq_chars, ref_chars)
# }
get_read_chars_mismatch_profile = function(seq_chars, ref_chars) {
  output = rep("-", length(ref_chars))
  mismatchpos = (seq_chars != ref_chars)
  output[mismatchpos] = seq_chars[mismatchpos]
  return(str_c(output, collapse=""))
}

str_to_chars = function(s) {
  strsplit(s,"")[[1]]
}

is_dna_letter = function(s_chars) {
  return(s_chars %in% c('A', 'C', 'G', 'T'))
}

get_mismatch_count = function(s, ref) {
  s_chars = str_to_chars(s)
  ref_chars = str_to_chars(ref)
  letterPos = is_dna_letter(s_chars)
  sum(s_chars[letterPos] != ref_chars[letterPos])
}

get_mismatch_chars_count = function(s_chars, ref_chars) {
  letterPos = is_dna_letter(s_chars)
  sum(s_chars[letterPos] != ref_chars[letterPos])
}

udp_ratio_estimate = function(udp.counts.df, replicates.df, numerator = "num_reads", denominator = "num_wt_reads") {
  num_reads = udp.counts.df[, numerator, drop = T]
  num_wt_reads = udp.counts.df[, denominator, drop = T]
  udp.counts.df$num_reads = num_reads
  udp.counts.df$num_wt_reads = num_wt_reads

  if (sum(udp.counts.df$num_reads) == 0) {
    # We can't run a model when all counts are zeroes
    return(tibble(effect = 0, effect_sd = NA, effect_confint_lo = NA, effect_confint_hi = NA, df_estimated = NA, pval = NA, method = NA))
  }

  numWtZero = sum(udp.counts.df$num_wt_reads == 0)
  if (numWtZero > 0) {
    udp.counts.df = udp.counts.df %>% filter(num_wt_reads != 0)
    warning(sprintf("In udp_ratio_estimate: removing %d replicates with zero WT reads.", numWtZero))
  }
  udp_ratio_estimate_welch(cDNA_UDP_counts = udp.counts.df$num_reads[udp.counts.df$type == "cDNA"],
                           cDNA_wt_counts = udp.counts.df$num_wt_reads[udp.counts.df$type == "cDNA"],
                           gDNA_UDP_counts = udp.counts.df$num_reads[udp.counts.df$type == "gDNA"],
                           gDNA_wt_counts = udp.counts.df$num_wt_reads[udp.counts.df$type == "gDNA"])
}

udp_ratio_estimate_welch = function(cDNA_UDP_counts, cDNA_wt_counts, gDNA_UDP_counts, gDNA_wt_counts) {
  N_cDNA = length(cDNA_wt_counts)
  N_gDNA = length(gDNA_wt_counts)

  # There are two ways that we could calculate the UDP:WT ratio (for both cDNA and gDNA).
  # One is to compute it for each replicate and take the average. The other is to sum
  # UDP counts for replicate, and divide by the summed WT counts for replicate. This latter
  # method effectively weights the ratio by the number of reads in each replicate. The main
  # problem with this method is that it gives a value that will differ from the one we
  # are assuming in our statistical test further down (t test). For example, our confidence
  # interval may include 1 (no effect), but our p value gives a significant result. For this
  # reason I use the unweighted UDP:WT ratios. It has the downside of not using all of the
  # reads across the replicates "efficiently", but if outliers are excluded properly, then
  # it should be a more robust measure of the true effect.
  cDNA_ratio = mean(cDNA_UDP_counts / cDNA_wt_counts)
  gDNA_ratio = mean(gDNA_UDP_counts / gDNA_wt_counts)
  sd_cDNA_ratio = sd(cDNA_UDP_counts / cDNA_wt_counts)
  sd_gDNA_ratio = sd(gDNA_UDP_counts / gDNA_wt_counts)
  SE_cDNA_ratio = sd_cDNA_ratio / sqrt(N_cDNA)
  SE_gDNA_ratio = sd_gDNA_ratio / sqrt(N_gDNA)

  effect = cDNA_ratio / gDNA_ratio
  effect_se = ratio_uncertainty_sd(cDNA_ratio, SE_cDNA_ratio, gDNA_ratio, SE_gDNA_ratio, 0)

  # Use a Welch two-sample T test to see whether the UDP:WT ratio differs
  # in cDNA vs. gDNA. This test accounts for there being different variance
  # in cDNA and gDNA.
  pvalue = NaN
  tryres = tryCatch({
    res = t.test(cDNA_UDP_counts / cDNA_wt_counts,
                 gDNA_UDP_counts / gDNA_wt_counts,
                 alternative = "two.sided", var.equal = F)
    pvalue = res$p.value
  }, error = function(e) e, warning = function(w) w)
  if (is(tryres, "error")) {
    warning("Error in t.test")
  }
  # Use a t distribution to determine the confidence interval. We estimate
  # the degrees of freedom using the Welch T test approximation. This gives
  # a wider confidence interval than if we assumed a normal distribution or
  # a standard T distribution.
  df_estimated = welch_df_estimate(sd_cDNA_ratio, N_cDNA, sd_gDNA_ratio, N_gDNA)
  tScore_alpha_0.05 = qt(0.975, df = df_estimated, lower.tail=T)
  effect_confint = tScore_alpha_0.05 * effect_se

  return(tibble(effect = effect,
                effect_sd = effect_se,
                effect_confint_lo = effect - tScore_alpha_0.05 * effect_se,
                effect_confint_hi = effect + tScore_alpha_0.05 * effect_se,
                df_estimated = df_estimated,
                pval = pvalue))
}

welch_df_estimate = function(sdA, nA, sdB, nB) {
  round((sdA^2/nA + sdB^2/nB)^2 / ( (sdA^2/nA)^2 / (nA - 1) + (sdB^2/nB)^2 / (nB - 1) ))
}

ratio_uncertainty_sd = function(meanA, sdA, meanB, sdB, covAB) {
  abs(meanA / meanB) * sqrt((sdA / meanA)^2 + (sdB / meanB)^2 - (2*covAB / (meanA*meanB)))
}

print_list = function(l, prefix = "    ") {
  list.df = tibble(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}

is_empty = function(s) {
  ifelse(is.null(s), T, s == "")
}

check_is_del_result = function(del_result) {
  if (!(hasName(del_result, "type") && del_result$type == "deletion_analysis")) {
    stop("Input del_result object is not of type deletion_analysis.")
  }
}

check_is_grep_result = function(grep_result) {
  if (!(hasName(grep_result, "type") && grep_result$type == "grep_analysis")) {
    stop("Input grep_result object is not of type grep_analysis.")
  }
}

output = function(opts, str) {
  if (!opts$quiet) {
    cat(str)
  }
}


#' Given results for rgenie grep or deletion analysis of a set of regions,
#' merges together tables across regions.
#'
#' @param results A list of rgenie results from either grep or deletion analyses.
#' @return Returns a list containing the same tables as in an individual result,
#' but concatenated across regions.
#' @examples
#' del_results = deletion_analysis(regions, replicates)
#' del_tables = bind_results(del_results)
#' @seealso \code{\link{grep_analysis}}
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
bind_results = function(results) {
  if (is.null(results)) {
    return(NULL)
  }
  if (length(results) < 1) {
    return(results)
  }
  if (results[[1]]$type == "grep_analysis") {
    new_res = list()
    new_res$replicate_stats = bind_rows(lapply(results, function(res) res$replicate_stats))
    new_res$region_stats = bind_rows(lapply(results, function(res) res$region_stats))
    return(new_res)
  }
  if (results[[1]]$type == "deletion_analysis") {
    new_res = list()
    new_res$replicate_qc = bind_rows(lapply(results, function(res) res$replicate_qc))
    new_res$replicate_alleles = bind_rows(lapply(results, function(res) res$replicate_alleles))
    new_res$region_alleles = bind_rows(lapply(results, function(res) res$region_alleles))
    new_res$allele_effect = bind_rows(lapply(results, function(res) res$allele_effect))
    new_res$site_profiles = bind_rows(lapply(results, function(res) res$site_profiles))
    new_res$mismatch_profiles = bind_rows(lapply(results, function(res) res$mismatch_profiles))
    new_res$replicate_stats = bind_rows(lapply(results, function(res) res$replicate_stats))
    new_res$region_stats = bind_rows(lapply(results, function(res) res$region_stats))
    return(new_res)
  }
  stop("Unknown results type.")
}


#' Downloads example data for rgenie.
#'
#' The example data is a set of BAM files for GenIE replicates.
#'
#' @param dir Directory where example data should be put.
#' @return Returns a list containing the same tables as in an individual result,
#' but concatenated across regions.
#' @examples
#' download_example(dir = "~/genie_example", name = "MUL1")
#' # Data are downloaded and we can run an rgenie analysis
#' setwd("~/genie_example")
#' regions = readr::read_tsv("mul1.genie_regions.tsv")
#' replicates = readr::read_tsv("mul1.genie_replicates.tsv")
#' grep_results = grep_analysis(regions, replicates)
#' @seealso \code{\link{grep_analysis}}
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
download_example = function(dir = "~/genie_example", name = "MUL1", overwrite = F, quiet = F) {
  valid_examples = c("MUL1", "ABHD4")
  if (!(name %in% valid_examples)) {
    stop(sprintf("Name %s is not one of the available examples: {%s}", name, paste(valid_examples, collapse = ", ")))
  }
  # If directory doesn't exist, then create it
  ex_dir = file.path(dir, name)
  if (dir.exists(ex_dir) && !overwrite) {
    message(sprintf("Example data for %s is already present. To overwrite it, set overwrite = TRUE.", name))
  } else {
    dir.create(ex_dir, recursive = TRUE, showWarnings = FALSE)
    fpath = sprintf("https://raw.githubusercontent.com/Jeremy37/rgenie_example/master/file_list.%s.txt", name)
    file_list = readr::read_csv(url(fpath), col_names = "path")$path

    for (fpath in file_list) {
      path_parts = strsplit(fpath, "/", fixed = T)[[1]]
      fname = path_parts[length(path_parts)]
      destfile = file.path(ex_dir, fname)
      download.file(fpath, destfile, quiet = quiet)
    }
    message(sprintf("Downloaded data for example '%s' to %s", name, ex_dir))
  }
}

