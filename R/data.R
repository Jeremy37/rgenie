#' GenIE regions for the MUL1 example.
#'
#' A data.frame with GenIE regions details for the MUL1 example.
#'
#' @format A data frame with 1 row and 9 variables:
#' \describe{
#'   \item{name}{region name}
#'   \item{sequence_name}{chromosome or amplicon name}
#'   \item{start}{start coordinate of the region}
#'   \item{end}{end coordinate of the region}
#'   \item{highlight_site}{coordinate of the SNP site}
#'   \item{cut_site}{coordinate of the cut site}
#'   \item{hdr_allele_profile}{HDR allele profile}
#'   \item{wt_allele_profile}{WT allele profile}
#'   \item{ref_sequence}{sequence of the amplicon region}
#' }
"mul1_regions"


#' GenIE replicates for the MUL1 example.
#'
#' A data.frame with GenIE replicate details for the MUL1 example.
#'
#' @format A data frame with 12 rows and 5 variables:
#' \describe{
#'   \item{name}{region name}
#'   \item{replicate}{short name of the replicate}
#'   \item{type}{cDNA" or "gDNA"}
#'   \item{vp_extraction}{which cDNA/gDNA extraction the replicate is from}
#'   \item{bam}{path to the BAM file of data for the replicate}
#' }
"mul1_replicates"


#' GenIE grep results list.
#'
#' A list with a single result object from GenIE grep analysis for the MUL1 example.
#'
#' @format A one-item list with a grep_analysis result object.
"mul1_grep_results"


#' GenIE alignment analysis results for the MUL1 example.
#'
#' A list with a single result object from GenIE alignment analysis for the MUL1 example.
#'
#' @format A one-item list with a alignment_analysis result object.
"mul1_del_results"


#' GenIE regions for the MUL1 example.
#'
#' A data.frame with GenIE regions details for the MUL1 example.
#'
#' @format A data frame with 1 row and 9 variables:
#' \describe{
#'   \item{name}{region name}
#'   \item{sequence_name}{chromosome or amplicon name}
#'   \item{start}{start coordinate of the region}
#'   \item{end}{end coordinate of the region}
#'   \item{highlight_site}{coordinate of the SNP site}
#'   \item{cut_site}{coordinate of the cut site}
#'   \item{hdr_allele_profile}{HDR allele profile}
#'   \item{wt_allele_profile}{WT allele profile}
#'   \item{ref_sequence}{sequence of the amplicon region}
#' }
"atac_regions"


#' GenIE replicates for the MUL1 example.
#'
#' A data.frame with GenIE replicate details for the MUL1 example.
#'
#' @format A data frame with 12 rows and 5 variables:
#' \describe{
#'   \item{name}{region name}
#'   \item{replicate}{short name of the replicate}
#'   \item{type}{cDNA" or "gDNA"}
#'   \item{vp_extraction}{which cDNA/gDNA extraction the replicate is from}
#'   \item{bam}{path to the BAM file of data for the replicate}
#' }
"atac_replicates"


#' GenIE ATAC example rs7729529 grep results.
#'
#' A list with a single result object from GenIE grep analysis for the ATAC_rs7729529 example.
#'
#' @format A one-item list with a grep_analysis result object.
"atac_grep_results"


#' GenIE ATAC example rs7729529 alignment analysis results.
#'
#' A list with a single result object from GenIE alignment analysis for the ATAC_rs7729529 example.
#'
#' @format A one-item list with an alignment_analysis result object.
"atac_del_results"

