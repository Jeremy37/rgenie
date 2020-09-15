

#' Returns a list with default options for all rgenie plots, useful in calls to deletion_plots().
#'
#' @examples
#' \dontrun{
#' del_result = deletion_analysis(regions, replicates)
#' opts = genie_plot_options()
#' plot_list = deletion_plots(del_result, opts)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
genie_plot_options = function() {
  default_opts = list(
    viewing_window = 40,
    outlier_threshold = NA,
    allele_plot_max_alleles = 40,
    variance_analysis_min_count = 100,
    variance_analysis_min_fraction = 0.001,
    variance_plot_split_by_fraction = F,
    deletion_alleles_plot_color_by = "window",
    plower_plot_min_count = 100,
    plower_plot_WT_fraction = NA)
  return(default_opts)
}


#' Plots a summary of genie results across multiple regions.
#'
#' @param grep_results Result from a call to grep_analysis (or NULL).
#' @param del_results Result from a call to deletion_analysis (or NULL).
#' @return Returns a ggplot object.
#'
#' @examples
#' \dontrun{
#' grep_results = grep_analysis(regions, replicates)
#' del_results = deletion_analysis(regions, replicates)
#' experiment_summary_plot(grep_results, del_results)
#' }
#' @seealso \code{\link{grep_analysis}}
#' @seealso \code{\link{deletion_analysis}}
#' @seealso \code{\link{bind_results}}
#' @export
#'
experiment_summary_plot = function(grep_results, del_results) {
  getSignificanceStr = function(pval) {
    if (is.na(pval) | pval >= 0.01) { "p >= 0.01" }
    else if (pval < 0.001) { "p < 0.001" }
    else { "p < 0.01" }
  }

  p.grep.effect = NULL
  p.stats.effect = NULL
  effect_size_theme = theme_bw(10) + theme(axis.text.x = element_blank(),
                                           legend.title = element_blank(),
                                           axis.title.x = element_blank(),
                                           plot.margin = unit(c(0.1, 0, 0.1 ,1), "cm"))
  if (!is.null(grep_results)) {
    grep_dfs = bind_results(grep_results)
    exp_names = sapply(grep_dfs$region_stats$name, FUN = function(s) strsplit(s, ",", T)[[1]][1])
    if (any(duplicated(exp_names))) {
      exp_names = grep_dfs$region_stats$name
    }
    grep_dfs$region_stats$name = factor(as.character(exp_names), levels=exp_names)
    grep_dfs$region_stats$hdr_significance = factor(sapply(grep_dfs$region_stats$pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))

    p.grep.effect = ggplot(grep_dfs$region_stats, aes(x=name, y=effect, fill=hdr_significance)) +
      geom_bar(stat = "identity", width=0.5) +
      geom_errorbar(aes(ymin = effect_confint_lo, ymax = effect_confint_hi),
                    width = 0.2, col = "grey30") +
      geom_hline(yintercept = 1, col = "red") +
      effect_size_theme +
      ylab("HDR effect") + ggtitle("HDR effect size - grep analysis") +
      scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3")) +
      coord_cartesian(ylim = c(0, max(1, max(grep_dfs$region_stats$effect * 1.05, na.rm = T))))

    hdr.df = grep_dfs$region_stats
  }

  if (!is.null(del_results)) {
    del_dfs = bind_results(del_results)
    exp_names = sapply(del_dfs$region_stats$name, FUN = function(s) strsplit(s, ",", T)[[1]][1])
    if (any(duplicated(exp_names))) {
      exp_names = del_dfs$region_stats$name
    }
    del_dfs$region_stats$name = factor(as.character(exp_names), levels=exp_names)
    del_dfs$region_stats$hdr_significance = factor(sapply(del_dfs$region_stats$hdr_pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))

    p.stats.effect = ggplot(del_dfs$region_stats, aes(x=name, y=hdr_effect, fill=hdr_significance)) +
      geom_bar(stat = "identity", width=0.5) +
      geom_errorbar(aes(ymin = hdr_effect_confint_lo, ymax = hdr_effect_confint_hi),
                    width = 0.2, col = "grey30") +
      geom_hline(yintercept = 1, col = "red") +
      effect_size_theme +
      ylab("HDR effect") + ggtitle("HDR effect size - alignment analysis") +
      scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3")) +
      coord_cartesian(ylim = c(0, max(1, max(del_dfs$region_stats$hdr_effect * 1.05, na.rm = T))))

    del_dfs$region_stats$del_significance = factor(sapply(del_dfs$region_stats$del_pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))
    p.stats.del = ggplot(del_dfs$region_stats, aes(x=name, y=del_effect, fill=del_significance)) +
      geom_bar(stat = "identity", width=0.5) +
      geom_errorbar(aes(ymin = del_effect_confint_lo, ymax = del_effect_confint_hi),
                    width = 0.2, col = "grey30") +
      geom_hline(yintercept = 1, col = "red") +
      effect_size_theme +
      ylab("Del effect") + ggtitle("Deletion effect size") +
      scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3")) +
      coord_cartesian(ylim = c(0, min(3, max(del_dfs$region_stats$del_effect * 1.2, na.rm = T))))

    plot.df = del_dfs$region_stats %>% dplyr::select(name, HDR=hdr_rate_gDNA, NHEJ=del_rate_gDNA) %>%
      tidyr::gather(key = "type", value = "value", -name)
    plot.df$type = factor(as.character(plot.df$type), levels = c("NHEJ", "HDR"))
    p.stats.editing = ggplot(plot.df, aes(x=name, y=value*100, fill=type)) +
      geom_bar(stat = "identity", position = position_stack(), width=0.5) +
      theme_bw(10) + theme(axis.text.x = element_blank(),
                           legend.title = element_blank(),
                           axis.title.x = element_blank(),
                           plot.margin = unit(c(0.1, 0, 0.1 ,1), "cm")) +
      ylab("% editing") + ggtitle("Editing rates") +
      scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue"))

    hdr.df = del_dfs$region_stats
  }

  hdr.df$type = "HDR"
  p.stats.hdr = ggplot(hdr.df, aes(x=name, y=hdr_rate_gDNA * 100, fill=type)) +
    geom_bar(stat = "identity", position = position_stack(), width=0.5) +
    theme_bw(10) + theme(axis.text.x = element_text(angle = 37, hjust = 1, size=7),
                         legend.title = element_blank(),
                         axis.title.x = element_blank(),
                         plot.margin = unit(c(0.1, 0, 0.1 ,1), "cm")) +
    ylab("% HDR") + ggtitle("HDR rates") +
    scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue"))

  p.title = cowplot::ggdraw() + cowplot::draw_label("Experiment summary", fontface='bold')
  if (!is.null(p.grep.effect) & !is.null(p.stats.effect)) {
    p.res = egg::ggarrange(p.title, p.grep.effect, p.stats.effect, p.stats.del, p.stats.editing, p.stats.hdr, ncol=1, heights=c(1.2,2.4,2.4,2.4,2.4,2.4), draw = F)
  } else if (!is.null(p.grep.effect)) {
    p.res = egg::ggarrange(p.title, p.grep.effect, p.stats.hdr, ncol=1, heights=c(1.2,3,3), draw = F)
  } else {
    p.res = egg::ggarrange(p.title, p.stats.effect, p.stats.del, p.stats.editing, p.stats.hdr, ncol=1, heights=c(1.2,3,3,3,3), draw = F)
  }
  p.res
}


#' Plots a summary of genie results from a grep analysis.
#' @param grep_result Result from a call to grep_analysis.
#' @return Returns a ggplot object.
#'
#' @examples
#' \dontrun{
#' grep_results = grep_analysis(regions, replicates)
#' grep_summary_plot(grep_results[[1]])
#' }
#' @seealso \code{\link{grep_analysis}}
#' @export
#'
grep_summary_plot = function(grep_result) {
  check_is_grep_result(grep_result)
  counts.gDNA = grep_result$replicate_stats %>% filter(type == "gDNA")
  counts.cDNA = grep_result$replicate_stats %>% filter(type == "cDNA")

  if (is.null(grep_result$region_stats) && (nrow(counts.gDNA) < 2 || nrow(counts.cDNA) < 2)) {
    summary.left = sprintf("Unable to calculate stats with < 2 replicates")
    summary.right = ""
  } else {
    confIntervalString = function(ratioRes) {
      sprintf("95%% CI: (%.3g, %.3g)", ratioRes$effect_confint_lo, ratioRes$effect_confint_hi)
    }
    hdr.conf.interval.str = confIntervalString(grep_result$region_stats)
    hdr.summary = sprintf("cDNA:gDNA ratio (HDR/WT): %.3g\n%s,    p = %.3g",
                          grep_result$region_stats$effect, hdr.conf.interval.str, grep_result$region_stats$pval)

    hdr.frac.str = sprintf("Mean HDR frac gDNA: %.2g%%,  cDNA: %.2g%%", mean(counts.gDNA$HDR_frac) * 100, mean(counts.cDNA$HDR_frac) * 100)
    wt.frac.str = sprintf("Mean WT frac gDNA: %.2g%%,  cDNA: %.2g%%", mean(counts.gDNA$WT_frac) * 100, mean(counts.cDNA$WT_frac) * 100)
    summary.left = paste(hdr.frac.str, wt.frac.str, sep = "\n")
    summary.right = hdr.summary
  }

  # Convert to strings for nice printing (with 3 significant digits)
  stats.plot.df = grep_result$replicate_stats %>%
    dplyr::select(replicate, type, num_reads, "HDR reads" = num_hdr_reads, "WT reads" = num_wt_reads, HDR_WT_ratio, HDR_frac, WT_frac)
  stats.plot.df$HDR_WT_ratio = sprintf("%.3g", stats.plot.df$HDR_WT_ratio)
  stats.plot.df$HDR_frac = sprintf("%.2f%%", 100 * stats.plot.df$HDR_frac)
  stats.plot.df$WT_frac = sprintf("%.2f%%", 100 * stats.plot.df$WT_frac)
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())

  mycex = min(0.70, 0.75 * 11 / nrow(stats.plot.df))
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = mycex)),
                                            colhead = list(fg_params=list(cex = mycex)),
                                            rowhead = list(fg_params=list(cex = mycex)),
                                            padding = unit(c(2, 4), "mm"))
  plot_title = sprintf("%s grep summary", grep_result$region$name)
  p.stats = ggplot(tibble(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(gridExtra::tableGrob(t(stats.plot.df), theme = myTableTheme), xmin=1.2, xmax=10, ymin=1, ymax=8.5) +
    annotate("text", x=1, y=10, label = summary.left, vjust = 1, hjust = 0, size = 2.9) +
    annotate("text", x=9.5, y=10, label = summary.right, vjust = 1, hjust = 1, size = 2.9)

  color_values = c(`cDNA`="firebrick1", `cDNA outlier`="orange1", `gDNA`="dodgerblue3", `gDNA outlier`="turquoise1")
  p.num_reads = ggplot(grep_result$replicate_stats, aes(x=replicate, y=num_reads, fill=type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%d", num_reads)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), legend.title=element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = color_values) +
    ylab("Number of reads") +
    ggtitle("Number of reads")

  p.hdr_wt = ggplot(grep_result$replicate_stats, aes(x=replicate, y=HDR_WT_ratio, fill=type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3g", HDR_WT_ratio)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank()) +
    scale_fill_manual(values = color_values, guide=F) +
    ylab("HDR:WT ratio") +
    ggtitle("HDR:WT ratio")

  p.title = cowplot::ggdraw() + cowplot::draw_label(plot_title, fontface='bold')
  p.res = egg::ggarrange(p.title, p.stats, p.num_reads, p.hdr_wt, ncol=1, heights=c(0.5, 4.5, 1.5, 1.5), draw = F)
  return(p.res)
}


#' Returns all main plots for a single deletion analysis result.
#'
#' @param del_result Result from a call to deletion_analysis.
#' @param opts A list with all options needed for rgenie plotting functions.
#' @param variance_components_plot If TRUE, then variance_components_plot() is called.
#' @param power_plots If TRUE, then power_plots() is called.
#' @return Returns a list of ggplot objects.
#'
#' @examples
#' \dontrun{
#' del_result = deletion_analysis(regions, replicates)
#' deletion_plots(del_result, genie_plot_options())
#' genie_plot_options
#' }
#' @seealso \code{\link{genie_plot_options}}
#' @seealso \code{\link{deletion_analysis}}
#' @seealso \code{\link{deletion_summary_plot}}
#' @seealso \code{\link{deletion_alleles_plot}}
#' @seealso \code{\link{deletion_profile_plot}}
#' @seealso \code{\link{replicate_summary_plot}}
#' @seealso \code{\link{replicate_qc_plot}}
#' @seealso \code{\link{allele_effect_plot}}
#' @seealso \code{\link{variance_components_plot}}
#' @seealso \code{\link{power_plots}}
#' @export
#'
deletion_plots = function(del_result,
                          opts = genie_plot_options(),
                          variance_components_plot = F,
                          power_plots = F) {
  check_is_del_result(del_result)
  plot_list = list( deletion_summary = deletion_summary_plot(del_result),
                    deletion_alleles = deletion_alleles_plot(del_result, viewing_window = opts$viewing_window, color_by = opts$deletion_alleles_plot_color_by),
                    deletion_profile = deletion_profile_plot(del_result, viewing_window = opts$viewing_window),
                    replicate_summary = replicate_summary_plot(del_result, outlier_threshold = opts$outlier_threshold),
                    replicate_qc = replicate_qc_plot(del_result, outlier_threshold = opts$outlier_threshold),
                    allele_effect = allele_effect_plot(del_result, viewing_window = opts$viewing_window, max_alleles = opts$allele_plot_max_alleles) )

  if (variance_components_plot) {
    vc = get_variance_components(del_result, replicates,
                                 allele_min_reads = opts$variance_analysis_min_count,
                                 allele_min_fraction = opts$variance_analysis_min_fraction)
    vc_plot = variance_components_plot(vc, split_by_fraction = opts$variance_plot_split_by_fraction)
    plot_list$variance_components = vc_plot
  }
  if (power_plots) {
    pwr = power_plots(del_result,
                      allele_min_reads = opts$plower_plot_min_count,
                      WT_fraction = opts$plower_plot_WT_fraction)

    plot_list = c(plot_list, pwr)
  }
  return(plot_list)
}


#' Plots a summary of deletion analysis results for a single region.
#'
#' @param del_result Result from a call to deletion_analysis.
#' @return Returns a ggplot object with a summary of deletion analysis results for a single region.
#'
#' @examples
#' \dontrun{
#' del_result = deletion_analysis(regions, replicates)
#' deletion_summary_plot(del_result)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
deletion_summary_plot = function(del_result) {
  check_is_del_result(del_result)
  stats.plot.df = del_result$replicate_stats %>%
    dplyr::select(replicate, type, num_udps, HDR_WT_ratio, DEL_WT_ratio,
                  HDR_rate, DEL_rate, editing_rate, WT_rate, num_reads, num_hdr_reads, num_wt_reads,
                  num_deletion_reads, num_insertion, reads_excluded_for_minoverlap, reads_excluded_for_mismatches,
                  reads_excluded_nonspanning, reads_excluded_for_multiple_deletions)
  stats.plot.df = stats.plot.df %>% dplyr::rename("HDR reads" = num_hdr_reads,
                                                  "WT reads" = num_wt_reads,
                                                  "Deletion reads" = num_deletion_reads,
                                                  "excluded-insertion" = num_insertion,
                                                  "excluded-minoverlap" = reads_excluded_for_minoverlap,
                                                  "excluded-mismatches" = reads_excluded_for_mismatches,
                                                  "excluded-nonspanning" = reads_excluded_nonspanning,
                                                  "excluded-mult.deletions" = reads_excluded_for_multiple_deletions)
  # Convert to strings for nice printing (with 3 significant digits)
  stats.plot.df$HDR_WT_ratio = sprintf("%.3g", stats.plot.df$HDR_WT_ratio)
  stats.plot.df$DEL_WT_ratio = sprintf("%.3g", stats.plot.df$DEL_WT_ratio)
  stats.plot.df$HDR_rate = sprintf("%.2f%%", 100 * stats.plot.df$HDR_rate)
  stats.plot.df$DEL_rate = sprintf("%.2f%%", 100 * stats.plot.df$DEL_rate)
  stats.plot.df$editing_rate = sprintf("%.2f%%", 100 * stats.plot.df$editing_rate)
  stats.plot.df$WT_rate = sprintf("%.2f%%", 100 * stats.plot.df$WT_rate)

  stats.gDNA = del_result$replicate_stats %>% filter(type == "gDNA")
  stats.cDNA = del_result$replicate_stats %>% filter(type == "cDNA")
  hdr.rate = sprintf("Mean HDR frac gDNA: %.2g%%,  cDNA: %.2g%%", mean(stats.gDNA$HDR_rate) * 100, mean(stats.cDNA$HDR_rate) * 100)
  del.rate = sprintf("Mean DEL frac gDNA: %.2g%%,  cDNA: %.2g%%", mean(stats.gDNA$DEL_rate) * 100, mean(stats.cDNA$DEL_rate) * 100)
  wt.rate = sprintf("Mean WT frac gDNA: %.2g%%,  cDNA: %.2g%%", mean(stats.gDNA$WT_rate) * 100, mean(stats.cDNA$WT_rate) * 100)

  if (is.null(del_result$region_stats)) {
    summary.left = "Unable to calculate stats with < 2 replicates"
    summary.right = ""
  } else {
    confIntervalString = function(effect_confint_lo, effect_confint_hi) {
      sprintf("95%% CI: (%.3g, %.3g)", effect_confint_lo, effect_confint_hi)
    }
    hdr.conf.interval.str = confIntervalString(del_result$region_stats$hdr_effect_confint_lo, del_result$region_stats$hdr_effect_confint_hi)
    hdr.summary = sprintf("cDNA:gDNA ratio (HDR/WT): %.3g\n%s,    p = %.3g",
                          del_result$region_stats$hdr_effect, hdr.conf.interval.str, del_result$region_stats$hdr_pval)

    summary.left = paste(hdr.rate, del.rate, wt.rate, "", hdr.summary, sep = "\n")

    del.conf.interval.str = confIntervalString(del_result$region_stats$del_effect_confint_lo, del_result$region_stats$del_effect_confint_hi)
    del.summary = sprintf("cDNA:gDNA ratio (DEL/WT): %.3g\n%s,    p = %.3g",
                          del_result$region_stats$del_effect, del.conf.interval.str, del_result$region_stats$del_pval)

    region_length = del_result$region$end - del_result$region$start + 1
    window_start = min(region_length, max(1, del_result$region$highlight_site + del_result$opts$del_span_start))
    window_end = min(region_length, max(1, del_result$region$highlight_site + del_result$opts$del_span_end))
    del_window.conf.interval.str = confIntervalString(del_result$region_stats$del_window_effect_confint_lo, del_result$region_stats$del_window_effect_confint_hi)
    del.summary.window = sprintf("cDNA:gDNA ratio (DEL/WT) [%d-%d]: %.3g\n%s,    p = %.3g",
                                 window_start, window_end,
                                 del_result$region_stats$del_window_effect, del_window.conf.interval.str, del_result$region_stats$del_window_pval)

    summary.right = paste(del.summary, "", del.summary.window, sep = "\n")
  }

  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  mycex = min(0.70, 0.75 * 11 / nrow(stats.plot.df))
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = mycex)),
                                            colhead = list(fg_params=list(cex = mycex)),
                                            rowhead = list(fg_params=list(cex = mycex)),
                                            padding = unit(c(2, 3), "mm"))
  p.stats = ggplot(tibble(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(gridExtra::tableGrob(t(stats.plot.df), theme = myTableTheme), xmin=1.2, xmax=10, ymin=1, ymax=8) +
    annotate("text", x=5, y=10, label = sprintf("%s analysis summary", del_result$region$name), vjust = 1, fontface = 2, size = 5) +
    annotate("text", x=1, y=9.4, label = summary.left, vjust = 1, hjust = 0, size = 2.9) +
    annotate("text", x=9.5, y=9.4, label = summary.right, vjust = 1, hjust = 1, size = 2.9)

  return(p.stats)
}


#' Plots a summary of deletion analysis replicates.
#'
#' @param del_result Result from a call to deletion_analysis.
#' @param outlier_threshold A numeric threshold for the outlier score, above which replicates will be colored differently.
#' @return Returns a ggplot object with a summary of deletion analysis replicates.
#'
#' @examples
#' \dontrun{
#' del_results = deletion_analysis(regions, replicates)
#' replicate_summary_plot(del_results)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
replicate_summary_plot = function(del_result,
                                  outlier_threshold = NA) {
  check_is_del_result(del_result)
  stats.df = del_result$replicate_stats
  color_values = c(`cDNA`="firebrick1", `cDNA outlier`="orange1", `gDNA`="dodgerblue3", `gDNA outlier`="turquoise1")
  stats.df$replicate = forcats::fct_reorder(stats.df$replicate, as.integer(factor(stats.df$type)))
  stats.df = stats.df %>% group_by(type) %>%
    rowwise() %>%
    mutate(type2 = if_else(!is.na(outlier_threshold) && !is.na(outlier_score_knn) && outlier_score_knn > outlier_threshold, paste(type, "outlier"), type)) %>%
    arrange(replicate)
  p1 = ggplot(stats.df, aes(x=replicate, y=num_reads, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%d", num_reads)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), legend.title=element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = color_values) +
    ylab("Number of reads") +
    ggtitle("Number of reads")

  p2 = ggplot(stats.df, aes(x=replicate, y=num_udps, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%d", num_udps)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), legend.title=element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = color_values, guide=F) +
    ylab("Number of UDPs") +
    ggtitle("Number of UDPs")

  p3 = ggplot(stats.df, aes(x=replicate, y=HDR_WT_ratio, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3g", HDR_WT_ratio)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), legend.title=element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = color_values, guide=F) +
    ylab("HDR:WT ratio") +
    ggtitle("HDR:WT ratio")

  p4 = ggplot(stats.df, aes(x=replicate, y=DEL_WT_ratio, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3g", DEL_WT_ratio)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank()) +
    scale_fill_manual(values = color_values, guide=F) +
    ylab("Del:WT ratio") +
    ggtitle("Deletion:WT ratio")

  p.title = cowplot::ggdraw() + cowplot::draw_label(sprintf("%s replicate summary", del_result$region$name), fontface='bold')
  #p.replicate_qc = cowplot::plot_grid(p.udp_fractions, p.udp_avg_deviation, ncol=1)
  #p.res = cowplot::plot_grid(p.title, p.replicate_qc, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  p.res = egg::ggarrange(p.title, p1, p2, p3, p4, ncol=1, heights=c(1,3,3,3,3), draw = F)
  p.res
}


#' Plots quality control metrics for deletion analysis replicates.
#'
#' @param del_result Result from a call to deletion_analysis.
#' @param outlier_threshold A numeric threshold for the outlier score, above which replicates will be colored differently.
#' @return Returns a ggplot object with quality control metrics for deletion analysis replicates.
#'
#' @examples
#' \dontrun{
#' del_results = deletion_analysis(regions, replicates)
#' replicate_qc_plot(del_results)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
replicate_qc_plot = function(del_result,
                             outlier_threshold = NA) {
  check_is_del_result(del_result)
  stats.df = del_result$replicate_stats
  replicate.udp.fractions = del_result$replicate_allele_fractions
  replicate.udp.fractions$replicate = forcats::fct_reorder(replicate.udp.fractions$replicate, as.integer(factor(replicate.udp.fractions$type)))

  color_values = c(`cDNA`="firebrick1", `cDNA outlier`="orange1", `gDNA`="dodgerblue3", `gDNA outlier`="turquoise1")
  shapes = c(16,17,15,3,12,8,6,5,0,1,11,10,18,7,9,2,3,4,13,14)
  num_gDNA_reps = length(unique(replicate.udp.fractions %>% dplyr::filter(type == "gDNA") %>% .$replicate))
  num_cDNA_reps = length(unique(replicate.udp.fractions %>% dplyr::filter(type == "cDNA") %>% .$replicate))
  shape_values = c(shapes[1:min(20, num_cDNA_reps)], shapes[1:min(20, num_gDNA_reps)])
  if (length(unique(replicate.udp.fractions$replicate)) <= 20) {
    p.udp_fractions = ggplot(replicate.udp.fractions, aes(x=udp_id, y=udp_fraction, group=replicate, col=type, shape=replicate))
  } else {
    p.udp_fractions = ggplot(replicate.udp.fractions, aes(x=udp_id, y=udp_fraction, group=replicate, col=type))
  }
  num_udps = length(unique(replicate.udp.fractions$udp))
  p.udp_fractions = p.udp_fractions +
    geom_line() + geom_point(size=3, alpha=0.5) +
    scale_color_manual(values = color_values, guide = F) +
    scale_shape_manual(values = shape_values) +
    theme_bw(10) + xlab("UDP") + ylab("UDP fraction") +
    scale_y_log10() +
    ggtitle(sprintf("UDP fractions (min %.2g%%, %d UDPs)", del_result$opts$min_avg_allele_fraction*100, num_udps))

  stats.df = stats.df %>% group_by(type) %>%
    rowwise() %>%
    mutate(type2 = if_else(!is.na(outlier_threshold) && !is.na(outlier_score_knn) && outlier_score_knn > outlier_threshold, paste(type, "outlier"), type)) %>%
    arrange(replicate)
  stats.df$replicate = forcats::fct_reorder(stats.df$replicate, as.integer(factor(stats.df$type)))

  p.udp_avg_deviation = ggplot(stats.df %>% filter(!is.na(avg_udp_deviation)),
                               aes(x=replicate, y=avg_udp_deviation, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2g%%", avg_udp_deviation*100)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(), legend.position="none") +
    scale_fill_manual(values = color_values) +
    ylab("UDP frac deviation") +
    ggtitle("Mean UDP fraction deviation")

  p.udp_avg_deviation_gDNA = ggplot(stats.df %>% filter(!is.na(avg_udp_deviation_from_gDNA)),
                                    aes(x=replicate, y=avg_udp_deviation_from_gDNA, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2g%%", avg_udp_deviation_from_gDNA*100)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.title=element_blank()) +
    scale_fill_manual(values = color_values) +
    ylab("UDP frac deviation") +
    ggtitle("Mean UDP fraction deviation (compared to gDNA)")

  # p.udp_avg_deviation_relative = ggplot(stats.df, aes(x=replicate, y=avg_udp_rel_deviation_from_gDNA, fill=type2)) +
  #   geom_bar(stat = "identity", alpha = 0.8) +
  #   geom_text(aes(label = sprintf("%.0f%%", avg_udp_rel_deviation_from_gDNA*100)), size = 2.4, position = position_stack(vjust = 0.5)) +
  #   theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank()) +
  #   scale_fill_manual(values = color_values) +
  #   ylab("UDP rel. deviation") +
  #   ggtitle("Mean UDP fraction relative deviation (rel. to gDNA UDP mean)")

  outlier.df = stats.df
  if (all(is.na(outlier.df$outlier_score_knn))) {
    p.outlier_scores = ggplot(outlier.df, aes(x=replicate, y=outlier_score_knn)) +
      geom_blank() +
      theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank(), legend.position="none") +
      xlab("Replicate") + ylab("Outlier score") +
      ggtitle("KNN outlier score") +
      annotate("text", x=outlier.df$replicate[1], y=NA, label="Outlier scores are only calculated with sufficient UDPs and >= 4 replicates", hjust = 0)
  } else {
    p.outlier_scores = ggplot(outlier.df, aes(x=replicate, y=outlier_score_knn, fill=type2)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.3g", outlier_score_knn)), size = 2.4, position = position_stack(vjust = 0.5)) +
      theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank(), legend.position="none") +
      scale_fill_manual(values = color_values) +
      xlab("Replicate") + ylab("Outlier score") +
      ggtitle("KNN outlier score")
  }

  p.title = cowplot::ggdraw() + cowplot::draw_label(sprintf("%s replicate QC", del_result$region$name), fontface='bold')
  #p.replicate_qc = cowplot::plot_grid(p.udp_fractions, p.udp_avg_deviation, ncol=1)
  #p.replicate_qc = cowplot::plot_grid(p.title, p.replicate_qc, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  p.replicate_qc = egg::ggarrange(p.title, p.udp_fractions, p.udp_avg_deviation, p.udp_avg_deviation_gDNA, p.outlier_scores, ncol=1, heights=c(1.5,4,4,4,4), draw = F)

  return( p.replicate_qc )
}


#' Plots estimated effect sizes and confidence intervals for top alleles from a deletion analysis.
#'
#' @param del_result Result from a call to deletion_analysis.
#' @param viewing_window Window on either size of the CRISPR cut site to show in the plot.
#' @param max_alleles The maximum number of alleles to show in the plot.
#' @return Returns a ggplot object plotting effect sizes and confidence intervals for top alleles
#'  from a deletion analysis.
#' Top alleles are in decreasing order of their total read count in gDNA across replicates.
#'
#' @examples
#' \dontrun{
#' del_results = deletion_analysis(regions, replicates)
#' allele_effect_plot(del_results)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
allele_effect_plot = function(del_result,
                              viewing_window = 40,
                              max_alleles = 40) {
  check_is_del_result(del_result)
  if (viewing_window < 1) {
    warning("viewing_window should be a positive integer.")
    return(NULL)
  }
  if (max_alleles > 200) {
    warning("Setting max_alleles to 200.")
    max_alleles = 200
  }
  cut_site = del_result$region$cut_site
  view_start = max(1, cut_site - viewing_window)
  view_end = cut_site + viewing_window

  plot_title = sprintf("%s effect estimates", del_result$region$name)

  udp.dels.df = del_result$allele_effect %>%
    dplyr::filter(!is_wt_allele) %>%
    dplyr::arrange(!is_hdr_allele, desc(total_count))

  numUDPs = nrow(udp.dels.df)
  if (numUDPs < 2) {
    return(egg::ggarrange(text_plot("No UDPs pass the thresholds for min read counts."), top=plot_title, draw = F))
  }
  if (!is.na(max_alleles) & numUDPs > max_alleles) {
    udp.dels.df = udp.dels.df %>% .[1:max_alleles,]
    numUDPs = max_alleles
  }

  # matrix containing the deletion binary code
  M = as.matrix(do.call(rbind, lapply(as.list(udp.dels.df$udp), udp_to_binary)))
  M_chars = as.matrix(do.call(rbind, lapply(as.list(udp.dels.df$udp), str_to_chars)))

  model <- hclust(dist(M))
  dhc <- as.dendrogram(model)
  ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
  dendro_span = max(ddata$segments$y, ddata$segments$yend) - min(ddata$segments$y, ddata$segments$yend)

  udp.order.map = tibble(y = ggdendro::label(ddata)$x, udp = udp.dels.df[model$order,]$udp)
  udp.dels.df = udp.dels.df[model$order,] %>% dplyr::left_join(udp.order.map, by="udp")

  plot.df = as_tibble(M_chars[model$order,])
  colnames(plot.df) = as.character(c(1:ncol(plot.df)))
  profileSpan = ncol(plot.df)
  plot.df$id = c(1:nrow(plot.df))
  plot.df[,"HDR allele"] = udp.dels.df$is_hdr_allele

  plot.gather.df = plot.df %>% tidyr::gather(key = "position", value = "udpchar", 1:profileSpan)
  plot.gather.df$pos = as.numeric(plot.gather.df$position)

  p.dendro = ggplot() +
    geom_segment(aes(x = ggdendro::segment(ddata)$x, y = ggdendro::segment(ddata)$y, xend = ggdendro::segment(ddata)$xend, yend = ggdendro::segment(ddata)$yend)) +
    coord_cartesian(xlim=c(min(plot.gather.df$id), max(plot.gather.df$id)), ylim=c(view_start, view_end), default = T) +
    ylab("") + scale_y_continuous(expand = c(0, 0), trans = "reverse") +
    coord_flip() +
    theme_bw() + theme(axis.title.y = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line = element_blank(),
                       panel.border = element_blank(),
                       panel.grid = element_blank(),
                       plot.margin = unit(c(0,0,0.05,0), "cm"))

  numBases = view_end - view_start + 1
  dash_size = 2
  dot_size = 0.5
  if (numBases > 60) {
    dash_size = 1.2
    dot_size = 0.4
  }
  if (numBases > 101) {
    dash_size = 0.8
    dot_size = 0.25
  }
  if (any(plot.gather.df$`HDR allele`)) {
    p.udp = ggplot() +
      geom_rect(aes(xmin=min(pos), xmax=max(pos), ymin=(id-0.5), ymax=(id+0.5)), fill = "palegreen", data = plot.gather.df[plot.gather.df$`HDR allele`,]) +
      geom_point(aes(x = pos, y = id), size = dot_size, shape = 19, alpha = 0.7, data = plot.gather.df[plot.gather.df$udpchar == '-' & plot.gather.df$`HDR allele`,]) +
      geom_point(aes(x = pos, y = id), size = dash_size, shape = '-', alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '*' & plot.gather.df$`HDR allele`,])
  } else {
    p.udp = ggplot()
  }
  p.udp = p.udp +
    geom_point(aes(x = pos, y = id), size = dot_size, shape = 19, alpha = 0.7, data = plot.gather.df[plot.gather.df$udpchar == '-',]) +
    geom_point(aes(x = pos, y = id), size = dash_size, shape = '-', alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '*',]) +
    scale_color_manual(values = c("grey20", "green4")) +
    geom_point(aes(x = pos, y = id, shape = udpchar), color = "black", size = 2.5, data = plot.gather.df[!(plot.gather.df$udpchar == '*' | plot.gather.df$udpchar == '-'),]) + scale_shape_identity() +
    coord_cartesian(ylim=c(min(plot.gather.df$id), max(plot.gather.df$id)), xlim=c(view_start, view_end)) +
    xlab("Position") + scale_x_continuous(expand = c(0.01, 0)) +
    theme_bw() + theme(legend.position = "none",
                       axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.line = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit(c(0,0,0.05,0), "cm"))
  # Used to colour the HDR allele text, but now we show a rectangle around HDR alleles
  #geom_point(aes(x = pos, y = id, colour = `HDR allele`), size = dot_size, shape = 19, alpha = 0.7, data = plot.gather.df[plot.gather.df$udpchar == '-',]) +

  if (!is.na(cut_site)) {
    p.udp = p.udp + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }

  # Left side of the plot, showing UDP-normalised scores (UNS)
  uns.plot.df = data.frame(y = ggdendro::label(ddata)$x,
                           x = 0,
                           uns = udp.dels.df$uns,
                           cDNA_lower = udp.dels.df$uns < 1,
                           gDNA = udp.dels.df$gDNA_total_count,
                           uns_conf_hi = udp.dels.df$uns_confint_hi,
                           uns_conf_lo = udp.dels.df$uns_confint_lo) %>%
    rowwise() %>%
    mutate(uns_conf_hi = if_else(uns_conf_hi < 0.001, NaN, uns_conf_hi),
           uns_conf_lo = if_else(uns_conf_lo < 0.001, NaN, uns_conf_lo))

  maxDotSize = 5.5
  if (numUDPs > 60) {
    maxDotSize = 3.5
  } else if (numUDPs > 30) {
    maxDotSize = 4.5
  }

  min_uns_threshold = 0.25
  uns.plot.df$uns = sapply(uns.plot.df$uns, FUN = function(x) max(x, min_uns_threshold))
  uns.plot.df$uns_conf_lo = sapply(uns.plot.df$uns_conf_lo, FUN = function(x) max(x, 0.01))
  max_uns_display = max(1.5, ceiling(max(uns.plot.df$uns, na.rm = T)))
  min_uns_display = min(uns.plot.df$uns, na.rm = T)

  # Code testing different visualisation for representing confidence in UNS value
  #log2_uns_range = log2(max(uns.plot.df$uns)) - log2(min(uns.plot.df$uns))
  #max_log_gDNA = max(log10(uns.plot.df$gDNA))
  #uns.plot.df$confidence = sapply(1:nrow(uns.plot.df), FUN=function(i) min(log10(max(uns.plot.df[i,]$gDNA-10, 10)) / max_log_gDNA, 1) * abs(log2(uns.plot.df[i,]$uns)) / log2_uns_range)
  #uns.plot.df$confidence = (uns.plot.df$alpha - min(uns.plot.df$alpha)) / (max(uns.plot.df$alpha) - min(uns.plot.df$alpha))
  uns.plot.df$cDNA_expr = "higher"
  uns.plot.df$cDNA_expr[uns.plot.df$cDNA_lower] = "lower"

  p.uns = ggplot(uns.plot.df) +
    #annotate("segment", x = uns.plot.df$x, xend = uns.plot.df$uns, y = uns.plot.df$y, yend = uns.plot.df$y, colour = 'gray') +
    geom_errorbarh(mapping = aes(y = y, xmin = uns_conf_lo, xmax = uns_conf_hi, height = 0.4),
                   data = uns.plot.df %>% filter(!is.nan(uns_conf_hi), !is.nan(uns_conf_lo)),
                   colour = "grey70") +
    geom_point(aes(x = uns, y = y, colour = cDNA_expr, size = log10(gDNA), alpha = log10(gDNA))) +
    scale_size(range = c(1, maxDotSize)) +
    scale_color_manual(values=c(higher="red", lower="blue"), name="cDNA expr") +
    scale_alpha_continuous(range = c(0.2, 1)) +
    scale_x_continuous(trans="log2", breaks = c(0.5, 1, 2, 4, 8)) +
    #scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, limits = c(0.1, 4)) +
    #geom_point(data = uns.replicates.plot.df, mapping = aes(y = y, x = uns, size = 2, alpha = 0.6)) +
    scale_y_continuous(position = "right") +
    coord_cartesian(xlim = c(min_uns_display, max_uns_display),
                    ylim=c(min(plot.gather.df$id), max(plot.gather.df$id))) +
    ylab("Unique deletion profile") + xlab("Effect") +
    geom_vline(xintercept = 1.0, alpha = 0.25, linetype = "longdash") +
    theme_bw() + theme(legend.position = "right",
                       #axis.line = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit(c(0,0,0.05,0), "cm"))

  p.udp_profile = egg::ggarrange(p.dendro, p.udp, p.uns, nrow=1, ncol=3, widths=c(0.6,4,1), top = plot_title, draw = F)
  #p.udp_profile = cowplot::plot_grid(p.uns, p.udp, p.dendro, nrow=1, ncol=3, rel_widths = c(2,4,1))
  return(p.udp_profile)
}


variance_partition_plot = function(vp.df, residuals=T, pointColor = NA) {
  vp.df = vp.df %>%
    tidyr::gather(variable, variance, -udp, -name, -frac) %>%
    mutate(pctvariance = variance*100)

  if (!residuals) {
    vp.df = vp.df %>% dplyr::filter(variable != "Residuals")
  }

  medians = vp.df %>% group_by(variable) %>% summarise(median=median(pctvariance)) %>%
    arrange(-median) %>% filter(variable != "Residuals")
  vp.df$variable = factor(as.character(vp.df$variable),
                          levels=c(medians$variable, "Residuals"))
  numVariables = length(levels(vp.df$variable))

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  }
  getFractionCategory = function(x) {
    if (x < 0.005) { "< 0.5%" }
    else if (x < 0.02) { "< 2%" }
    else { "> 2%" }
  }
  vp.df$fraction = sapply(vp.df$frac, getFractionCategory)

  plotaes = aes(color = "blue")
  if (!is.na(pointColor) && pointColor == "fraction") {
    plotaes = aes(color = fraction)
  }
  ptSize = 1.8
  ptAlpha = 0.6
  if (nrow(vp.df) > 300) {
    ptSize = 1
    ptAlpha = 0.35
  }
  p = ggplot(vp.df, aes(x=variable, y=pctvariance, fill="grey")) +
    geom_violin(scale="width") +
    geom_boxplot(width=0.2, fill="grey80", outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = ptAlpha, size=ptSize, mapping = plotaes) +
    theme_bw() + guides(fill=FALSE) +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(angle=20, hjust = 1), axis.title.x=element_blank(), panel.grid=element_blank()) +
    ylab("Variance explained (%)") +
    ylim(c(-0.01,100.01)) +
    scale_fill_manual(values = c("gray95"))
  #  scale_fill_manual(values = c(gg_color_hue(numVariables-1), "gray85"))
  if (is.na(pointColor)) {
    p = p + scale_color_manual(values = c(blue="blue"), guide = F)
  } else if (pointColor == "fraction") {
    #p = p + scale_color_manual(values = c(rep(color, numVariables)))
    p = p + scale_color_manual(values = c("< 0.5%"="chartreuse", "< 2%"="blue", "> 2%"="red2"))
  }
  return(p)
}


#' Plots variance components estimates for all unique alleles.
#'
#' @param varcomp Result from a call to get_variance_components.
#' @param split_by_fraction If TRUE, then points are colored by allele fraction.
#' @return Returns a ggplot object.
#'
#' @examples
#' \dontrun{
#' del_results = deletion_analysis(regions, replicates)
#' vc = get_variance_components(del_result[[1]], replicates)
#' variance_components_plot(vc)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @seealso \code{\link{get_variance_components}}
#' @export
#'
variance_components_plot = function(varcomp,
                                    split_by_fraction = F) {
  pointColor = NA
  if (split_by_fraction) {
    pointColor = "fraction"
  }
  if (nrow(varcomp$vp_cDNA) <= 0) {
    p.cDNA = text_plot("Cannot make variance components plot.")
    cDNA.summary.df = data.frame(x="empty")
    cDNA_note = ""
  } else {
    p.cDNA = variance_partition_plot(varcomp$vp_cDNA, pointColor=pointColor) + ggtitle("cDNA")
    cDNA.summary.df = variance_partition_summary_table(varcomp$vp_cDNA %>% dplyr::select(-name, -udp, -frac))
    cDNA_note = sprintf("Based on %d alleles with >%d reads, >%.1g%% allele fraction",
                        nrow(varcomp$vp_cDNA), varcomp$opts$allele_min_reads, varcomp$opts$allele_min_fraction*100)
  }

  if (nrow(varcomp$vp_gDNA) <= 0) {
    p.gDNA = text_plot("Cannot make variance components plot.")
    gDNA.summary.df = data.frame(x="empty")
    gDNA_note = ""
  } else {
    p.gDNA = variance_partition_plot(varcomp$vp_gDNA, pointColor=pointColor) + ggtitle("gDNA") + theme(legend.position="left")
    gDNA.summary.df = variance_partition_summary_table(varcomp$vp_gDNA %>% dplyr::select(-name, -udp, -frac))
    gDNA_note = sprintf("Based on %d alleles with >%d reads, >%.1g%% alleles fraction",
                        nrow(varcomp$vp_gDNA), varcomp$opts$allele_min_reads, varcomp$opts$allele_min_fraction*100)
  }

  ggThemeBlank = theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                    axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
                                    panel.border = element_blank())
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.75)),
                                            colhead = list(fg_params=list(cex = 0.75)),
                                            rowhead = list(fg_params=list(cex = 0.75)))

  plot_title = sprintf("%s Variance components", varcomp$name)
  p.stats = ggplot(tibble(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(gridExtra::tableGrob(cDNA.summary.df, theme = myTableTheme), xmin=1, xmax=5, ymin=2, ymax=8.5) +
    annotation_custom(gridExtra::tableGrob(gDNA.summary.df, theme = myTableTheme), xmin=6, xmax=10, ymin=2, ymax=8.5) +
    # annotate("text", x=3.5, y=8, label = "cDNA", hjust = 0.5, fontface = 2, size = 5) +
    # annotate("text", x=7.5, y=8, label = "gDNA", hjust = 0.5, fontface = 2, size = 5) +
    annotate("text", x=5, y=10, label = plot_title, vjust = 1, fontface = 2, size = 5) +
    annotate("text", x=1, y=1, label = cDNA_note, hjust = 0, vjust = 0, size = 3) +
    annotate("text", x=6, y=1, label = gDNA_note, hjust = 0, vjust = 0, size = 3)

  return( cowplot::plot_grid(p.stats,
                             cowplot::plot_grid(p.cDNA, p.gDNA,
                                                nrow = 1, rel_widths=c(0.46, 0.54)),
                             nrow = 2, rel_heights = c(0.4, 0.6)) )
}

variance_partition_summary_table = function(vp.df) {
  # Get some summary stats for each column of the variance partition data.frame
  col = colnames(vp.df)[1]
  summary_df = summary(vp.df[, col, drop=T])
  for (col in colnames(vp.df)[-1]) {
    summary_df = rbind(summary_df, summary(vp.df[, col, drop=T]))
  }
  summary_df = as.data.frame(summary_df)
  rownames(summary_df) = colnames(vp.df)
  for (col in colnames(summary_df)) {
    summary_df[,col] = sapply(summary_df[,col,drop=T], FUN = function(x) sprintf("%.3g", x))
  }
  summary_df[, c("1st Qu.", "Median", "3rd Qu.")]
}


text_plot = function(text, title = NA, fontface = "plain", size = 4) {
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  p = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) +
    geom_blank() + ggThemeBlank +
    annotate("text", x=5.5, y=5.5, label = text, fontface = "bold", size = size, hjust = "center")
  if (!is.na(title)) {
    p = p + ggtitle(title)
  }
  return(p)
}


#' Plots unique deletion alleles and their "pileup" count profile separately for cDNA and gDNA.
#'
#' @param del_result Result from a call to deletion_analysis.
#' @param viewing_window Window on either size of the CRISPR cut site to show in the plot.
#' @param color_by A string with one of the values: "none" (default), "window", or "sharing".
#' @return Returns a ggplot object.
#'
#' @examples
#' \dontrun{
#' del_results = deletion_analysis(regions, replicates)
#' deletion_alleles_plot(del_results)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
deletion_alleles_plot = function(del_result,
                                 viewing_window = 40,
                                 color_by="window") {
  check_is_del_result(del_result)
  stopifnot(color_by %in% c("none", "sharing", "window"))
  p.udp_gDNA = udp_plot(del_result$region_alleles %>% dplyr::filter(type == "gDNA"), plot_title="gDNA", del_result$region$cut_site, del_result$region$highlight_site, viewing_window, color_by)
  p.udp_cDNA = udp_plot(del_result$region_alleles %>% dplyr::filter(type == "cDNA"), plot_title="cDNA", del_result$region$cut_site, del_result$region$highlight_site, viewing_window, color_by)
  p.title = cowplot::ggdraw() + cowplot::draw_label(sprintf("%s deletion alleles", del_result$region$name), fontface='bold')
  p.cDNA_gDNA = cowplot::plot_grid(p.udp_gDNA, p.udp_cDNA, nrow=1)
  p.merged_UDP = cowplot::plot_grid(p.title, p.cDNA_gDNA, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  return(p.merged_UDP)
}


udp_plot = function(udp.df, plot_title, cut_site, highlight_site, viewing_window = 40, color_by="none") {
  if (viewing_window < 1) {
    warning("udp_plot: viewing_window should be a positive integer.")
    return(NULL)
  }
  view_start = max(1, cut_site - viewing_window)
  view_end = cut_site + viewing_window

  # We allow there to be up to 2 deletions in the UDP
  udp.df = udp.df %>% filter(has_crispr_deletion, !is.na(deletion_start))
  if (nrow(udp.df) == 0) {
    return(egg::ggarrange(text_plot("No UDPs to plot."), top=plot_title, draw = F))
  }
  if (color_by == "none") {
    udp.df$udp_sharing = "both"
  }
  udp.dels.df = udp.df %>% dplyr::select(deletion_start, deletion_end, deletion2_start, deletion2_end, has_deletion_in_window, sharing=udp_sharing)
  udp.dels.df = udp.dels.df %>% arrange(deletion_start, deletion2_start)
  udp.dels.df$y = 1:nrow(udp.dels.df)
  udp.plot.df = bind_rows(udp.dels.df %>% dplyr::select(y, deletion_start, deletion_end, has_deletion_in_window, sharing),
                          udp.dels.df %>% dplyr::select(y, deletion_start=deletion2_start, deletion_end=deletion2_end, has_deletion_in_window, sharing))
  udp.plot.df = udp.plot.df %>% dplyr::filter(!is.na(deletion_start))

  xmax = nchar(udp.df$udp[1])
  segment_size = 0.5
  if (nrow(udp.dels.df) > 200) {
    segment_size = 0.3
  }

  simple_theme = theme_bw() + theme(panel.grid.minor = element_line(colour="grey90", size=0.1), panel.grid.major = element_line(colour="grey90", size=0.1))
  highlight_color = "orange"
  udp_colors = c("cDNA only"=highlight_color, "gDNA only"=highlight_color, "both"="dodgerblue3", "unclear"="grey", "in_window"="red")
  if (color_by == "window") {
    highlight_color = "red"
    p.udp_dist = ggplot(udp.plot.df) +
      geom_segment(aes(x = deletion_start, xend = deletion_end, y = y, yend = y, color = has_deletion_in_window), size = segment_size) +
      simple_theme + coord_cartesian(xlim = c(view_start, view_end)) + scale_y_reverse() +
      scale_color_manual(values=c("TRUE"=highlight_color, "FALSE"="dodgerblue3"), name="in window") +
      theme(legend.position = c(0.25, 0.2), legend.spacing.y = unit(0.1, "cm"), legend.key.height = unit(0.45, "cm")) +
      xlab("Nucleotide position") + ylab("Cumulative UDP count")
  } else {
    p.udp_dist = ggplot(udp.plot.df) +
      geom_segment(aes(x = deletion_start, xend = deletion_end, y = y, yend = y, color = sharing), size = segment_size) +
      simple_theme + coord_cartesian(xlim = c(view_start, view_end)) + scale_y_reverse() +
      xlab("Nucleotide position") + ylab("Cumulative UDP count")
    if (color_by == "sharing") {
      p.udp_dist = p.udp_dist + scale_color_manual(values = udp_colors) +
        theme(legend.position = c(0.25, 0.2), legend.spacing.y = unit(0.1, "cm"), legend.key.height = unit(0.45, "cm"))
    } else {
      p.udp_dist = p.udp_dist + scale_color_manual(values = udp_colors, guide=F)
    }
  }

  count.plot.df = tibble(x = 1:xmax)
  udp.char.matrix = str_split_fixed(udp.df$udp, "", n = nchar(udp.df$udp[1]))
  # isPositionDel = function(i) {
  #   sapply(udp.df$udp, FUN=function(x) substr(x, i, i) == '*')
  # }
  isPositionDel = function(i) {
    (udp.char.matrix[,i] == '*')
  }
  getDelCount = function(i) {
    sum(udp.char.matrix[,i] == '*')
  }
  getDelReadCount = function(i, udp.df) {
    sum(isPositionDel(i) * udp.df$num_reads)
  }

  count.plot.df$udp_count = sapply(count.plot.df$x, FUN=getDelCount)
  count.plot.df$read_count = sapply(count.plot.df$x, FUN=getDelReadCount, udp.df)

  count.plot.df$udp_count_in_window = 0
  count.plot.df$read_count_in_window = 0
  if (color_by == "window") {
    udp.in_window.df = udp.df %>% dplyr::filter(has_deletion_in_window)
    udp.char.matrix = str_split_fixed(udp.in_window.df$udp, "", n = nchar(udp.in_window.df$udp[1]))
    count.plot.df$udp_count_in_window = sapply(count.plot.df$x, FUN=getDelCount)
    count.plot.df$read_count_in_window = sapply(count.plot.df$x, FUN=getDelReadCount, udp.in_window.df)
  } else if (color_by == "sharing") {
    udp.unshared.df = udp.df %>% dplyr::filter(udp_sharing != "both")
    udp.char.matrix = str_split_fixed(udp.unshared.df$udp, "", n = nchar(udp.unshared.df$udp[1]))
    count.plot.df$udp_count_in_window = sapply(count.plot.df$x, FUN=getDelCount)
    count.plot.df$read_count_in_window = sapply(count.plot.df$x, FUN=getDelReadCount, udp.unshared.df)
  }

  p.udp_count = ggplot(count.plot.df, aes(x=x, y=udp_count)) +
    geom_bar(stat="identity", fill="dodgerblue2") +
    geom_bar(mapping = aes(x=x, y=udp_count_in_window), stat="identity", fill=highlight_color) +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("UDP count") +
    coord_cartesian(xlim = c(view_start, view_end))
  p.read_count = ggplot(count.plot.df, aes(x=x, y=read_count)) +
    geom_bar(stat="identity", fill="dodgerblue2") +
    geom_bar(mapping = aes(x=x, y=read_count_in_window), stat="identity", fill=highlight_color) +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("Read count") +
    coord_cartesian(xlim = c(view_start, view_end))

  if (!is.na(highlight_site)) {
    p.udp_dist = p.udp_dist + geom_vline(xintercept = highlight_site, color="darkgreen", alpha=0.5)
    p.udp_count = p.udp_count + geom_vline(xintercept = highlight_site, color="darkgreen", alpha=0.5)
    p.read_count = p.read_count + geom_vline(xintercept = highlight_site, color="darkgreen", alpha=0.5)
  }
  if (!is.na(cut_site)) {
    p.udp_dist = p.udp_dist + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
    p.udp_count = p.udp_count + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
    p.read_count = p.read_count + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  p.full = egg::ggarrange(p.read_count, p.udp_count, p.udp_dist, ncol=1, heights=c(1,1,4),
                          top = plot_title, draw = F)
  return(p.full)
}


#' Plots the deletion "pileup" profile separately for each replicate.
#'
#' @param del_result Result from a call to deletion_analysis.
#' @param viewing_window Window on either size of the CRISPR cut site to show in the plot.
#' @return Returns a ggplot object.
#'
#' @examples
#' \dontrun{
#' del_results = deletion_analysis(regions, replicates)
#' deletion_profile_plot(del_results)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
deletion_profile_plot = function(del_result,
                                 viewing_window = 40) {
  check_is_del_result(del_result)
  delprofile.udp.df = del_result$replicate_alleles
  delprofile.udp.df$count_udp = 1
  # Old code that would show only deletion reads in the window of interest
  if (FALSE) {
    delprofile.udp.df = replicate.udp.df %>% dplyr::mutate(count_udp = ifelse(has_deletion_in_window, 1, 0))
  }
  p.del_profile_pct = del_profile_internal(delprofile.udp.df,
                                           del_result$region$cut_site, del_result$region$highlight_site, viewing_window,
                                           plot_title = "Relative to all reads",
                                           show_average = F, show_replicates = T, ratioToWT = F)
  p.del_profile_pct = p.del_profile_pct + theme(axis.title.x=element_blank())
  p.del_profile_wtratio = del_profile_internal(delprofile.udp.df,
                                               del_result$region$cut_site, del_result$region$highlight_site, viewing_window,
                                               plot_title = "Relative to WT",
                                               show_average = F, show_replicates = T, ratioToWT = T)
  del_profile_title = sprintf("%s deletion profile", del_result$region$name)
  p.del_profile = egg::ggarrange(cowplot::ggdraw() + cowplot::draw_label(del_profile_title, fontface='bold'),
                                 p.del_profile_pct, p.del_profile_wtratio,
                                 ncol=1, heights=c(0.1, 0.5, 0.5), draw = F)
  return(p.del_profile)
}


del_profile_internal = function(replicate.udp.df, cut_site, highlight_site, viewing_window = 40, plot_title = NA, show_average = T, show_replicates = T, ratioToWT = F) {
  if (viewing_window < 1) {
    warning("udp_plot: viewing_window should be a positive integer.")
    return(NULL)
  }
  view_start = max(1, cut_site - viewing_window)
  view_end = cut_site + viewing_window

  if (!show_average & !show_replicates) {
    stop("getDeletionProfilePlot: One of show_average and show_replicates should be true.")
  }
  if (nrow(replicate.udp.df) == 0) {
    return(egg::ggarrange(text_plot("No alleles to plot."), top=plot_title, draw = F))
  }
  xmax = nchar(replicate.udp.df$udp[1])
  #udp.char.matrix = NULL
  #cur.df = NULL
  # Get the deletion profile for each replicate separately
  isPositionDel = function(udp.char.matrix, i) {
    (udp.char.matrix[,i] == '*')
  }
  getDelReadCount = function(cur.df, udp.char.matrix, i, filterUdps) {
    if (filterUdps) {
      # only count deletions in UDPs where count_udp == 1
      sum(isPositionDel(udp.char.matrix, i) * cur.df$num_reads * cur.df$count_udp)
    } else {
      sum(isPositionDel(udp.char.matrix, i) * cur.df$num_reads)
    }
  }
  getDelpctDataframe = function(unique_reps.df, filterUdps = T, ratioToWT = F) {
    counts.list = list()
    for (i in 1:nrow(unique_reps.df)) {
      curType = unique_reps.df[i,]$type
      curReplicate = unique_reps.df[i,]$replicate
      count.df = tibble(x=1:xmax, type = curType, replicate = curReplicate)
      cur.df = replicate.udp.df %>% dplyr::filter(replicate == curReplicate, type == curType)
      udp.char.matrix = str_split_fixed(cur.df$udp, "", n = nchar(cur.df$udp[1]))
      if (ratioToWT) {
        wtCount = sum(cur.df$is_wt_allele * cur.df$num_reads)
        count.df$del_pct = sapply(count.df$x, FUN=function(i) {getDelReadCount(cur.df, udp.char.matrix, i, filterUdps)}) / wtCount
      } else {
        count.df$del_pct = 100 * sapply(count.df$x, FUN=function(i) {getDelReadCount(cur.df, udp.char.matrix, i, filterUdps)}) / sum(cur.df$num_reads)
      }
      counts.list = c(counts.list, list(count.df))
    }
    # Combine deletion profiles for replicates
    bind_rows(counts.list)
  }

  unique_reps.df = unique(replicate.udp.df %>% dplyr::select(type, replicate))
  delpct.plot.df = getDelpctDataframe(unique_reps.df, filterUdps=T, ratioToWT=ratioToWT)
  delpct.plot.df$replicate = as.character(delpct.plot.df$replicate)
  delpct.plot.df$type = paste(delpct.plot.df$type, "replicate")

  if (show_average) {
    # Merge gDNA replicates together (and similarly for cDNA), and
    # determine deletion profiles
    # merged.udp.df = summarise(replicate.udp.df %>% group_by(udp, type),
    #                           num_reads = sum(num_reads)) %>%
    #   arrange(-num_reads)
    merged.udp.df = summarise(replicate.udp.df %>% group_by(type, udp),
                              num_reads = sum(num_reads), count_udp = first(count_udp), is_wt_allele = first(is_wt_allele)) %>%
      arrange(-num_reads)
    unique_reps.df = unique(merged.udp.df %>% ungroup() %>% dplyr::select(type))
    counts.list = list()
    for (i in 1:nrow(unique_reps.df)) {
      curType = unique_reps.df[i,]$type
      count.df = tibble(x=1:xmax, type = curType)
      cur.df = merged.udp.df %>% dplyr::filter(type == curType)
      udp.char.matrix = str_split_fixed(cur.df$udp, "", n = nchar(cur.df$udp[1]))
      if (ratioToWT) {
        wtCount = sum(cur.df$is_wt_allele * cur.df$num_reads)
        count.df$del_pct = sapply(count.df$x, FUN=function(i) {getDelReadCount(cur.df, udp.char.matrix, i, filterUdps=T)}) / wtCount
      } else {
        count.df$del_pct = 100 * sapply(count.df$x, FUN=function(i) {getDelReadCount(cur.df, udp.char.matrix, i, filterUdps=T)}) / sum(cur.df$num_reads)
      }
      counts.list = c(counts.list, list(count.df))
    }
    merged.delpct.plot.df = bind_rows(counts.list)
    merged.delpct.plot.df$type[merged.delpct.plot.df$type == "gDNA"] = "gDNA average"
    merged.delpct.plot.df$type[merged.delpct.plot.df$type == "cDNA"] = "cDNA average"
    merged.delpct.plot.df$replicate = merged.delpct.plot.df$type
  }

  if (any(!replicate.udp.df$count_udp)) {
    # Make a deletion profile for unfiltered UDPs
    unique_reps.df = unique(replicate.udp.df %>% dplyr::select(type, replicate))
    delpct.plot.unfiltered.df = getDelpctDataframe(unique_reps.df, filterUdps=F)
    delpct.plot.unfiltered.df$replicate = as.character(delpct.plot.unfiltered.df$replicate)
    delpct.plot.unfiltered.df$type = paste(delpct.plot.unfiltered.df$type, "rep unflt")
    delpct.plot.df = bind_rows(delpct.plot.df, delpct.plot.unfiltered.df)
  }

  alpha_values = c(`cDNA average`=0.4, `gDNA average`=0.4, `cDNA replicate`=0.9, `gDNA replicate`=0.9, `cDNA rep unflt`=0.9, `gDNA rep unflt`=0.9)
  color_values = c(`cDNA average`="darkorange", `gDNA average`="cyan3", `cDNA replicate`="red", `gDNA replicate`="blue", `cDNA rep unflt`="orange", `gDNA rep unflt`="turquoise4")
  size_values = c(`cDNA average`=1.8, `gDNA average`=1.8, `cDNA replicate`=0.3, `gDNA replicate`=0.3, `cDNA rep unflt`=0.3, `gDNA rep unflt`=0.3)
  if (show_average & show_replicates) {
    plot.df = bind_rows(delpct.plot.df, merged.delpct.plot.df)
  } else if (show_average) {
    plot.df = merged.delpct.plot.df
    size_values["cDNA average"] = size_values["gDNA average"] = 1.8
    alpha_values["cDNA average"] = alpha_values["gDNA average"] = 0.8
  } else {
    plot.df = delpct.plot.df
    size_values["cDNA replicate"] = size_values["gDNA replicate"] = 0.4
  }
  plot.df$type = factor(plot.df$type, levels = c("cDNA average", "gDNA average", "cDNA replicate", "gDNA replicate", "cDNA rep unflt", "gDNA rep unflt"))
  plot.df$lty_dash = plot.df$type %in% c("cDNA rep unflt", "gDNA rep unflt")
  simple_theme = theme_bw() + theme(panel.grid.minor = element_line(colour="grey90", size=0.1), panel.grid.major = element_line(colour="grey90", size=0.1))
  p.gDNA_cDNA = ggplot() +
    geom_line(aes(x=x, y=del_pct, color=type, size=type, alpha=type, linetype=lty_dash, group=paste(type, replicate)), data = plot.df) +
    simple_theme + xlab("Nucleotide position") + ylab("Deletion frequency (%)") +
    coord_cartesian(xlim = c(view_start, view_end)) +
    scale_color_manual(values = color_values) +
    scale_size_manual(values = size_values) +
    scale_alpha_manual(values = alpha_values) +
    scale_linetype_manual(values = c("TRUE"="dashed", "FALSE"="solid"), guide=F)

  if (!is.na(plot_title)) {
    p.gDNA_cDNA = p.gDNA_cDNA + ggtitle(plot_title)

    num_udps = length(unique(replicate.udp.df$udp))
    included_udps = length(unique(replicate.udp.df %>% filter(count_udp != 0) %>% .$udp))
    if (included_udps < num_udps) {
      p.gDNA_cDNA = p.gDNA_cDNA + labs(subtitle=sprintf("%d UDPs included out of %d", included_udps, num_udps))
    }
  }
  if (!is.na(highlight_site)) {
    p.gDNA_cDNA = p.gDNA_cDNA + geom_vline(xintercept = highlight_site, color="darkgreen", alpha=0.5)
  }
  if (!is.na(cut_site)) {
    p.gDNA_cDNA = p.gDNA_cDNA + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  if (ratioToWT) {
    p.gDNA_cDNA = p.gDNA_cDNA + ylab("Del:WT ratio")
  }
  return(p.gDNA_cDNA)
}

udp_to_binary = function(udp) {
  udp_array = strsplit(udp,"")[[1]]
  binary_array = rep(0, length(udp_array))
  binary_array[udp_array == '*'] = 1
  return(binary_array)
}


#' Returns a set of plots summarising the expected power to detect a significant
#' effect given various allele fraction and effect size combinations.
#'
#' @param del_result Result from a call to deletion_analysis.
#' @param allele_min_reads The minimum number of reads that a deletion allele must have across all replicates to be included.
#' @param WT_fraction If specified, then the model will assume this fraction of WT reads
#' @return Returns a list with three ggplot objects:
#' \itemize{
#'   \item cv_plot
#'   \item power_plot
#'   \item replicate_allocation_plot
#' }
#'
#' @examples
#' \dontrun{
#' del_results = deletion_analysis(regions, replicates)
#' power_plots(del_results)
#' }
#' @seealso \code{\link{deletion_analysis}}
#' @export
#'
power_plots = function(del_result,
                       allele_min_reads = 100,
                       WT_fraction = NA) {
  pwr = power_analysis(del_result, allele_min_reads, WT_fraction)

  # Make some plots that show the relationship between UDP read fraction and CV
  # Below I plot the UDP read % (fraction of total reads for each UDP) against CV
  # I've also plotted UDP read count against CV, but read fraction makes more sense
  # since it controls for the number of reads in each replicate. (For very low UDP
  # read counts it might make sense to use read count though.)
  p.udpfrac_cv = ggplot(pwr$udp_stats, aes(x=udp_frac_mean * 100, y=udp_frac_cv, col = type, shape = allele_type, size = allele_type, alpha = allele_type)) +
    geom_point() +
    scale_x_log10(breaks = c(0.1, 1, 10)) +
    theme_bw() + xlab("Mean UDP read %") + ylab("UDP coeff of variation") +
    stat_function(fun = function(x) 1/sqrt(x / 100 * sum(pwr$gDNA_stats$udp_mean)), color = "grey50", alpha = 0.5, size = 1) + # poisson rate of CV
    stat_function(fun = function(x) {fit_UDP_CV(pwr$coefs$cDNA, x / 100)}, color = "red", alpha = 0.5, size = 1.3) +
    stat_function(fun = function(x) {fit_UDP_CV(pwr$coefs$gDNA, x / 100)}, color = "blue", alpha = 0.5, size = 1.3) +
    scale_color_manual(values = c(cDNA = "red", gDNA = "blue")) +
    scale_shape_manual(values = c(HDR = 17, WT = 15, Del = 16), name = "Allele type") +
    scale_size_manual(values = c(HDR = 3, WT = 3, Del = 2), name = "Allele type") +
    scale_alpha_manual(values = c(HDR = 1, WT = 1, Del = 0.7), name = "Allele type")

  # Add a label showing the equations of the fit
  fit_label = sprintf("cDNA: %.3g + %.3g/sqrt(x)\ngDNA: %.3g + %.3g/sqrt(x)",
                      pwr$coefs$cDNA["a"], pwr$coefs$cDNA["b"],
                      pwr$coefs$gDNA["a"], pwr$coefs$gDNA["b"])
  p.udpfrac_cv = p.udpfrac_cv + annotate("text", Inf, Inf, label = fit_label, hjust = 1, vjust = 1, size = 2.8, color = "grey30")

  p.cv_density = ggplot(pwr$udp_stats, aes(x = udp_frac_cv, color = type)) + geom_density() +
    scale_color_manual(values = c(cDNA = "red", gDNA = "blue")) +
    theme_bw() + xlab("UDP frac CV") + ylab("Density")
  p.cv_plots = egg::ggarrange(p.udpfrac_cv, p.cv_density, ncol=1, heights = c(2,1),
                              top = sprintf("%s variance summary", del_result$region$name), draw = F)

  n_cDNA_rep = length(del_result$replicates %>% dplyr::filter(type == "cDNA") %>% .$replicate)
  n_gDNA_rep = length(del_result$replicates %>% dplyr::filter(type == "gDNA") %>% .$replicate)

  if (!is.na(WT_fraction)) {
    cvString = sprintf("fitted CV of %.3f in cDNA, %.3f in gDNA (%.1f%% WT)", pwr$cDNA_WT_CV, pwr$gDNA_WT_CV, WT_fraction*100)
  } else {
    cvString = sprintf("actual WT CV of %.3f in cDNA, %.3f in gDNA", pwr$cDNA_WT_CV, pwr$gDNA_WT_CV)
  }

  df = tibble(x=seq(0.001, 0.999, by=0.001))
  # ggplot(df, aes(x = x*100)) +
  #   scale_x_log10(breaks = c(0.1, 1, 10)) +
  #   stat_function(fun = function(x) sd_ratio_cDNA_gDNA(x/100, x/100, n_cDNA_rep, n_gDNA_rep), color = "blue") +
  #   theme_bw() + theme(panel.grid.major = element_line(size = 1, colour = "grey95")) +
  #   xlab("UDP read %") + ylab("SE of cDNA:gDNA ratio")
  power_plot_title = sprintf("%s power estimate", del_result$region$name)

  p.power = ggplot(df, aes(x = x*100)) +
    scale_x_log10(breaks = c(0.1, 1, 10)) +
    stat_function(fun = function(x) calc_power(pwr, 1.05, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "1.05x"), size=1) +
    stat_function(fun = function(x) calc_power(pwr, 1.1, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "1.1x"), size=1) +
    stat_function(fun = function(x) calc_power(pwr, 1.2, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "1.2x"), size=1) +
    stat_function(fun = function(x) calc_power(pwr, 1.5, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "1.5x"), size=1) +
    stat_function(fun = function(x) calc_power(pwr, 2, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "2.0x"), size=1) +
    theme_bw() + xlab("UDP read %") + ylab("Power") +
    scale_color_discrete(name = "Effect size") +
    ggtitle(label = power_plot_title,
            subtitle = sprintf("Based on %d cDNA and %d gDNA replicates, %s", n_cDNA_rep, n_gDNA_rep, cvString)) +
    theme(panel.grid.major = element_line(size = 1, colour = "grey95"),
          plot.title = element_text(size=14, hjust=0.5, face="bold"),
          plot.subtitle = element_text(size=10, hjust=0.5))

  # We also make a plot that shows the optimal allocation of replicates to cDNA vs. gDNA.
  get_power_df = function(power_res, nRep, udp_fraction) {
    df1.1 = tibble(num_cDNA = seq(1:(nRep-1)), effectSize = 1.1, udp_fraction = udp_fraction)
    df1.1$power = calc_power(power_res, df1.1$effectSize, df1.1$udp_fraction, df1.1$udp_fraction, df1.1$num_cDNA, nRep - df1.1$num_cDNA)
    df1.2 = tibble(num_cDNA = seq(1:(nRep-1)), effectSize = 1.2, udp_fraction = udp_fraction)
    df1.2$power = calc_power(power_res, df1.2$effectSize, df1.2$udp_fraction, df1.2$udp_fraction, df1.2$num_cDNA, nRep - df1.2$num_cDNA)
    df = rbind(df1.1, df1.2) %>%
      mutate(nReps = nRep, cDNA_rep_fraction = num_cDNA / nRep) %>%
      select(nReps, num_cDNA, cDNA_rep_fraction, everything())
    df$effectSize = sprintf("%.1fx", df$effectSize)
    df
  }
  replicate_allocation_plot = function(pwr, udp_fraction, showTitle=T) {
    df.fraction = tibble(cDNA_rep_fraction=seq(0.01, 0.99, by=0.01))
    df6 = get_power_df(pwr, 6, udp_fraction)
    df10 = get_power_df(pwr, 10, udp_fraction)
    df15 = get_power_df(pwr, 15, udp_fraction)
    df25 = get_power_df(pwr, 25, udp_fraction)
    nReps = factor(as.character(1:25), levels = as.character(1:25))
    #p = ggplot(df, aes(x = x)) +
    p = ggplot(df.fraction, aes(x = cDNA_rep_fraction)) +
      geom_point(data = df6,  mapping = aes(x = cDNA_rep_fraction, y = power, color = effectSize, shape = "6"), size = 2.5) +
      geom_point(data = df10, mapping = aes(x = cDNA_rep_fraction, y = power, color = effectSize, shape = "10"), size = 2.5) +
      geom_point(data = df15, mapping = aes(x = cDNA_rep_fraction, y = power, color = effectSize, shape = "15"), size = 2.5) +
      geom_point(data = df25, mapping = aes(x = cDNA_rep_fraction, y = power, color = effectSize, shape = "25"), size = 2.5) +
      stat_function(fun = function(cDNA_rep_fraction) calc_power(pwr, 1.1, udp_fraction, udp_fraction, 6*cDNA_rep_fraction, 6*(1-cDNA_rep_fraction)),   mapping = aes(color = "1.1x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(cDNA_rep_fraction) calc_power(pwr, 1.1, udp_fraction, udp_fraction, 10*cDNA_rep_fraction, 10*(1-cDNA_rep_fraction)), mapping = aes(color = "1.1x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(cDNA_rep_fraction) calc_power(pwr, 1.1, udp_fraction, udp_fraction, 15*cDNA_rep_fraction, 15*(1-cDNA_rep_fraction)), mapping = aes(color = "1.1x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(cDNA_rep_fraction) calc_power(pwr, 1.1, udp_fraction, udp_fraction, 25*cDNA_rep_fraction, 25*(1-cDNA_rep_fraction)), mapping = aes(color = "1.1x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(cDNA_rep_fraction) calc_power(pwr, 1.2, udp_fraction, udp_fraction, 6*cDNA_rep_fraction, 6*(1-cDNA_rep_fraction)),   mapping = aes(color = "1.2x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(cDNA_rep_fraction) calc_power(pwr, 1.2, udp_fraction, udp_fraction, 10*cDNA_rep_fraction, 10*(1-cDNA_rep_fraction)), mapping = aes(color = "1.2x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(cDNA_rep_fraction) calc_power(pwr, 1.2, udp_fraction, udp_fraction, 15*cDNA_rep_fraction, 15*(1-cDNA_rep_fraction)), mapping = aes(color = "1.2x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(cDNA_rep_fraction) calc_power(pwr, 1.2, udp_fraction, udp_fraction, 25*cDNA_rep_fraction, 25*(1-cDNA_rep_fraction)), mapping = aes(color = "1.2x"), size=0.8, alpha = 0.8) +
      theme_bw() + xlab("Fraction of replicates cDNA") + ylab("Power") +
      scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
      scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
      scale_color_discrete(name = "Effect size") +
      scale_shape_manual(name = "N Replicates", values = c("6"=16, "10"=17, "15"=15, "25"=8), breaks=c("6", "10", "15", "25")) +
      theme(panel.grid.major = element_line(size = 1, colour = "grey95"),
            plot.title = element_text(size=14, hjust=0.5, face="bold"),
            plot.subtitle = element_text(size=12, hjust=0.5))
    subtitleStr = sprintf("UDP fraction %d%%", as.integer(udp_fraction*100))
    if (showTitle) {
      p = p + ggtitle(label = power_plot_title, subtitle = subtitleStr)
    } else {
      p = p + ggtitle(label = NULL, subtitle = subtitleStr)
    }
    p
  }

  p.replicate_allocation = egg::ggarrange(replicate_allocation_plot(pwr, 0.01, showTitle=T),
                                          replicate_allocation_plot(pwr, 0.1, showTitle=F),
                                          ncol=1, heights = c(1,1), draw = F)

  replicate_allocation.df = expand.grid(udp_fraction = c(0.01, 0.1), nReps = c(6, 10, 15, 25)) %>%
    group_by(nReps, udp_fraction) %>%
    group_map(~ get_power_df(pwr, .y$nReps, .y$udp_fraction) %>% select(-nReps, -udp_fraction))

  plot_list = list(cv_plot = p.cv_plots, power_plot = p.power, replicate_allocation_plot = p.replicate_allocation)
  return(plot_list)
}

