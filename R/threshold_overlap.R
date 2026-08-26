#' Optimal threshold for the overlap of two ranked gene lists
#'
#' Implements the threshold-optimisation method of Muzzi et al. (2025):
#' instead of imposing an arbitrary cut-off (an adjusted p-value, a
#' fold-change floor, a top-k), the cut-off is estimated from the two
#' rankings themselves. Every candidate set size `n` is evaluated: the two
#' lists are cut at their top-`n`, the intersection `q(n)` is scored with a
#' hypergeometric tail probability (evaluated in log space, so genome-scale
#' lists do not saturate at the double-precision underflow limit), and that
#' score is regularised by the divergence `d(n) = n - q(n)` through
#' `lambda(n) = d(n) / max d(n)`. The selected threshold minimises the
#' objective `lambda(n) * log10 p_adj(n)`, which is degenerate at neither
#' extreme of the sweep.
#'
#' The function prints the common universe (also as a percentage of each
#' input list), the optimal threshold and the observed intersection, and
#' draws the diagnostic figure of the source publication: the two ranked
#' lists on top with the selected prefix highlighted, and the objective
#' function below with the intersection size on a secondary axis.
#'
#' @param x,y Named numeric vectors: names are feature identifiers (for
#'   example gene symbols) and values are the ranking statistic (for example
#'   a log2 fold change). The input order is irrelevant; features are ranked
#'   by decreasing absolute value. Plain character vectors already in rank
#'   order are also accepted.
#' @param label.x,label.y Names of the two experiments, used in the printed
#'   summary and inside the figure. Default: the names of the objects passed
#'   by the caller, as in [plot()].
#'
#' @return (Invisibly) a list with `threshold` (the optimal set size),
#'   `intersection` (features shared by the two top-`threshold` sets),
#'   `universe` (size of the common universe), `consensus` (the shared
#'   feature identifiers at the optimum), `curve` (a data frame with the
#'   whole sweep: `n`, `intersection`, `log10_p`, `log10_p_adjusted`,
#'   `divergence`, `objective`) and `plot` (the ggplot object).
#'
#' @details
#' The regularisation is what makes the optimum well defined: the raw
#' criterion alone returns `p = 0` at every size when a list is compared
#' with itself, and is minimised only at the extremes of the sweep when a
#' list is compared with its own reverse. The divergence vanishes exactly
#' where the raw criterion degenerates. The tail probability is Bonferroni
#' adjusted over the sweep, in log space, before the regularisation, as in
#' the source publication.
#'
#' The whole sweep costs O(N log N): the intersection curve is a cumulative
#' count over each feature's larger rank.
#'
#' @references
#' Muzzi JCD, Ruggiero C, Doghman-Bouguerra M, et al. Steroidogenic factor-1
#' regulates a core set of target genes to promote malignancy in
#' adrenocortical carcinoma. *European Journal of Endocrinology*
#' 193(1):135-145, 2025. \doi{10.1093/ejendo/lvaf138}
#'
#' @examples
#' data(sf1_overexpression)
#' lalli  <- setNames(sf1_overexpression$doghman$logFC,
#'                    sf1_overexpression$doghman$Symbol)
#' ferraz <- setNames(sf1_overexpression$ferraz$logFC,
#'                    sf1_overexpression$ferraz$Symbol)
#' fit <- threshold_overlap(lalli, ferraz)
#' fit$threshold      # 3221, as published
#'
#' @importFrom utils globalVariables
#' @export
threshold_overlap <- function(x, y, label.x = NULL, label.y = NULL) {
  if (is.null(label.x)) label.x <- paste(deparse(substitute(x)), collapse = "")
  if (is.null(label.y)) label.y <- paste(deparse(substitute(y)), collapse = "")

  xi <- .as_ranked(x, "x"); yi <- .as_ranked(y, "y")

  common <- intersect(xi$id, yi$id)
  N <- length(common)
  if (N < 3L) stop("The common universe has fewer than 3 features.",
                   call. = FALSE)
  ## each list is filtered and ranked in ITS OWN original row order, so tie
  ## breaking matches the source implementation; ranks are aligned afterwards
  xi <- xi[xi$id %in% common, ]; yi <- yi[yi$id %in% common, ]
  r1 <- rank(-abs(xi$stat), ties.method = "first")
  r2 <- rank(-abs(yi$stat), ties.method = "first")
  common <- xi$id
  r2 <- r2[match(common, yi$id)]
  yi <- yi[match(common, yi$id), ]

  ## the whole intersection curve in one cumulative count
  mx <- pmax(r1, r2)
  n <- seq_len(N - 1L)
  q <- cumsum(tabulate(mx, nbins = N))[n]

  lp <- stats::phyper(q, m = n, n = N - n, k = n,
                      lower.tail = FALSE, log.p = TRUE) / log(10)
  lp[!is.finite(lp)] <- -322
  lp_adj <- pmin(lp + log10(length(n)), 0)       # Bonferroni, in log space
  d <- n - q
  lambda <- d / max(d)
  objective <- lambda * lp_adj

  best <- which.min(objective)
  threshold <- n[best]
  consensus <- common[mx <= threshold]

  cat(sprintf("Common universe      : %d genes (%.1f%% of %s, %.1f%% of %s)\n",
              N, 100 * N / .n_in(x), label.x, 100 * N / .n_in(y), label.y))
  cat(sprintf("Optimal threshold    : %d\n", threshold))
  cat(sprintf("Observed intersection: %d\n", q[best]))

  curve <- data.frame(n = n, intersection = q, log10_p = lp,
                      log10_p_adjusted = lp_adj, divergence = d,
                      objective = objective)
  p <- .paper_figure(curve, threshold, q[best],
                     data.frame(rank = r1, stat = xi$stat),
                     data.frame(rank = r2, stat = yi$stat),
                     c(label.x, label.y))
  print(p)

  invisible(list(threshold = threshold, intersection = q[best],
                 universe = N, consensus = consensus, curve = curve,
                 plot = p))
}

#' Optimal threshold for two differential-expression tables
#'
#' Convenience wrapper around [threshold_overlap()] for the data frames
#' produced by differential-expression pipelines. Defaults match
#' `DESeq2::results()`: feature identifiers in the row names and the ranking
#' statistic in the `log2FoldChange` column; rows with a missing statistic
#' are dropped, and duplicated identifiers keep their first occurrence.
#'
#' @param df1,df2 Data frames (or objects coercible to data frames, such as
#'   the `DESeqResults` object of `DESeq2::results()`).
#' @param id_col Name of the identifier column. `NULL` (default) uses the
#'   row names, which is where `DESeq2` keeps the gene identifiers.
#' @param stat_col Name of the ranking-statistic column. Default
#'   `"log2FoldChange"`.
#' @inheritParams threshold_overlap
#'
#' @return The same (invisible) list as [threshold_overlap()].
#'
#' @examples
#' data(sf1_overexpression)
#' lalli  <- data.frame(log2FoldChange = sf1_overexpression$doghman$logFC,
#'                      row.names = sf1_overexpression$doghman$Symbol)
#' ferraz <- data.frame(log2FoldChange = sf1_overexpression$ferraz$logFC,
#'                      row.names = sf1_overexpression$ferraz$Symbol)
#' fit <- threshold_overlap_deg(lalli, ferraz)
#'
#' @export
threshold_overlap_deg <- function(df1, df2, id_col = NULL,
                                  stat_col = "log2FoldChange",
                                  label.x = NULL, label.y = NULL) {
  if (is.null(label.x)) label.x <- paste(deparse(substitute(df1)), collapse = "")
  if (is.null(label.y)) label.y <- paste(deparse(substitute(df2)), collapse = "")
  to_vector <- function(d, nm) {
    d <- as.data.frame(d)
    ids <- if (is.null(id_col)) rownames(d) else as.character(d[[id_col]])
    if (!stat_col %in% names(d))
      stop("column '", stat_col, "' not found in ", nm,
           " (DESeq2 results use 'log2FoldChange').", call. = FALSE)
    v <- stats::setNames(d[[stat_col]], ids)
    v <- v[!is.na(v) & nzchar(names(v))]
    v[!duplicated(names(v))]
  }
  threshold_overlap(to_vector(df1, "df1"), to_vector(df2, "df2"),
                    label.x = label.x, label.y = label.y)
}

## ---- internals -------------------------------------------------------------

# ggplot2 non-standard evaluation: these are data-frame columns, not globals
if (getRversion() >= "2.15.1")
  utils::globalVariables(c("n", "objective", "stat", "sel"))


.as_ranked <- function(v, nm) {
  if (is.character(v)) {
    return(data.frame(id = v, stat = rev(seq_along(v)),
                      stringsAsFactors = FALSE))
  }
  if (is.null(names(v)) || !is.numeric(v))
    stop("'", nm, "' must be a named numeric vector (names = identifiers) ",
         "or a character vector in rank order.", call. = FALSE)
  v <- v[!is.na(v) & nzchar(names(v))]
  v <- v[!duplicated(names(v))]
  data.frame(id = names(v), stat = unname(v), stringsAsFactors = FALSE)
}

.n_in <- function(v) if (is.character(v)) length(v) else
  length(v[!is.na(v) & nzchar(names(v))])

.paper_figure <- function(curve, threshold, intersection, d1, d2, labels) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("cowplot", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE))
    stop("Packages 'ggplot2', 'cowplot' and 'scales' are required for the ",
         "figure.", call. = FALSE)

  lim1 <- range(curve$objective, finite = TRUE)
  lim2 <- range(curve$intersection)
  to_left  <- function(v) (v - lim2[1]) / diff(lim2) * diff(lim1) + lim1[1]
  to_right <- function(v) (v - lim1[1]) / diff(lim1) * diff(lim2) + lim2[1]

  g_bottom <- ggplot2::ggplot(curve, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = objective)) +
    ggplot2::geom_line(ggplot2::aes(y = to_left(intersection)),
                       colour = "red") +
    ggplot2::annotate("point", x = threshold,
                      y = curve$objective[curve$n == threshold],
                      colour = "green") +
    ggplot2::scale_y_continuous(
      name = "Objective function",
      sec.axis = ggplot2::sec_axis(
        ~ to_right(.), name = "Genes in the intersection (x1000)",
        labels = scales::number_format(scale = 1 / 1000))) +
    ggplot2::labs(x = "Gene set size") +
    ggplot2::theme(
      axis.title.y.right = ggplot2::element_text(colour = "red"),
      axis.text.y.right  = ggplot2::element_text(colour = "red"),
      axis.ticks.y.right = ggplot2::element_line(colour = "red"),
      plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5, "pt"))

  rank_panel <- function(d, lab, last) {
    dd <- data.frame(n = d$rank, stat = abs(d$stat))
    dd$sel <- ifelse(dd$n <= threshold, "Yes", "No")
    ggplot2::ggplot(dd, ggplot2::aes(x = n)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = stat, fill = sel),
                           show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = c(No = "grey50", Yes = "green")) +
      ggplot2::annotate("text", x = mean(dd$n), y = 1, label = lab,
                        hjust = 0.5, vjust = 0.5) +
      ggplot2::labs(
        x = if (last) "Position in the ranked list of genes" else NULL,
        y = "Absolute log2 FC") +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(if (last) 0 else 5.5, 5.5,
                                      if (last) 5.5 else 0, 5.5, "pt"))
  }

  cowplot::plot_grid(rank_panel(d1, labels[1], FALSE),
                     rank_panel(d2, labels[2], TRUE),
                     g_bottom, ncol = 1, rel_heights = c(0.6, 0.6, 1.6),
                     align = "v", axis = "lr")
}
