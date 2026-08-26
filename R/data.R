#' SF-1 overexpression rankings in adrenocortical carcinoma cells
#'
#' Complete differential-expression rankings of two independent SF-1/NR5A1
#' overexpression experiments in H295R adrenocortical carcinoma cells, as
#' re-analysed in the source publication. These are the data behind the
#' published optimal threshold of 3221 features with 1096 in the
#' intersection, which the package examples and tests reproduce exactly.
#'
#' @format A named list of two data frames, each with columns Symbol (gene
#'   symbol) and logFC (log2 fold change of overexpression versus control):
#'   \describe{
#'     \item{doghman}{17948 genes, Doghman et al. (2007), the Lalli
#'       laboratory.}
#'     \item{ferraz}{17974 genes, Ferraz-de-Souza et al. (2011).}
#'   }
#'
#' @source Re-analysis reported in Muzzi JCD, et al., *European Journal of
#'   Endocrinology* 193(1):135-145, 2025. \doi{10.1093/ejendo/lvaf138}
#'
#' @examples
#' data(sf1_overexpression)
#' vapply(sf1_overexpression, nrow, integer(1))
"sf1_overexpression"
