# thresholdopt

<!-- badges: start -->
[![R-CMD-check](https://github.com/jmuzzi/thresholdopt.v1/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jmuzzi/thresholdopt.v1/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
<!-- badges: end -->

Estimation of the cut-off at which two independently ranked gene lists agree
most strongly.

Comparisons between differential-expression studies are confounded by the
cut-off each study applied: an adjusted p-value threshold, a fold-change
floor or a fixed top-*k* are not comparable across experiments, and the size
of the reported overlap follows the choice rather than the biology.
`thresholdopt` estimates that cut-off from the two rankings themselves and
returns the consensus feature set at the estimated optimum, with the
diagnostic figure of the source publication.

The method is the one published in

> Muzzi JCD, Ruggiero C, Doghman-Bouguerra M, Colodel ME, Magno JM, Resende
> JSS, Durand N, de Moura JF, Alvarenga LM, Cavalli LR, Figueiredo BC,
> Lalli E, Castro MAA. **Steroidogenic factor-1 regulates a core set of
> target genes to promote malignancy in adrenocortical carcinoma.**
> *European Journal of Endocrinology* 193(1):135–145, 2025.
> [doi:10.1093/ejendo/lvaf138](https://doi.org/10.1093/ejendo/lvaf138)

## Installation

```r
# install.packages("remotes")
remotes::install_github("jmuzzi/thresholdopt.v1", build_vignettes = TRUE)
```

`build_vignettes = TRUE` is what makes `vignette("thresholdopt")` available;
without it `remotes` installs the code but skips the vignette, and the
vignette is where the worked example lives. It requires `knitr` and
`rmarkdown`.

Only base R is required to run the functions; `ggplot2`, `cowplot` and
`scales` are needed for the figure.

## Usage

Two named numeric vectors — names are gene identifiers, values are the
ranking statistic (only `abs(value)` is used for ranking):

```r
library(thresholdopt)
data(sf1_overexpression)

lalli  <- setNames(sf1_overexpression$doghman$logFC,
                   sf1_overexpression$doghman$Symbol)
ferraz <- setNames(sf1_overexpression$ferraz$logFC,
                   sf1_overexpression$ferraz$Symbol)

fit <- threshold_overlap(lalli, ferraz)
#> Common universe      : 16931 genes (97.1% of lalli, 94.2% of ferraz)
#> Optimal threshold    : 3221
#> Observed intersection: 1096
```

The printed threshold and intersection are the values reported in the source
publication, and the figure is its diagnostic plot: the two ranked lists
with the selected prefix highlighted, and the objective function with the
intersection size on a secondary axis.

For differential-expression tables (defaults match `DESeq2::results()`:
identifiers in row names, statistic in `log2FoldChange`, `NA` rows dropped):

```r
fit <- threshold_overlap_deg(res_treatment, res_validation)
```

## How it works

Every candidate set size `n` is evaluated: the intersection `q(n)` of the
two top-`n` sets is scored with a hypergeometric tail probability computed
in log space (no underflow plateau at genome scale) and Bonferroni-adjusted
over the sweep; the score is regularised by the divergence
`d(n) = n − q(n)` through `λ(n) = d(n)/max d(n)`, and the threshold
minimises `λ(n) · log10 p_adj(n)`. The regularisation is what makes the
optimum well defined: the raw score degenerates for identical lists and for
reversed lists, exactly where the divergence vanishes. The whole sweep is
`O(N log N)`.

The optimiser always returns a threshold, including for unrelated lists —
the objective is a selection score, not a calibrated test of the selected
result. Verdicts about the overlap (convergent, divergent, inconclusive)
are the subject of the development version of this package.

See `vignette("thresholdopt")` for a worked example — or read it online in
[`vignettes/thresholdopt.Rmd`](vignettes/thresholdopt.Rmd) without
installing anything.

## Getting help

Bug reports and questions: <https://github.com/jmuzzi/thresholdopt.v1/issues>,
ideally with a minimal reproducible example and `sessionInfo()`.

## Citation

Please cite the paper that introduced the method (above);
`citation("thresholdopt")` gives the BibTeX entry.

## License

MIT © João C. D. Muzzi. See [LICENSE.md](LICENSE.md).
