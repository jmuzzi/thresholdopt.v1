# thresholdopt 0.1.0

First release of the minimal package: the threshold-optimisation estimator
of Muzzi et al. (2025) behind a two-function interface.

* `threshold_overlap()` — two named numeric vectors in, and out: the common
  universe (printed as a percentage of each list), the optimal threshold,
  the observed intersection, the consensus identifiers, the full sweep and
  the publication's diagnostic figure.
* `threshold_overlap_deg()` — the same, for differential-expression tables;
  defaults match `DESeq2::results()`.
* The test suite pins the published values (threshold 3221, intersection
  1096, universe 16931) on the bundled `sf1_overexpression` data.

Inference about the selected overlap — whether the shared membership is
significant, and whether the shared response is directionally coherent — is
a separate problem, addressed in the development line of the package and not
part of this release.
