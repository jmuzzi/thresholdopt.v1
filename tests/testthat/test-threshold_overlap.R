test_that("the published optimum is reproduced exactly", {
  data(sf1_overexpression)
  lalli  <- setNames(sf1_overexpression$doghman$logFC,
                     sf1_overexpression$doghman$Symbol)
  ferraz <- setNames(sf1_overexpression$ferraz$logFC,
                     sf1_overexpression$ferraz$Symbol)
  pdf(NULL)
  fit <- suppressMessages(threshold_overlap(lalli, ferraz))
  dev.off()
  expect_equal(fit$threshold, 3221L)
  expect_equal(fit$intersection, 1096L)
  expect_equal(fit$universe, 16931L)
  expect_length(fit$consensus, 1096L)
})

test_that("the DESeq2-style wrapper gives the identical result", {
  data(sf1_overexpression)
  d1 <- data.frame(log2FoldChange = sf1_overexpression$doghman$logFC,
                   row.names = sf1_overexpression$doghman$Symbol)
  d2 <- data.frame(log2FoldChange = sf1_overexpression$ferraz$logFC,
                   row.names = sf1_overexpression$ferraz$Symbol)
  pdf(NULL)
  fit <- suppressMessages(threshold_overlap_deg(d1, d2))
  dev.off()
  expect_equal(fit$threshold, 3221L)
  expect_equal(fit$intersection, 1096L)
})

test_that("NA statistics and duplicated identifiers are handled", {
  data(sf1_overexpression)
  d1 <- data.frame(log2FoldChange = sf1_overexpression$doghman$logFC,
                   row.names = sf1_overexpression$doghman$Symbol)
  d2 <- data.frame(log2FoldChange = sf1_overexpression$ferraz$logFC,
                   row.names = sf1_overexpression$ferraz$Symbol)
  d1$log2FoldChange[1:5] <- NA          # as DESeq2 results contain
  pdf(NULL)
  expect_no_error(fit <- suppressMessages(threshold_overlap_deg(d1, d2)))
  dev.off()
  expect_true(fit$universe <= 16931L)
})

test_that("labels default to the caller's object names", {
  data(sf1_overexpression)
  meu_exp1 <- setNames(sf1_overexpression$doghman$logFC,
                       sf1_overexpression$doghman$Symbol)
  meu_exp2 <- setNames(sf1_overexpression$ferraz$logFC,
                       sf1_overexpression$ferraz$Symbol)
  pdf(NULL)
  out <- capture.output(threshold_overlap(meu_exp1, meu_exp2))
  dev.off()
  expect_true(any(grepl("meu_exp1", out)))
  expect_true(any(grepl("meu_exp2", out)))
})

test_that("degenerate input is refused", {
  expect_error(threshold_overlap(c(a = 1), c(a = 1)), "fewer than 3")
  expect_error(threshold_overlap(1:5, 1:5), "named numeric")
})
