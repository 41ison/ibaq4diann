## Helper: build a minimal mock DIA-NN report tibble
make_mock_report <- function() {
  tibble::tibble(
    Protein.Group = c("P001", "P001", "P002", "P003"),
    Run = c("S1", "S1", "S1", "S1"),
    Stripped.Sequence = c("PEPTIDEK", "SAMPLEK", "UNIQUEK", "ANOTHER"),
    Precursor.Normalised = c(1e6, 2e6, 3e6, 4e6),
    Q.Value = c(0.001, 0.005, 0.02, 0.001), # P002 > cutoff
    PG.Q.Value = c(0.001, 0.005, 0.001, 0.001),
    Lib.PG.Q.Value = c(0.001, 0.005, 0.001, 0.001),
    Proteotypic = c(1L, 1L, 1L, 0L), # P003 non-proteotypic
    PG.MaxLFQ.Quality = c(0.9, 0.8, 0.9, 0.9)
  )
}

test_that("filter_diann_report removes rows above Q.Value cutoff", {
  report <- make_mock_report()
  out <- filter_diann_report(report, q_value_cutoff = 0.01)
  expect_true(all(out$Q.Value <= 0.01))
  expect_equal(nrow(out), 3L)
})

test_that("filter_diann_report removes non-proteotypic rows", {
  report <- make_mock_report()
  out <- filter_diann_report(report, q_value_cutoff = NULL, proteotypic_only = TRUE)
  expect_true(all(out$Proteotypic == 1))
  expect_equal(nrow(out), 3L)
})

test_that("filter_diann_report skips missing columns gracefully", {
  # Remove Q.Value column – should not error
  report <- make_mock_report()
  report$Q.Value <- NULL
  expect_no_error(
    filter_diann_report(report, q_value_cutoff = 0.01)
  )
})

test_that("filter_diann_report NULL cutoffs skip corresponding filters", {
  report <- make_mock_report()
  out <- filter_diann_report(
    report,
    q_value_cutoff = NULL,
    pg_q_value_cutoff = NULL,
    lib_pg_q_value_cutoff = NULL,
    proteotypic_only = FALSE,
    lfq_quality_cutoff = NULL
  )
  # All rows should pass
  expect_equal(nrow(out), nrow(report))
})

test_that("filter_diann_report applies LFQ quality filter", {
  report <- make_mock_report()
  report$PG.MaxLFQ.Quality[1] <- 0.3 # below 0.5 threshold
  out <- filter_diann_report(
    report,
    q_value_cutoff = NULL,
    pg_q_value_cutoff = NULL,
    lib_pg_q_value_cutoff = NULL,
    proteotypic_only = FALSE,
    lfq_quality_cutoff = 0.5
  )
  expect_false(any(out$PG.MaxLFQ.Quality < 0.5))
})
