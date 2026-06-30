test_that("count_theoretical_peptides handles basic trypsin cleavage", {
  # MASTK -> [MASTK] (no internal K/R before end)
  # PEPTIDER -> [PEPTIDER]
  # PEPTIDEK -> [PEPTIDEK]
  # MASTERPIECE: M-A-S-T-E-R | P-I-E-C-E  → "MASTER" (6 aa) + "PIECE" (5 aa)
  # With min_len=6, max_len=30: only "MASTER" qualifies
  seq_simple <- "MASTERPIECER"
  result <- count_theoretical_peptides(seq_simple, min_len = 6, max_len = 30)
  expect_type(result, "integer")
  expect_gte(result, 0L)
})

test_that("count_theoretical_peptides returns NA for empty/NA input", {
  expect_identical(count_theoretical_peptides(""),  NA_integer_)
  expect_identical(count_theoretical_peptides(NA_character_), NA_integer_)
})

test_that("count_theoretical_peptides respects length limits", {
  # Sequence: AAAAAKBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBK (too long C-term frag)
  # AAAAAK (6 aa) qualifies; second fragment too long
  seq <- paste0("AAAAAK", paste(rep("B", 35), collapse = ""), "K")
  n_strict <- count_theoretical_peptides(seq, min_len = 6, max_len = 30)
  expect_equal(n_strict, 1L)

  # With relaxed max_len both fragments qualify
  n_relaxed <- count_theoretical_peptides(seq, min_len = 6, max_len = 50)
  expect_equal(n_relaxed, 2L)
})

test_that("count_theoretical_peptides handles no-K/R sequence", {
  # No K or R -> single fragment
  seq <- "AAAAAAAAAA"  # 10 aa, no trypsin sites
  result <- count_theoretical_peptides(seq, min_len = 6, max_len = 30)
  expect_equal(result, 1L)
})

test_that("count_theoretical_peptides respects proline rule", {
  # KP should NOT be cut (proline rule)
  # MASTKPEPTIDE: no cut between K and P -> one fragment "MASTKPEPTIDE" (12 aa)
  seq_kp <- "MASTKPEPTIDE"
  result <- count_theoretical_peptides(seq_kp, min_len = 6, max_len = 30)
  expect_equal(result, 1L)

  # MASTKEPTIDE: K followed by E -> should cut -> "MASTK" (5 aa) + "EPTIDE" (6 aa)
  seq_ke <- "MASTKEPTIDE"
  result_ke <- count_theoretical_peptides(seq_ke, min_len = 6, max_len = 30)
  expect_equal(result_ke, 1L)  # only "EPTIDE" passes min_len=6; "MASTK"=5 fails
})

test_that("count_theoretical_peptides handles missed cleavages", {
  # PEPTIDEK|AAAAAK -> fragments: "PEPTIDEK" (8aa), "AAAAAK" (6aa)
  # mc=0: both qualify -> 2
  # mc=1: also "PEPTIDEKAAAAK" -> 3
  seq <- "PEPTIDEKAAAAAK"
  n0 <- count_theoretical_peptides(seq, min_len = 6, max_len = 30, max_missed = 0)
  n1 <- count_theoretical_peptides(seq, min_len = 6, max_len = 30, max_missed = 1)
  expect_equal(n0, 2L)
  expect_equal(n1, 3L)
})
