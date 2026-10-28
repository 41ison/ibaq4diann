## Helpers: a mock report with all four DIA-NN identifier columns, and a
## matching FASTA. Two protein groups (P002, P003) share the gene HBA1.
write_mock_inputs <- function() {
  report <- tibble::tibble(
    Protein.Group = c("P001", "P001", "P002", "P003", "P004;cRAP-P004", "P005"),
    Protein.Ids = c("P001", "P001", "P002", "P003", "cRAP-P004;P004", "P005"),
    Protein.Names = c(
      "ALBU_HUMAN",
      "ALBU_HUMAN",
      "HBA_HUMAN",
      "HBA3_HUMAN",
      "CATA_HUMAN;cRAP-CATA_HUMAN",
      "NOGN_HUMAN"
    ),
    Genes = c("ALB", "ALB", "HBA1", "HBA1", "CAT;cRAP-CAT", ""),
    Run = c("S1", "S2", "S1", "S1", "S1", "S1"),
    Stripped.Sequence = c("PEPA", "PEPA", "PEPB", "PEPC", "PEPD", "PEPE"),
    Precursor.Normalised = c(100, 200, 300, 400, 500, 600),
    Q.Value = 0.001,
    PG.Q.Value = 0.001,
    Lib.PG.Q.Value = 0.001,
    Proteotypic = 1L,
    PG.MaxLFQ.Quality = 0.9
  )
  parquet <- tempfile(fileext = ".parquet")
  arrow::write_parquet(report, parquet)

  fasta <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      # 2 tryptic peptides >= 6 aa each
      ">sp|P001|ALBU_HUMAN Albumin OS=Homo sapiens GN=ALB PE=1 SV=1",
      "AAAAAAKBBBBBBK",
      ">sp|P002|HBA_HUMAN Hb alpha OS=Homo sapiens GN=HBA1 PE=1 SV=1",
      "AAAAAAK",
      ">sp|P003|HBA3_HUMAN Hb alpha 3 OS=Homo sapiens GN=HBA1 PE=1 SV=1",
      "AAAAAAKBBBBBBKCCCCCCK",
      ">sp|P004|CATA_HUMAN Catalase OS=Homo sapiens GN=CAT PE=1 SV=1",
      "AAAAAAKBBBBBBKCCCCCCKDDDDDDK",
      ">sp|P005|NOGN_HUMAN No gene OS=Homo sapiens PE=1 SV=1",
      "AAAAAAK"
    ),
    fasta
  )
  list(parquet = parquet, fasta = fasta)
}

run_quiet <- function(...) {
  suppressMessages(suppressWarnings(compute_ibaq(..., log2_transform = FALSE)))
}

test_that("compute_ibaq with Protein.Group matches accessions and leading members", {
  inp <- write_mock_inputs()
  out <- run_quiet(inp$parquet, inp$fasta, id_col = "Protein.Group")

  expect_equal(
    out$protein_id,
    c("P001", "P002", "P003", "P004;cRAP-P004", "P005")
  )
  expect_equal(out$n_theoretical_peptides, c(2L, 1L, 3L, 4L, 1L))
  expect_equal(out$S1, c(100 / 2, 300, 400 / 3, 500 / 4, 600))
  expect_equal(out$S2, c(200 / 2, NA, NA, NA, NA))
  expect_named(
    out,
    c(
      "protein_id",
      "Protein.Names",
      "Genes",
      "n_theoretical_peptides",
      "S1",
      "S2"
    )
  )
})

test_that("compute_ibaq with Protein.Ids finds a non-leading member in the FASTA", {
  inp <- write_mock_inputs()
  out <- run_quiet(inp$parquet, inp$fasta, id_col = "Protein.Ids")

  expect_true("cRAP-P004;P004" %in% out$protein_id)
  expect_equal(
    out$n_theoretical_peptides[out$protein_id == "cRAP-P004;P004"],
    4L
  )
  expect_true("Protein.Group" %in% names(out))
})

test_that("compute_ibaq with Protein.Names matches FASTA entry names", {
  inp <- write_mock_inputs()
  out <- run_quiet(inp$parquet, inp$fasta, id_col = "Protein.Names")

  expect_equal(nrow(out), 5L)
  expect_equal(
    out$n_theoretical_peptides[out$protein_id == "ALBU_HUMAN"],
    2L
  )
  expect_false("Protein.Names" %in% names(out))
  expect_true(all(c("Protein.Group", "Genes") %in% names(out)))
})

test_that("compute_ibaq with Genes merges groups, drops empty genes", {
  inp <- write_mock_inputs()
  expect_message(
    out <- suppressWarnings(
      compute_ibaq(
        inp$parquet,
        inp$fasta,
        id_col = "Genes",
        log2_transform = FALSE
      )
    ),
    "empty 'Genes'"
  )

  # P005 (empty gene) dropped; P002 + P003 merged into HBA1
  expect_equal(sort(out$protein_id), c("ALB", "CAT;cRAP-CAT", "HBA1"))

  hba1 <- out[out$protein_id == "HBA1", ]
  # First FASTA entry for HBA1 is P002 (1 peptide); intensities summed
  expect_equal(hba1$n_theoretical_peptides, 1L)
  expect_equal(hba1$S1, 300 + 400)
  expect_equal(hba1$Protein.Group, "P002;P003")
  expect_equal(hba1$Protein.Names, "HBA_HUMAN;HBA3_HUMAN")
})

test_that("compute_ibaq rejects unsupported id_col", {
  inp <- write_mock_inputs()
  expect_error(compute_ibaq(inp$parquet, inp$fasta, id_col = "Run"))
})

test_that("compute_ibaq errors when nothing matches the FASTA", {
  inp <- write_mock_inputs()
  other_fasta <- tempfile(fileext = ".fasta")
  writeLines(c(">sp|X1|X_MOUSE x OS=Mus GN=X", "AAAAAAK"), other_fasta)
  expect_error(
    suppressMessages(compute_ibaq(inp$parquet, other_fasta)),
    "None of the"
  )
})
