## Helper: write a small UniProt-style FASTA to a temp file
write_mock_fasta <- function() {
  path <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">sp|P001|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2",
      "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFK",
      ">sp|P002|HBA_HUMAN Hemoglobin alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2",
      "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK",
      ">tr|P003|P003_HUMAN Isoform of HBA1 OS=Homo sapiens OX=9606 GN=HBA1 PE=3 SV=1",
      "MVLSPADKTNVKAAWGK",
      ">sp|P004|NOGN_HUMAN No gene field OS=Homo sapiens OX=9606 PE=1 SV=1",
      "MKKKRRRPEPTIDEK",
      ">PLAIN001 a plain header",
      "MPEPTIDEKPEPTIDER"
    ),
    path
  )
  path
}

test_that("parse_fasta keys sequences by accession by default", {
  seqs <- suppressMessages(parse_fasta(write_mock_fasta()))
  expect_named(seqs, c("P001", "P002", "P003", "P004", "PLAIN001"))
})

test_that("parse_fasta keys sequences by entry name", {
  seqs <- suppressMessages(parse_fasta(
    write_mock_fasta(),
    id_type = "entry_name"
  ))
  expect_named(
    seqs,
    c("ALBU_HUMAN", "HBA_HUMAN", "P003_HUMAN", "NOGN_HUMAN", "PLAIN001")
  )
})

test_that("parse_fasta keys by gene, dropping entries without GN= and duplicates", {
  expect_message(
    expect_message(
      seqs <- parse_fasta(write_mock_fasta(), id_type = "gene"),
      "no 'gene' identifier"
    ),
    "keeping the first occurrence"
  )
  # P003 shares HBA1 with P002 -> first occurrence (P002's sequence) kept
  expect_named(seqs, c("ALB", "HBA1"))
  expect_equal(
    unname(seqs[["HBA1"]]),
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK"
  )
})

test_that("parse_fasta rejects unknown id_type", {
  expect_error(parse_fasta(write_mock_fasta(), id_type = "foo"))
})

test_that("match_theoretical_peptides uses the first member found in the FASTA", {
  theo <- tibble::tibble(
    protein_id = c("P001", "P002"),
    n_theoretical_peptides = c(10L, 5L)
  )
  out <- match_theoretical_peptides(
    c("P001", "cRAP-P002;P002", "P003;P001", "P999"),
    theo
  )
  expect_equal(out$fasta_id, c("P001", "P002", "P001", NA))
  expect_equal(out$n_theoretical_peptides, c(10L, 5L, 10L, NA))
})
