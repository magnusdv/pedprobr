
test_that("theta works in two-generation pedigrees", {
  x = linearPed(2)
  m = marker(x, `5` = "A/A", afreq = c(A = 0.3, B = 0.7))

  expect_equal(likelihood(x, m), likTheta(x, m, theta = 0))
  expect_equal(likelihood(x, m, theta = 0.01), 0.0921)

  y = ancestralPed(2)
  m2 = marker(y, `7` = "A/A", afreq = c(A = 0.3, B = 0.7))
  expect_equal(likelihood(y, m2), likTheta(y, m2, theta = 0))
  expect_equal(likelihood(y, m2, theta = 0.01), 0.0921)
})



test_that("theta-correction gives error if inbred founders", {
  x = singleton(1) |> addMarker(`1` = "1/1", alleles = 1:2) |>
    setFounderInbreeding(1, 0.5)

  expect_error(likelihood(x, m, theta = 0.1), "Theta correction cannot be used in pedigrees with inbred founders")
})


test_that("X-linked pedigree theta correction is rejected", {
  x = nuclearPed(1) |>
    addMarker(geno = c("1", "1/1", "1"), alleles = 1:2, chrom = "X")
  expect_error(likelihood(x, theta = 0.1), "X-linked")
})
