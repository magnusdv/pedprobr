context("inbred founders")

test_that("likelihood of inbred singleton agrees with expanded pedigree, and theory", {
  p = 0.1; q = 1-p
  x = addChildren(nuclearPed(2, sex=1:2), 3, 4, 1)
  mhom = marker(x, '5'=1, alleles=1:2, afreq=c(p,q), name="hom")
  mhet = marker(x, '5'=1:2, alleles=1:2, afreq=c(p,q), name="het")
  x = setMarkers(x, list(mhom, mhet))

  y = subset(x, 5)
  founder_inbreeding(y, 5) = 0.25

  lik_y_hom = likelihood(y, "hom")
  lik_y_het = likelihood(y, "het")

  # agree with expanded
  expect_equal(lik_y_hom, likelihood(x, "hom", verbose = F))
  expect_equal(lik_y_het, likelihood(x, "het", verbose = F))

  # agree with theory
  expect_equal(lik_y_hom, 0.25 * p + 0.75 * p^2)
  expect_equal(lik_y_het, 0.25 * 0 + 0.75 * 2*p*q)
})


test_that("likelihood of nuc with inbred founder is correct", {
  p = 0.1; q = 1-p
  x = addSon(addChildren(nuclearPed(2, sex=1:2), 3, 4, 1), 5)
  mhom = marker(x, '5'=1, alleles=1:2, afreq=c(p,q), name="hom")
  mhet = marker(x, '5'=1:2, alleles=1:2, afreq=c(p,q), name="het")
  x = setMarkers(x, list(mhom, mhet))

  y = branch(x, 5)
  founder_inbreeding(y, 5) = 0.25

  lik_y_hom = likelihood(y, "hom")
  lik_y_het = likelihood(y, "het")

  # agree with expanded
  expect_equal(lik_y_hom, likelihood(x, "hom", verbose = F))
  expect_equal(lik_y_het, likelihood(x, "het", verbose = F))

  # agree with theory
  expect_equal(lik_y_hom, 0.25 * p + 0.75 * p^2)
  expect_equal(lik_y_het, 0.25 * 0 + 0.75 * 2*p*q)
})

