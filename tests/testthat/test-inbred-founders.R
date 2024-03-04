
test_that("likelihood of inbred singleton agrees with expanded pedigree, and theory", {
  p = 0.1; q = 1-p
  x = nuclearPed(2, sex=1:2) |> addSon(3:4)
  mhom = marker(x, '5'=1, alleles=1:2, afreq=c(p,q), name="hom")
  mhet = marker(x, '5'=1:2, alleles=1:2, afreq=c(p,q), name="het")
  x = setMarkers(x, list(mhom, mhet))

  y = subset(x, 5)
  founderInbreeding(y, 5) = 0.25

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
  founderInbreeding(y, 5) = 0.25

  lik_y_hom = likelihood(y, "hom")
  lik_y_het = likelihood(y, "het")

  # agree with expanded
  expect_equal(lik_y_hom, likelihood(x, "hom", verbose = F))
  expect_equal(lik_y_het, likelihood(x, "het", verbose = F))

  # agree with theory
  expect_equal(lik_y_hom, 0.25 * p + 0.75 * p^2)
  expect_equal(lik_y_het, 0.25 * 0 + 0.75 * 2*p*q)
})

test_that("founder inb raises error in likelihood of linked markers", {
  x = nuclearPed(1)
  m = marker(x, "1" = 1, alleles = 1:2)
  founderInbreeding(x, 1) = 1

  # ped method
  expect_error(likelihood2(x, m, m, rho = 0.1),
               "Likelihood of linked markers is not implemented in pedigrees with founder inbreeding.")

  # singleton method
  s = subset(x, 1)
  expect_error(likelihood2(s, m, m, rho = 0.1),
               "Likelihood of linked markers is not implemented in pedigrees with founder inbreeding.")
})

test_that("complete inbreeding + heterozygosity = 0", {
  s = singleton(1) |> addMarker(geno = "1/2")
  founderInbreeding(s, 1) = 1
  expect_equal(likelihood(s), 0)

  x = nuclearPed(1) |> addMarker("1" = "1/2")
  founderInbreeding(x, 1) = 1
  expect_equal(likelihood(x), 0)

  y = list(x,s)
  expect_equal(likelihood(y), 0)

})

