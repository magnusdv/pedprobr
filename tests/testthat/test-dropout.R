test_that("likelihood handles dropout for singletons", {
  p = 0.1
  q = 1 - p
  d = 0.2
  afr = c("1" = p, "2" = q)

  het = singleton() |> addMarker(geno = "1/2", afreq = afr)
  hom = singleton() |> addMarker(geno = "1/1", afreq = afr)
  ped = nuclearPed() |> addMarker(`3` = "1/1", afreq = afr)

  expect_equal(likelihood(het, dropout = d), 2*p*q * (1 - d)^2)
  expect_equal(likelihood(hom, dropout = d), p^2 * (1 - d^2) + 2*p*q*d*(1 - d))
  expect_equal(likelihood(het, dropout = 0), 2*p*q)
  expect_error(likelihood(ped, dropout = d), "only implemented for singletons")
})
