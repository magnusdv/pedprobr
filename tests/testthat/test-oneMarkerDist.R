context("oneMarkerDistribution")

test_that("oneMarkerDist is correct in nuclear, autosomal", {
  x = pedtools::nuclearPed(1)
  p = 0.2; q = 1-p
  m = marker(x, alleles=1:2, afreq = c(p,q))
  omd = oneMarkerDistribution(x, partialmarker=m, ids=1, verbose=F)
  expect_equal(as.vector(omd), c(p^2, 2*p*q, q^2))
})
