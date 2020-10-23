
test_that("pedlist with empty markers give likelihood 1", {
  x = list(singleton("NN"), nuclearPed(1))
  x = setMarkers(x, locusAttributes = list(alleles = 1:2, name = "m1"))
  expect_equal(likelihood(x, "m1"), 1)
})

test_that("pedlist with HW markers give correct likelihood", {
  p = 0.9; q= 1-p
  x = nuclearPed(1)
  y = singleton("nn")
skip("likelihood-list tests need modifications")
  x = setMarkers(x, marker(x, '3'=1:2, alleles=1:2, afreq=c(p,q)))
  y = setMarkers(y, marker(y, nn=1, alleles=1:2, afreq=c(p,q)))
  liks = likelihood(list(x,y), 1, total=FALSE)
  expect_equal(liks, c(2*p*q, p^2))

  liktot = likelihood(list(x,y), list(mx, my))
  expect_equal(liktot, 2*p*q * p^2)
})

