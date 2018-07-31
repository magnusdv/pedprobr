context("likelihood of pedlist")

test_that("pedlist with empty markers give likelihood 1", {
  x = nuclearPed(1)
  y = singleton("nn")
  liks = likelihood(list(x,y), list(marker(x), marker(y)), total=FALSE)
  liks = likelihood(list(y,x), list(marker(y), marker(x)), total=FALSE)
  expect_equal(liks, c(1,1))
})

test_that("pedlist with HW markers give correct likelihood", {
  p = 0.9; q= 1-p
  x = nuclearPed(1)
  y = singleton("nn")

  mx = marker(x, '3'=1:2, alleles=1:2, afreq=c(p,q))
  my = marker(y, nn=1, alleles=1:2, afreq=c(p,q))
  liks = likelihood(list(x,y), list(mx,my), total=FALSE)
  expect_equal(liks, c(2*p*q, p^2))

  liktot = likelihood(list(x,y), list(mx, my))
  expect_equal(liktot, 2*p*q * p^2)
})

