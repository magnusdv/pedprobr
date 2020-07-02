context("likelihoods of linked markers")

liktest = function(x, m1, m2) {
  th_0.0 = likelihood(x, m1, m2, rho=0, verbose=F)
  th_0.25 = likelihood(x, m1, m2, rho=0.25, verbose=F)
  th_0.5 = likelihood(x, m1, m2, rho=0.5, verbose=F)
  c(th_0.0, th_0.25, th_0.5)
}

test_that("linked empty markers give likelihood 1", {
  x1 = nuclearPed(1)
  m1.1 = marker(x1)
  expect_equal(liktest(x1, m1.1, m1.1), c(1,1,1))

  m1.2 = marker(x1, alleles=1:2)
  expect_equal(liktest(x1, m1.2, m1.2), c(1,1,1))

  x2 = cousinPed(1)
  m2 = marker(x2, alleles=1:3)
  expect_equal(liktest(x2, m2, m2), c(1,1,1))

  x3 = halfCousinPed(0, child=T)
  m3 = marker(x3, alleles=1:3)
  expect_equal(liktest(x3, m3, m3), c(1,1,1))
})

test_that("two linked HW-like markers are indep of rho", {
  p = 0.9; q = 1-p
  r = 0.8; s = 1-r

  x = nuclearPed(1)
  m1 = marker(x, `1`=1:2, alleles=1:2, afreq=c(p,q))
  m2 = marker(x, `1`=1:2, alleles=1:2, afreq=c(r,s))
  answ = 2*p*q * 2*r*s
  expect_equal(liktest(x, m1, m2), rep(answ,3))

  m3 = marker(x, `1`=1, `2`=1:2, alleles=1:2, afreq=c(p,q))
  m4 = marker(x, `1`=2, `2`=1:2, alleles=1:2, afreq=c(r,s))
  answ = p^2*2*p*q * s^2*2*r*s
  expect_equal(liktest(x, m3, m4), rep(answ,3))

  y = addSon(x, 3)
  m5 = marker(y, `5`=1:2, alleles=1:2, afreq=c(p,q))
  m6 = marker(y, `5`=1:2, alleles=1:2, afreq=c(r,s))
  answ = 2*p*q * 2*r*s
  expect_equal(liktest(y, m5, m6), rep(answ,3))
})

test_that("two X-linked HW-like markers are indep of rho", {
  p = 0.9; q = 1-p
  r = 0.8; s = 1-r

  # Male
  x = nuclearPed(1)
  m1 = marker(x, `1`=1, alleles=1:2, afreq=c(p,q), chrom=23)
  m2 = marker(x, `1`=2, alleles=1:2, afreq=c(r,s), chrom=23)
  answ = p*s
  expect_equal(liktest(x, m1, m2), rep(answ,3))

  # Female
  m3 = marker(x, `2`=1, alleles=1:2, afreq=c(p,q), chrom=23)
  m4 = marker(x, `2`=2, alleles=1:2, afreq=c(r,s), chrom=23)
  answ = p^2 * s^2
  expect_equal(liktest(x, m3, m4), rep(answ,3))

  y = addDaughter(x, 3)
  m5 = marker(y, `5`=1:2, alleles=1:2, afreq=c(p,q), chrom=23)
  m6 = marker(y, `5`=1:2, alleles=1:2, afreq=c(r,s), chrom=23)
  answ = 2*p*q * 2*r*s
  expect_equal(liktest(y, m5, m6), rep(answ,3))
})

