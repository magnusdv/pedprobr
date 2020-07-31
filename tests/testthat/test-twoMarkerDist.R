context("twoMarkerDistribution")

TMD = function(x, pm1, pm2, ...)
  twoMarkerDistribution(x, partialmarker1 = pm1, partialmarker2 = pm2, ..., verbose = F)

### Setup
p = 0.2; q = 1-p
g = c("1/1", "1/2", "2/2")
hw = c(p^2, 2*p*q, q^2)
PQ_ARRAY = array(c(p,q), dimnames = list(1:2))
HW_ARRAY = array(hw, dimnames = list(g))

test_that("twoMarkerDist works in empty autosomal examples", {

  ### Nuclear
  x = nuclearPed(1)
  m = marker(x, alleles=1:2, afreq = c(p,q))

  expect_equal(TMD(x, pm1=m, pm2=m, id=1, rho=0),
               HW_ARRAY %o% HW_ARRAY)
  expect_equal(TMD(x, pm1=m, pm2=m, id=3, rho=0.25),
               HW_ARRAY %o% HW_ARRAY)

  ### Singleton
  s = singleton("a")
  ms1 = marker(s, alleles=1:2, afreq = c(p,q))
  ms2 = marker(s, alleles=letters[1:2], afreq = c(p,q))
  HW_ab = setNames(HW_ARRAY, c("a/a","a/b","b/b"))

  expect_equal(TMD(s, pm1=ms1, pm2=ms2, id="a", rho=0),
               HW_ARRAY %o% HW_ab)
})

test_that("twoMarkerDist works in empty X-linked examples", {

  ### Nuclear
  x = nuclearPed(1, sex = 2)
  m = marker(x, alleles=1:2, afreq = c(p,q), chrom = 23)

  expect_equal(TMD(x, pm1=m, pm2=m, id=1, rho=0),
               PQ_ARRAY %o% PQ_ARRAY)
  expect_equal(TMD(x, pm1=m, pm2=m, id=3, rho=0.25),
               HW_ARRAY %o% HW_ARRAY)

  ### Singleton
  s1 = singleton(1)
  s2 = swapSex(s1, 1)
  ms1 = marker(s1, alleles=1:2, afreq = c(p,q), chrom = 23)
  ms2 = marker(s1, alleles=letters[1:2], afreq = c(p,q), chrom = 23)
  HW_ab = setNames(HW_ARRAY, c("a/a","a/b","b/b"))
  PQ_ab = setNames(PQ_ARRAY, c("a","b"))
  expect_equal(TMD(s1, pm1=ms1, pm2=ms2, id=1, rho=0),
               PQ_ARRAY %o% PQ_ab)
  expect_equal(TMD(s2, pm1=ms1, pm2=ms2, id=1, rho=0),
               HW_ARRAY %o% HW_ab)
})

test_that("twoMarkerDist works in conditional nuclear example", {

  x = nuclearPed(1, sex = 2)

  # Autosomal
  m1 = m2 = marker(x, `3` = 2, alleles=1:2, afreq = c(p,q))
  h = HW_ARRAY; h[] = c(0, p, q); res1 = h %o% h
  expect_equal(TMD(x, pm1=m1, pm2=m2, id=1, rho=0), res1)

  # X-linked
  mX1 = mX2 = marker(x, `3` = 2, alleles=1:2, afreq = c(p,q), chrom=23)
  res2 = PQ_ARRAY %o% PQ_ARRAY; res2[] = c(0,0,0,1)
  expect_equal(TMD(x, pm1=mX1, pm2=mX2, id=1, rho=0), res2)
  expect_equal(TMD(x, pm1=mX1, pm2=mX2, id=2, rho=0), res1)
})

test_that("recombination rate is recovered", {

  x = linearPed(2, sex = 2:1)
  rho = 0.15

  # Autosomal
  m1 = marker(x, `1`=1, `2`=2, `4`=1, `5`=1, alleles=1:2, afreq = c(p,q))
  m2 = marker(x, `1`=1, `2`=2, `4`=1, `5`=0, alleles=1:2, afreq = c(p,q))
  # plot(x, list(m1,m2))
  expect_equal(TMD(x, pm1=m1, pm2=m2, id=5, rho=rho)["1/1", "1/2"], rho)

  # X-linked
  mX1 = marker(x, `1`=1, `2`=2, `5`=1, alleles=1:2, afreq = c(p,q), chrom=23)
  mX2 = marker(x, `1`=1, `2`=2, `5`=0, alleles=1:2, afreq = c(p,q), chrom=23)
  # plot(x, list(mX1,mX2))
  expect_equal(TMD(x, pm1=mX1, pm2=mX2, id=5, rho=rho)["1", "2"], rho)

})

