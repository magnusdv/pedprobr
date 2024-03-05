
OMD = function(...) oneMarkerDistribution(..., verbose = F)

### Setup
p = 0.2; q = 1-p
g = c("1/1", "1/2", "2/2")
hw = c(p^2, 2*p*q, q^2)
PQ_ARRAY = array(c(p,q), dimnames = list(1:2))
HW_ARRAY = array(hw, dimnames = list(g))

test_that("oneMarkerDist works in empty autosomal examples", {

  ### Nuclear
  x = nuclearPed(1)
  m = marker(x, alleles=1:2, afreq = c(p,q))
  omd1 = OMD(x, partialmarker=m, ids=1)
  omd3 = OMD(x, partialmarker=m, ids=3)
  omd13 = OMD(x, partialmarker=m, ids=c(1,3))
  expect_equal(omd1, HW_ARRAY)
  expect_equal(omd3, HW_ARRAY)
  expect_equal(omd13, array(hw * c(p, .5*p, 0, q, .5, p, 0, .5*q, q),
                            dim = c(3,3), dimnames = list(g,g)))

  ### Singleton
  s = singleton("a")
  ms = marker(s, alleles=letters[1:2], afreq = c(p,q))
  omd = OMD(s, partialmarker=ms, ids="a")
  expect_equal(omd, setNames(HW_ARRAY, c("a/a","a/b","b/b")))
})

test_that("oneMarkerDist works in empty X-linked examples", {

  ### Nuclear
  x = nuclearPed(1, sex = 2)
  m = marker(x, alleles=1:2, afreq = c(p,q), chrom = 23)
  omd1 = OMD(x, partialmarker=m, ids=1)
  omd3 = OMD(x, partialmarker=m, ids=3)
  omd13 = OMD(x, partialmarker=m, ids=c(1,3))
  expect_equal(omd1, PQ_ARRAY)
  expect_equal(omd3, HW_ARRAY)
  expect_equal(omd13, array(c(p^2, 0, p*q, p*q, 0, q^2), dim = c(2,3),
                            dimnames = list(1:2, g)))

  ### Singleton
  s1 = singleton(1)
  s2 = singleton(2, sex=2)
  ms = marker(s1, alleles=letters[1:2], afreq = c(p,q), chrom = "X")
  omd.s1 = OMD(s1, partialmarker=ms, ids=1)
  omd.s2 = OMD(s2, partialmarker=ms, ids=2)
  expect_equal(omd.s1, setNames(PQ_ARRAY, letters[1:2]))
  expect_equal(omd.s2, setNames(HW_ARRAY, c("a/a","a/b","b/b")))
})

test_that("oneMarkerDist works in conditional nuclear example", {

  x = nuclearPed(1, sex = 2)

  m = marker(x, `3` = 2, alleles=1:2, afreq = c(p,q))
  mX = marker(x, `3` = 2, alleles=1:2, afreq = c(p,q), chrom=23)

  res1 = HW_ARRAY; res1[] = c(0, p, q)
  expect_equal(OMD(x, partialmarker=m, ids=1), res1)

  res2 = PQ_ARRAY; res2[] = c(0, 1)
  expect_equal(OMD(x, partialmarker=mX, ids=1), res2)
  expect_equal(OMD(x, partialmarker=mX, ids=2), res1)
})


test_that("oneMarkerDist on X is independent of male-male", {
  x = nuclearPed(1) |>
    addMarker("1" = 1, alleles = 1:2, afreq = c(p,q), chrom = "X")

  expect_equal(OMD(x, partialmarker = 1, ids=c(1,3))[1, ], c("1"=p, "2"=q))
})


