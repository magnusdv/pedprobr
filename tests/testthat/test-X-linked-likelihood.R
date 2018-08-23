context("likelihoods on X")

test_that("empty X-markers give likelihood 1", {
  x1 = nuclearPed(1)
  expect_equal(likelihood(x1, marker(x1, chrom=23)), 1)
  expect_equal(likelihood(x1, marker(x1, alleles=1:2, chrom=23)), 1)
  expect_equal(likelihood(x1, marker(x1, alleles=letters[1:3], chrom=23)), 1)

  x2 = cousinPed(1)
  expect_equal(likelihood(x2, marker(x2, chrom=23)), 1)

  x3 = halfCousinPed(0)
  expect_equal(likelihood(x3, marker(x3, chrom=23)), 1)

})

test_that("nuclear X-likelihoods are correct, HW-like", {
  x = nuclearPed(1)
  p = 0.1; q = 1-p
  m1 = marker(x, '1'=1, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m1), p)

  m2 = marker(x, '2'=1:2, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m2), 2*p*q)

  m3 = marker(x, '3'=1, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m3), p)

  m4 = marker(x, '1'=1, '2'=2, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m4), p*q^2)

  m5 = marker(x, '2'=2, '3'=2, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m5), q^2)

  m6 = marker(x, '1'=1, '2'=1, '3'=2, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m6), 0)

  m7 = marker(x, '2'=1, '3'=2, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m7), 0)
})

test_that("nuclear X-likelihoods are correct, multiallelic markers", {
  x = nuclearPed(1)
  p = 0.1; q = 1-p

  m1 = marker(x, alleles=1:3, afreq=c(p,q,0), chrom=23)
  expect_equal(likelihood(x, m1), 1)

  m2 = marker(x, '1'=1, alleles=1:3, afreq=c(p,q,0), chrom=23)
  expect_equal(likelihood(x, m2), p)

  m3 = marker(x, '2'=1, alleles=1:3, afreq=c(p,q,0), chrom=23)
  expect_equal(likelihood(x, m3), p^2)

  m4 = marker(x, '1'=2, alleles=1:3, afreq=c(p,p,q-p), chrom=23)
  expect_equal(likelihood(x, m4), p)

  m5 = marker(x, '1'=1, '2'=10, alleles=1:10, chrom=23)
  expect_equal(likelihood(x, m5), .1^3)
})


test_that("likelihoods are correct in trio X-SNP examples", {
  x = nuclearPed(1, sex=2)
  p = 0.2; q = 1 - p

  m1 = marker(x, '1'=1, '2'=1:2, '3'=1:2, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m1), p * 2*p*q * .5)

  # Same, without dad
  m2 = marker(x, '2'=1:2, '3'=1:2, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m2), p*q)
})

test_that("likelihoods are correct in sibling X-SNP examples", {
  x = nuclearPed(2)
  p = 0.2; q = 1 - p

  # Brothers
  m1 = marker(x, '3'=1, '4'=1, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m1), p * .5*(1+p))

  # Sisters
  x = swapSex(x, 3:4)
  m2 = marker(x, '3'=1, '4'=1, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m2), p^2 * .5*(1+p))
})

test_that("male-male inheritance does not occur on X.", {
  x = addSon(nuclearPed(1),3)
  p = 0.2; q = 0.3; r = 1-p-q
  # Both homozygous
  m1 = marker(x, '1'=1, '3'=2, '5'=3, alleles=1:3, afreq=c(p,q,r), chrom=23)
  expect_equal(likelihood(x, m1), p*q*r)

  m2 = marker(x, '1'=3, '3'=3, '5'=3, alleles=1:3, afreq=c(p,q,r), chrom=23)
  expect_equal(likelihood(x, m2), r^3)
})

test_that("X-likelihoods are correct in 3 gener. SNP example", {
  x = nuclearPed(1)
  x = addParents(x, 1)
  x = addParents(x, 2)
  p = 0.2; q = 1 - p
  skip("")
  # Both homozygous
  m = marker(x, '4'=1:2, '7'=1:2, '3'=1:2, alleles=1:2, afreq=c(p,q), chrom=23)
  expect_equal(likelihood(x, m), (2*p*q)^2 * (3/8 + p*q/2))
})

test_that("X-likelihoods are correct in looped ped", {
  x = fullSibMating(1)
  p = 0.1; q = 0.2; r = 1-p-q

  m1 = marker(x, '5'=1:2, '6'=1:2, alleles=1:3, afreq=c(p,q,r), chrom=23)
  m2 = marker(x, '5'=1:2, '6'=1:2, alleles=1:4, afreq=c(p,q,r/2,r/2), chrom=23)
  expect_equal(likelihood(x, m1, verbose=F),
               likelihood(x, m2, verbose=F))

  m1 = marker(x, '5'=2, '6'=2, alleles=1:3, afreq=c(p,q,r), chrom=23)
  m2 = marker(x, '5'=2, '6'=2, alleles=4:1, afreq=c(r/2,r/2,q,p), chrom=23)
  expect_equal(likelihood(x, m1, verbose=F),
               likelihood(x, m2, verbose=F))
})

