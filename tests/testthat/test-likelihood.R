context("likelihoods")

test_that("empty markers give likelihood 1", {
  x = nuclearPed(1)
  expect_equal(likelihood(x, marker(x)), 1)
  expect_equal(likelihood(x, marker(x, alleles=1:2)), 1)
  expect_equal(likelihood(x, marker(x, alleles=letters[1:3])), 1)

  x2 = cousinPed(1)
  expect_equal(likelihood(x2, marker(x2)), 1)

  x3 = halfCousinPed(0)
  expect_equal(likelihood(x3, marker(x3)), 1)

  #x4 = cousinPed(0, child=T)
  #expect_equal(likelihood(x4, marker(x4)), 1)
})

test_that("nuclear likelihoods are correct, HW-like", {
  x = nuclearPed(1)
  p = 0.1; q = 1-p
  m1 = marker(x, '1'=1, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m1), p^2)

  m2 = marker(x, '1'=1:2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m2), 2*p*q)

  m3 = marker(x, '3'=1, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m3), p^2)

  m4 = marker(x, '3'=2:1, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m4), 2*p*q)

  m5 = marker(x, '1'=1, '2'=2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m5), p^2*q^2)

  m6 = marker(x, '1'=1, '2'=2, '3'=1:2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m6), p^2*q^2)

  m7 = marker(x, '1'=1, '2'=1, '3'=2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m7), 0)

  m8 = marker(x, '1'=1, '3'=2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m8), 0)
})

test_that("nuclear likelihoods are correct, multiallelic markers", {
  x = nuclearPed(1)
  p = 0.1; q = 1-p

  m1 = marker(x, alleles=1:3, afreq=c(p,q,0))
  expect_equal(likelihood(x, m1), 1)

  m2 = marker(x, '1'=1, alleles=1:3, afreq=c(p,q,0))
  expect_equal(likelihood(x, m2), p^2)

  m3 = marker(x, '1'=2, alleles=1:3, afreq=c(p,p,q-p))
  expect_equal(likelihood(x, m3), p^2)

  #m4 = marker(x, '1'=2, alleles=1:3, afreq=c(p,p,q-p))
  #expect_equal(likelihood(x, m4), p^2)

  m5 = marker(x, '1'=1, '2'=10, alleles=1:10)
  expect_equal(likelihood(x, m5), .1^4)
})

test_that("likelihoods are correct in trio SNP examples", {
  x = nuclearPed(1)
  p = 0.2; q = 1 - p

  # All heterozygous
  m1 = marker(x, '1'=1:2, '2'=1:2, '3'=1:2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m1), 2*p^2*q^2)

  # Same, without mom
  m2 = marker(x, '1'=1:2, '3'=1:2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m2), p*q)
})

test_that("likelihoods are correct in sibling SNP example", {
  x = nuclearPed(2)
  p = 0.2; q = 1 - p

  # Both homozygous
  m = marker(x, '3'=2, '4'=2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m), 1/4*q^2*(1+q)^2)
})

test_that("likelihoods are correct in 3 generat. SNP example", {
  x = nuclearPed(1)
  x = addParents(x, 1)
  x = addParents(x, 2)
  p = 0.2; q = 1 - p

  # Both homozygous
  m = marker(x, '4'=1:2, '7'=1:2, '3'=1:2, alleles=1:2, afreq=c(p,q))
  expect_equal(likelihood(x, m), (2*p*q)^2 * (3/8 + p*q/2))
})

test_that("likelihoods are correct in looped ped", {
  x = fullSibMating(1)
  p = 0.1; q = 0.2; r = 1-p-q

  m1 = marker(x, '5'=1:2, '6'=1:2, alleles=1:3, afreq=c(p,q,r))
  m2 = marker(x, '5'=1:2, '6'=1:2, alleles=1:4, afreq=c(p,q,r/2,r/2))
  expect_equal(likelihood(x, m1, verbose=F),
               likelihood(x, m2, verbose=F))

  m3 = marker(x, '5'=2, '6'=2, alleles=1:3, afreq=c(p,q,r))
  m4 = marker(x, '5'=2, '6'=2, alleles=4:1, afreq=c(r/2,r/2,q,p))
  expect_equal(likelihood(x, m3, verbose=F),
               likelihood(x, m4, verbose=F))
})


