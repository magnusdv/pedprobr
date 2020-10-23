
test_that("empty marker (with mutations) give likelihood 1", {
  x = nuclearPed(1)
  m = marker(x, alleles=1:3, mutmod = "random")
  expect_equal(likelihood(x, m), 1)
})

test_that("nuclear - foundersOnly - mutation", {
  x = nuclearPed(1)
  p = 0.1; q = 1-p

  m1 = marker(x, '1'=1, alleles=1:2, afreq=c(p,q), mutmod="ra")
  expect_equal(likelihood(x, m1), p^2)

  m2 = marker(x, '1'=1:2, alleles=1:2, afreq=c(p,q), mutmod="ra")
  expect_equal(likelihood(x, m2), 2*p*q)

  m3 = marker(x, '1'=1, '2'=1:2, alleles=1:3, afreq=c(p,q,0), mutmod="ra")
  expect_equal(likelihood(x, m3), p^2*2*p*q)

})

test_that("trio - someTyped - mutation", {
  x = nuclearPed(1)
  p = 0.9; q = 1-p

  m1 = marker(x, '1'=1, '2'=1, '3'=2,  afreq=c('1'=p,'2'=q), mutmod="eq", rate=0.1)
  rr = mutmod(m1)$male[1,2]
  expect_equal(likelihood(x, m1), p^4 * rr^2)

  m2 = marker(x, '1'=1, '3'=2,  afreq=c('1'=p,'2'=q), mutmod="pr", rate=0.1)
  rr = mutmod(m2)$male[1,2]
  expect_equal(likelihood(x, m2), p^2 * q * rr)
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

test_that("loops - mutation", {
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

test_that("ancestral ped - non-stationary mut", {
  mutmat = matrix(c(0,0,1,1), ncol=2, dimnames = list(1:2, 1:2))
  mut = pedmut::mutationModel("custom", matrix = mutmat)

  x = relabel(addParents(linearPed(2), 4, verbose = F), 1:7)
  m = marker(x, '3'=2, '6'=2, alleles=1:2, mutmod=mut)
  # plot(x,m)
  expect_equal(likelihood(x,m), 1)
})

test_that("lumped mutation model is robust to allele ordering", {
  x = nuclearPed(2)
  m1 = m2 = marker(x, alleles = 1:4, afreq = c(.1,.2,.3,.4),
                   mutmod = "prop", rate = 0.5)
  genotype(m1, '3') = 3:4
  genotype(m2, '3') = 4:3
  expect_identical(likelihood(x, m1), likelihood(x, m2))

  genotype(m1, '3') = 3
  genotype(m1, '4') = 4

  genotype(m2, '3') = 4
  genotype(m2, '4') = 3
  expect_identical(likelihood(x, m1), likelihood(x, m2))
})
