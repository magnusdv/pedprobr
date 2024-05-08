
test_that("likelihood() catches input errors", {
  x = singleton(1)
  y = nuclearPed(1)
  m = marker(x)
  expect_error(likelihood(x, y), "Invalid input for argument `markers`")
  expect_error(likelihood(x, y, m), "Invalid input for argument `markers`")
  expect_error(likelihood(x, y, marker = m), "Invalid input for argument `peelOrder`")

  expect_error(likelihood(x, marker = list(1)), "Invalid input for argument `markers`.")
  expect_error(likelihood(x, marker = NA), "`markers` is a logical of length 1")
  expect_error(likelihood(x, marker = 1), "Marker index out of range: 1")
  expect_error(likelihood(x, marker = "a"), "Unknown marker name: a")

  peds = list(x,y)
  expect_error(likelihood(peds, m), fixed = T,
               "`likelihood.list()` requires `markers` to be a vector of marker names or indices.")
  expect_error(likelihood(peds, list()), fixed = T,
               "`likelihood.list()` requires `markers` to be a vector of marker names or indices.")
  expect_error(likelihood(peds, 1), "Marker index out of range: 1")

  peds2 = list(setMarkers(x,marker(x)), y)
  expect_error(likelihood(peds2), "The pedigree components have different number of markers attached")
  expect_error(likelihood(peds2,1), "Marker index out of range: 1")
})


test_that("empty markers give likelihood 1", {
  x = nuclearPed(1)
  expect_equal(likelihood(x, marker(x)), 1)
  expect_equal(likelihood(x, marker(x, alleles=1:2)), 1)
  expect_equal(likelihood(x, marker(x, alleles=letters[1:3])), 1)

  x2 = cousinPed(1)
  expect_equal(likelihood(x2, marker(x2)), 1)

  x3 = halfCousinPed(0)
  expect_equal(likelihood(x3, marker(x3)), 1)

  x4 = cousinPed(0, child=T)
  expect_equal(likelihood(x4, marker(x4)), 1)
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

  m4 = marker(x, '1'=2, alleles=1:3, afreq=c(p,p,q-p))
  expect_equal(likelihood(x, m4), p^2)

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

test_that("likelihoods is the same with different peelorder", {
  x = avuncularPed() |> addMarker(geno = c(NA,NA, "1/2", "1/1", NA,NA))
  xx = reorderPed(x, 6:1)
  expect_equal(likelihood(x), likelihood(xx))
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

test_that("likelihood is correct with partial genotypes", {
  p = 0.8; q = 1-p

  # Singleton
  x = singleton("a") |> addMarker(a = "2/-", alleles = 1:2, afreq= c(p, q))
  expect_equal(likelihood(x), q * (q + 2*p))

  # Nuclear
  y = nuclearPed(ch = "a") |> addMarker(a = "2/-", alleles = 1:2, afreq= c(p, q))
  expect_equal(likelihood(y), q * (q + 2*p))
})

# test_that("pedprobr and Familias give identical results", {
#   skip("")
#   library(forrel)
#   library(pedlikCompare)
#   testRandom = function() {
#     x = randomPed(8)
#     ids = sample(labels(x), size = sample(pedsize(x), 1))
#     nals = sample.int(5, size = 1) + 1
#     afr = runif(nals, .2, .8)
#     afr = afr/sum(afr)
#     sims = markerSim(x, N=10, ids = ids, alleles = sample.int(10, size = nals),
#                      afreq = afr, mutmod = "prop", rate = 0.1, verbose = FALSE)
#     sims = reorderPed(sims, 8:1)
#     for(i in 1:10) {
#       tes = compare(sims, i, verbose = FALSE, programs = c("pedprobr", "Familias"))
#       if(nrow(tes) < 2) break
#       if(!isTRUE(all.equal(tes[1,2], tes[2,2]))) {
#         message("MISMATCH!")
#         print(selectMarkers(x, i))
#       }
#     }
#   }
#   for(i in 1:100) testRandom()
# })
