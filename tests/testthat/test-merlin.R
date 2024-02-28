if(!checkMerlin())
  skip("Merlin not installed")

lm = function(..., verbose = FALSE)
  likelihoodMerlin(..., verbose = verbose)

expect_signif = function(a, b)
  expect_equal(signif(a, 3), signif(b, 3))

test_that("MERLIN gives correct likelihood in selfing ped linked markers", {
  x = setMarkers(selfingPed(1), alleleMatrix = data.frame(m1 = c("1/2", "1/1"), m2 = c("1/2", "1/1")))
  p = q = 0.5
  r = 0.5
  expect_signif(lm(x), .5*p^2*q^2*(r^2 + (1-r)^2))

  # Complete linkage (requires that the markers are one the same chrom)
  chrom(x) = 1
  r = 0
  expect_signif(lm(x, rho = 0), .5*p^2*q^2*(r^2 + (1-r)^2))

  # Complete  linkage (requires that the markers are one the same chrom)
  r = 0.1
  expect_signif(lm(x, rho = 0.1), .5*p^2*q^2*(r^2 + (1-r)^2))
})


test_that("MERLIN agrees with likelihood2() with two linked markers", {
  x = linearPed(2)
  m = marker(x, geno = c("1/1", NA, "1/2", NA, "1/1"))
  x = setMarkers(x, list(m, m))

  rho = 0.1
  expect_signif(likelihood2(x, 1, 2, rho = rho),
               lm(x, 1:2, rho = rho, verbose = F))

  # With marker names and pre-set chroms
  y = x
  name(y) = c("m1", "m2")
  chrom(y) = 2
  expect_signif(likelihood2(y, 1, 2, rho = rho),
                lm(y, 1:2, rho = rho))

  # With pre-set map
  map = data.frame(CHROM = 1, MARKER = name(y), CM = 1:2)
  expect_signif(likelihood2(y, 1, 2, rho = haldane(cM = 1)),
                lm(y, 1:2, linkageMap = map))

})

test_that("MINX agrees with likelihood2() with two linked markers", {
  x = linearPed(2, sex = 2:1) |>
    addMarker(name = "m1", geno = c("1", NA, NA, "1/2", "1"), chrom = "X") |>
    addMarker(name = "m2", geno = c("1", NA, NA, "1/2", "1"), chrom = "X")

  # Unlinked
  expect_signif(likelihood(x, 1)^2, lm(x, 1:2, rho = 0.5))

  # With marker names and pre-set chroms
  rho = 0.25
  expect_signif(likelihood2(x, "m1", "m2", rho = rho),
                lm(x, c("m1", "m2"), rho = rho))

  # With pre-set map
  map = data.frame(CHROM = "X", MARKER = name(x), CM = 1:2)
  expect_signif(likelihood2(x, 1, 2, rho = haldane(cM = 1)),
                lm(x, 1:2, linkageMap = map))
})

test_that("MINX gives correct likelihood with two markers", {
  m = data.frame(m1 = c("1/1", "1/2", "1/1"), m2 = c("1/1", "1/2", "1/1"), row.names = c(1,4,5))
  x = setMarkers(linearPed(2, sex=2:1), alleleMatrix = m)
  chrom(x) = "X"
  # plot(x, marker= 1:2)
  expect_signif(prod(likelihood(x)), lm(x))
})


