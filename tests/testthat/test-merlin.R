if(Sys.which("merlin.exe") == "")
  skip("Merlin not installed")

lm = function(..., verbose = FALSE)
  likelihoodMerlin(..., verbose = verbose)

test_that("MERLIN gives correct likelihood in selfing ped linked markers", {
  x = setMarkers(selfingPed(1), alleleMatrix = data.frame(m1 = c("1/2", "1/1"), m2 = c("1/2", "1/1")))
  p = q = 0.5
  r = 0.5
  expect_equal(lm(x), .5*p^2*q^2*(r^2 + (1-r)^2), tolerance = 0.0001)

  # Complete linkage (requires that the markers are one the same chrom)
  chrom(x) = 1
  r = 0
  expect_equal(lm(x, rho = 0), .5*p^2*q^2*(r^2 + (1-r)^2), tolerance = 0.0001)

  # Complete  linkage (requires that the markers are one the same chrom)
  r = 0.1
  expect_equal(lm(x, rho = 0.1), .5*p^2*q^2*(r^2 + (1-r)^2), tolerance = 0.0001)
})


test_that("MINX gives correct likelihood with two markers", {
  m = data.frame(m1 = c("1/1", "1/2", "1/1"), m2 = c("1/1", "1/2", "1/1"), row.names = c(1,4,5))
  x = setMarkers(linearPed(2, sex=2:1), alleleMatrix = m)
  chrom(x) = "X"
  # plot(x, 1:2)
  expect_equal(prod(likelihood(x)), lm(x), tolerance = 0.00001)
})


