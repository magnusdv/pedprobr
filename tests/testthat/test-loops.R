

test_that("several copies of one loop breaker give correct likelihoods", {
  p = c(0.2, 0.3)

  x = ped(
    id = 1:6,
    fid = c(0, 0, 1, 1, 1, 5),
    mid = c(0, 0, 2, 2, 3, 3),
    sex = c(1, 2, 2, 2, 1, 1)
  ) |>
    addMarker(`6` = "1/1", alleles = 1:2, afreq = c(p[1], 1 - p[1])) |>
    addMarker(`6` = "1/1", alleles = 1:2, afreq = c(p[2], 1 - p[2]))

  repeated = cbind(loopBreaker = c("3", "3"), child = c("5", "6"))
  alternative = cbind(loopBreaker = c("3", "5"), child = c("5", "6"))

  #repeated = c(3,3)
  #alternative = c(3,5)

  xr = breakLoops(x, repeated, verbose = FALSE)
  xa = breakLoops(x, alternative, verbose = FALSE)

  expected = p^2 + 3/8 * p * (1 - p)

  expect_equal(likelihood(xr), expected)
  expect_equal(likelihood(xa), expected)

  rho = c(0, 0.2, 0.5)
  lr = vapply(rho, \(r) likelihood2(xr, 1, 2, rho = r), 1)
  la = vapply(rho, \(r) likelihood2(xa, 1, 2, rho = r), 1)

  expect_equal(lr, la)
  expect_equal(lr[3], prod(expected))
})


test_that("founder loop breakers work for autosomal markers", {
  p = c(0.2, 0.3)

  x = ped(
    id = 1:6,
    fid = c(0, 0, 0, 1, 1, 4),
    mid = c(0, 0, 0, 2, 3, 5),
    sex = c(1, 2, 2, 1, 2, 1)
  ) |>
    addMarker(`6` = "1/1", alleles = 1:2, afreq = c(p[1], 1 - p[1])) |>
    addMarker(`6` = "1/2", alleles = 1:2, afreq = c(p[2], 1 - p[2]))

  founder = cbind(loopBreaker = "1", child = "5")
  alternative = cbind(loopBreaker = "5", child = "6")

  xf = breakLoops(x, founder, verbose = FALSE)
  xa = breakLoops(x, alternative, verbose = FALSE)

  expected = c(p[1]^2 + 1/8 * p[1] * (1 - p[1]),
               2 * p[2] * (1 - p[2]) * 7/8)

  expect_equal(likelihood(xf), expected)
  expect_equal(likelihood(xa), expected)

  rho = c(0, 0.2, 0.5)
  lf = vapply(rho, \(r) likelihood2(xf, 1, 2, rho = r), 1)
  la = vapply(rho, \(r) likelihood2(xa, 1, 2, rho = r), 1)

  expect_equal(lf, la)
})


test_that("female founder loop breakers work on X", {
  p = 0.2
  q = 1 - p

  x = ped(
    id = 1:6,
    fid = c(0, 0, 0, 2, 3, 4),
    mid = c(0, 0, 0, 1, 1, 5),
    sex = c(2, 1, 1, 1, 2, 1)
  ) |>
    addMarker(`4` = "1", `6` = "1", alleles = 1:2, afreq = c(p, q), chrom = "X")

  founder = cbind(loopBreaker = "1", child = "5")
  alternative = cbind(loopBreaker = "5", child = "6")

  xf = breakLoops(x, founder, verbose = FALSE)
  xa = breakLoops(x, alternative, verbose = FALSE)

  expect_equal(likelihood(xf), p^2 + p * q/4)
  expect_equal(likelihood(xa), p^2 + p * q/4)
})


test_that("theta counts loop-breaker copies only once", {
  p = 0.3
  theta = 0.05

  x = ped(
    id = 1:6,
    fid = c(0, 0, 0, 1, 1, 4),
    mid = c(0, 0, 0, 2, 3, 5),
    sex = c(1, 2, 2, 1, 2, 1)
  ) |>
    addMarker(`1` = "1/1", alleles = 1:2, afreq = c(p, 1 - p))

  plan = cbind(loopBreaker = "1", child = "5")
  xb = breakLoops(x, plan, verbose = FALSE)

  expect_equal(likelihood(xb), p^2)
  expect_equal(likelihood(xb, theta = theta),
               p * (theta + (1 - theta) * p))
})
