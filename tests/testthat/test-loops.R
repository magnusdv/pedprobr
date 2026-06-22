

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

  expected = p^2 + 3/8 * p * (1 - p)

  expect_equal(likelihood(x, loopBreakers = c(3,3)), expected)
  expect_equal(likelihood(x, loopBreakers = c(3,5)), expected)

  rho = c(0, 0.2, 0.5)
  lr = vapply(rho, \(r) likelihood2(x, 1, 2, rho = r, loopBreakers = c(3,3)), 1)
  la = vapply(rho, \(r) likelihood2(x, 1, 2, rho = r, loopBreakers = c(3,5)), 1)

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

  expected = c(p[1]^2 + 1/8 * p[1] * (1 - p[1]),
               2 * p[2] * (1 - p[2]) * 7/8)

  expect_equal(likelihood(x, loopBreakers = 1), expected)
  expect_equal(likelihood(x, loopBreakers = 5), expected)

  rho = c(0, 0.2, 0.5)
  lf = vapply(rho, \(r) likelihood2(x, 1, 2, rho = r, loopBreakers = 1), 1)
  la = vapply(rho, \(r) likelihood2(x, 1, 2, rho = r, loopBreakers = 5), 1)

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

  e = p^2 + p * q/4
  expect_equal(likelihood(x, loopBreakers = 1), e)
  expect_equal(likelihood(x, loopBreakers = 5), e)
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

  expect_equal(likelihood(x, loopBreakers = 1), p^2)
  expect_equal(likelihood(x, loopBreakers = 1, theta = theta),
               p * (theta + (1 - theta) * p))
})

test_that("verbose looped likelihood shows loop complexity", {
  x = ped(
    id = 1:6,
    fid = c(0, 0, 1, 1, 1, 5),
    mid = c(0, 0, 2, 2, 3, 3),
    sex = c(1, 2, 2, 2, 1, 1)
  ) |>
    addMarker(`6` = "1/1", alleles = 1:2)
  # plot(x)
  expect_message(likelihood(x, loopBreakers = c(3, 3), .diagnostics = TRUE),
                 "Genotype combinations")
})


test_that("automatic loop breaking uses X-specific scores", {
  x = ped(
    id = 1:6,
    fid = c(0, 0, 0, 2, 3, 4),
    mid = c(0, 0, 0, 1, 1, 5),
    sex = c(2, 1, 1, 1, 2, 1)
  ) |>
    addMarker(`6` = "1", alleles = 1:10, chrom = "X")

  y = pedprobr:::.breakLoops(x, verbose = FALSE)
  lb = y$ID[y$LOOP_BREAKERS[, "orig"]]

  expect_equal(lb, "4")
})
