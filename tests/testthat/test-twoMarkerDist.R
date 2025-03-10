
TMD = function(x, ...)  twoMarkerDistribution(x, ..., verbose = F)

### Setup
p = 0.2; q = 1-p
g = c("1/1", "1/2", "2/2")
hw = c(p^2, 2*p*q, q^2)
PQ_ARRAY = array(c(p,q), dimnames = list(1:2))
HW_ARRAY = array(hw, dimnames = list(g))

test_that("twoMarkerDist catches errors", {
  x = nuclearPed() |> addMarker() |> addMarker()
  expect_error(TMD(x, id = 1:2, rho = 0.25),
               "Argument `id` must have length 1: 1, 2")
  expect_error(TMD(x, marker1 = 3, id = 1, rho = 0.25),
               "Marker index out of range: 3")
  expect_error(TMD(x, 1:2, id = 1, rho = 0.25),
               "Argument `marker1` must have length 1: 1, 2")
  expect_error(TMD("foo", id = 1, rho = 0.25),
               "Input is not a pedigree")
  expect_error(TMD(x, id = 1, rho = -0.25),
               "Argument `rho` cannot be negative")
  expect_error(TMD(x, id = 1, rho = 0.75),
               "Argument `rho` cannot exceed 0.5")
  expect_error(TMD(x, id = 1), "Argument `rho` is missing")
})

test_that("twoMarkerDist works in empty autosomal examples", {

  ### Nuclear
  x = nuclearPed(1) |> addMarker(alleles = 1:2, afreq = c(p,q))

  expect_equal(TMD(x, 1,1, id = 1, rho = 0),
               HW_ARRAY %o% HW_ARRAY)
  expect_equal(TMD(x, 1,1, id=3, rho=0.25),
               HW_ARRAY %o% HW_ARRAY)

  ### Singleton
  s = singleton("a") |>
    addMarker(alleles = 1:2,          afreq = c(p,q)) |>
    addMarker(alleles = letters[1:2], afreq = c(p,q))

  HW_ab = setNames(HW_ARRAY, c("a/a","a/b","b/b"))

  expect_equal(TMD(s, id="a", rho=0), HW_ARRAY %o% HW_ab)
})

test_that("twoMarkerDist works in empty X-linked examples", {

  ### Nuclear
  x = nuclearPed(1, sex = 2) |>
    addMarker(alleles=1:2, afreq = c(p,q), chrom = 23)

  expect_equal(TMD(x, 1,1, id=1, rho=0),
               PQ_ARRAY %o% PQ_ARRAY)
  expect_equal(TMD(x, 1,1, id=3, rho=0.25),
               HW_ARRAY %o% HW_ARRAY)

  ### Singleton
  s = singletons(1:2, sex = 1:2) |>
    addMarker(alleles = 1:2, afreq = c(p,q), chrom = 23) |>
    addMarker(alleles = letters[1:2], afreq = c(p,q), chrom = 23)

  HW_ab = setNames(HW_ARRAY, c("a/a","a/b","b/b"))
  PQ_ab = setNames(PQ_ARRAY, c("a","b"))

  expect_equal(TMD(s, id = 1, rho = 0), PQ_ARRAY %o% PQ_ab)
  expect_equal(TMD(s, id = 2, rho = 0), HW_ARRAY %o% HW_ab)
})

test_that("twoMarkerDist works in conditional nuclear example", {

  x = nuclearPed(1, sex = 2) |>
    addMarker(`3` = 2, alleles = 1:2, afreq = c(p,q)) |>
    addMarker(`3` = 2, alleles = 1:2, afreq = c(p,q), chrom = 23)

  # Autosomal
  h = HW_ARRAY; h[] = c(0, p, q); res1 = h %o% h
  expect_equal(TMD(x, 1,1, id = 1, rho = 0), res1)

  # X-linked
  res2 = PQ_ARRAY %o% PQ_ARRAY; res2[] = c(0,0,0,1)
  expect_equal(TMD(x, 2,2, id = 1, rho = 0), res2)
  expect_equal(TMD(x, 2,2, id = 2, rho = 0), res1)
})

test_that("recombination rate is recovered", {
  rho = 0.15

  x = linearPed(2, sex = 2:1) |>
    addMarker(`1`=1, `2`=2, `3`=1, `5`=1, alleles = 1:2, afreq = c(p,q)) |>
    addMarker(`1`=1, `2`=2, `3`=1, `5`=0, alleles = 1:2, afreq = c(p,q))
  # plot(x, marker = 1:2)

  expect_equal(TMD(x, id = 5, rho = rho)["1/1", "1/2"], rho)

  ### Same on X
  y = setChrom(x, marker = 1:2, chrom = "X")
  # plot(y, marker = 1:2)

  expect_equal(TMD(y, id = 5, rho = rho)["1", "2"], rho)
})

