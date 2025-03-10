
OMD = function(...) oneMarkerDistribution(..., verbose = F)

### Setup
p = 0.2; q = 1-p
g = c("1/1", "1/2", "2/2")
hw = c(p^2, 2*p*q, q^2)
PQ_ARRAY = array(c(p,q), dimnames = list(1:2))
HW_ARRAY = array(hw, dimnames = list(g))

test_that("oneMarkerDist works in empty autosomal examples", {

  ### Nuclear
  x = nuclearPed(1) |>
   addMarker(alleles = 1:2, afreq = c(p,q))
  omd1 = OMD(x, ids = 1)
  omd3 = OMD(x, ids = 3)
  omd13 = OMD(x, ids = c(1,3))
  expect_equal(omd1, HW_ARRAY)
  expect_equal(omd3, HW_ARRAY)
  expect_equal(omd13, array(hw * c(p, .5*p, 0, q, .5, p, 0, .5*q, q),
                            dim = c(3, 3), dimnames = list(g, g)))

  ### Singleton
  s = singleton("a") |>
    addMarker(alleles=letters[1:2], afreq = c(p,q))
  omd = OMD(s, ids="a")
  expect_equal(omd, setNames(HW_ARRAY, c("a/a","a/b","b/b")))
})

test_that("oneMarkerDist works in empty X-linked examples", {

  ### Nuclear
  x = nuclearPed(1, sex = 2) |>
    addMarker(alleles = 1:2, afreq = c(p,q), chrom = 23)
  omd1 = OMD(x, ids = 1)
  omd3 = OMD(x, ids = 3)
  omd13 = OMD(x, ids = c(1,3))
  expect_equal(omd1, PQ_ARRAY)
  expect_equal(omd3, HW_ARRAY)
  expect_equal(omd13, array(c(p^2, 0, p*q, p*q, 0, q^2), dim = c(2,3),
                            dimnames = list(1:2, g)))

  ### Singletons
  s = singletons(1:2, sex = 1:2) |>
    addMarker(alleles = letters[1:2], afreq = c(p,q), chrom = "X")

  omdx1 = OMD(s, ids = 1)
  omdx2 = OMD(s, ids = 2)

  expect_equal(omdx1, setNames(PQ_ARRAY, letters[1:2]))
  expect_equal(omdx2, setNames(HW_ARRAY, c("a/a","a/b","b/b")))
})

test_that("oneMarkerDist works in conditional nuclear example", {

  x = nuclearPed(1, sex = 2) |>
    addMarker(`3` = 2, alleles = 1:2, afreq = c(p,q)) |>
    addMarker(`3` = 2, alleles = 1:2, afreq = c(p,q), chrom = 23)

  res1 = HW_ARRAY; res1[] = c(0, p, q)
  expect_equal(OMD(x, ids = 1, marker = 1), res1)

  res2 = PQ_ARRAY; res2[] = c(0, 1)
  expect_equal(OMD(x, ids = 1, marker = 2), res2)
  expect_equal(OMD(x, ids = 2, marker = 2), res1)
})


test_that("oneMarkerDist on X is independent of male-male", {
  x = nuclearPed(1) |>
    addMarker("1" = 1, alleles = 1:2, afreq = c(p,q), chrom = "X")

  omd = OMD(x, ids = 3) |> c() # convert array to vector
  expect_equal(omd, c("1"=p, "2"=q))
})

test_that("oneMarkerDist with multiple components", {
  # Trivial
  x = singletons(1:2) |> addMarker(alleles = 1:2, afreq = c(p,q))
  y = nuclearPed(3) |> transferMarkers(from = x, to = _)
  expect_equal(OMD(x, ids = 1:2), OMD(y, 1:2))

  # Non-trivial example
  x = list(nuclearPed(ch = 5), nuclearPed(fa = 3, mo = 4, ch = 7))
  y = addSon(x, 3:2) |> reorderPed()
  x = x |> addMarker(`1` = "1/-", `2`="1/1", `3` = "2/2", alleles = 1:2, afreq = c(p,q))
  y = y |> addMarker(`1` = "1/-", `2`="1/1", `3` = "2/2", alleles = 1:2, afreq = c(p,q))
  # plotPedList(list(x,y), marker = 1, cex = 1.5, widths = c(1,1,2))
  expect_equal(OMD(x, ids = c(1,5,4,7)), OMD(y, ids = c(1,5,4,7)))
})
