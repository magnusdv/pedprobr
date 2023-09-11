
test_that("singleton with dropout", {
  p = 0.2
  afr = c("1" = p, "2" = (1-p)/2, "3" = (1-p)/2)
  s = singleton(1) |> addMarker(geno = "1/1", afreq = afr)

  d = 0.1
  expect_equal(likelihood(s, dropout = 0), p^2)
  expect_equal(likelihood(s, dropout = d), (1-d^2)*p^2 + d*(1-d)*2*p*(1-p))

  f = 0.2
  Lf = (1-d^2)*(f*p + (1-f)*p^2) + d*(1-d)*(1-f)*2*p*(1-p)

  s = setFounderInbreeding(s, 1, f)
  expect_equal(likelihood(s, dropout = d), Lf)
})

test_that("trio with dropout", {
  x = nuclearPed(fa = "fa", mo = "mo", ch = "ch")
  p = 0.2
  afr = c("1" = p, "2" = p, "3" = 1-2*p)

  d = 0.1
  L = (1-d^2)*p^2 + d*(1-d)*2*p*(1-p)

  x1 = x |> addMarker(fa = "1/1", afreq = afr)
  expect_equal(likelihood(x1, dropout = d), L)
  expect_equal(likelihood(x1, dropout = d, lump = F), L)

  x2 = x |> addMarker(ch = "1/1", afreq = afr)
  expect_equal(likelihood(x2, dropout = d), L)
  expect_equal(likelihood(x2, dropout = c(fa = d)), p^2)

  x3 = x |> addMarker(fa = "1/1", mo = "2/2", afreq = afr)
  expect_equal(likelihood(x3, dropout = d), L^2)

  x4 = x |> addMarker(fa = "1/1", mo = "2/2", afreq = afr)
  expect_equal(likelihood(x4, dropout = c(fa = d)), L * p^2) # dropout only in father
})

test_that("inbreeding with dropout", {
  x = nuclearPed() |> addSon(2:3)

  p = 0.2
  afr = c("1" = p, "2" = p, "3" = 1-2*p)

  d = 0.1
  f = 0.25 # ribd::inbreeding(x, 4)
  L = (1-d^2)*p^2 + d*(1-d)*2*p*(1-p) # without inbr
  Lf = (1-d^2)*(f*p + (1-f)*p^2) + d*(1-d)*(1-f)*2*p*(1-p)

  x1 = x |> addMarker("3" = "1/1", afreq = afr)
  expect_equal(likelihood(x1, dropout = d), L)

  x2 = x |> addMarker("4" = "1/1", afreq = afr)
  expect_equal(likelihood(x2, dropout = d), Lf)
  expect_equal(likelihood(x2, dropout = d, lump = F), Lf)

  x2 = x |> addMarker("4" = "1/1", afreq = afr)
  expect_equal(likelihood(x2, dropout = d), Lf)

})

