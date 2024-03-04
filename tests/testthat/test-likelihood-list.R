
test_that("pedlist with empty markers give likelihood 1", {
  x = list(singleton("NN"), nuclearPed(1)) |> addMarker(alleles = 1:2, name = "m1")
  expect_equal(likelihood(x, "m1"), 1)
})

test_that("pedlist with HW markers gives correct likelihood", {
  p = 0.9; q= 1-p

  x = list(nuclearPed(1), singleton("nn")) |>
    addMarker("3" = "1/2", nn = "1/1", alleles = 1:2, afreq = c(p,q))
  # plot(x, mark = 1)

  expect_equal(likelihood(x), 2*p*q * p^2)
})

