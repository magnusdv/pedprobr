
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

test_that("pedlists preserve component-specific repeated loop breakers", {
  x1 = ped(
    id = 1:6,
    fid = c(0, 0, 1, 1, 1, 5),
    mid = c(0, 0, 2, 2, 3, 3),
    sex = c(1, 2, 2, 2, 1, 1)
  ) |>
    addMarker(`6` = "1/1", alleles = 1:2) |>
    addMarker(`6` = "1/1", alleles = 1:2)

  x2 = ped(
    id = 11:16,
    fid = c(0, 0, 11, 11, 11, 15),
    mid = c(0, 0, 12, 12, 13, 13),
    sex = c(1, 2, 2, 2, 1, 1)
  ) |>
    addMarker(`16` = "1/1", alleles = 1:2) |>
    addMarker(`16` = "1/1", alleles = 1:2)

  lb = c(3, 3, 13, 13)

  expect_equal(likelihood(list(x1, x2), loopBreakers = lb),
               likelihood(x1, loopBreakers = c(3, 3)) *
                 likelihood(x2, loopBreakers = c(13, 13)))

  expect_equal(likelihood2(list(x1, x2), 1, 2, rho = 0.25, loopBreakers = lb),
               likelihood2(x1, 1, 2, rho = 0.25, loopBreakers = c(3, 3)) *
                 likelihood2(x2, 1, 2, rho = 0.25, loopBreakers = c(13, 13)))
})
