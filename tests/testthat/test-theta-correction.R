
test_that("theta works in two-generation pedigrees", {
  x = linearPed(2)
  m = marker(x, `5` = "A/A", afreq = c(A = 0.3, B = 0.7))

  expect_equal(likelihood(x, m), likTheta(x, m, theta = 0))
  expect_equal(likelihood(x, m, theta = 0.01), 0.0921)

  y = ancestralPed(2)
  m2 = marker(y, `7` = "A/A", afreq = c(A = 0.3, B = 0.7))
  expect_equal(likelihood(y, m2), likTheta(y, m2, theta = 0))
  expect_equal(likelihood(y, m2, theta = 0.01), 0.0921)
})



test_that("theta-correction gives error if founders are inbred", {
  x = singleton(1)
  founderInbreeding(x, 1) = 0.5
  m = marker(x, afreq = c(A = 0.3, B = 0.7), "1" = "A/B")

  expect_error(likelihood(x, m, theta = 0.05), "Theta correction cannot be used in pedigrees with inbred founders")
})
