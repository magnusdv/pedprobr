
test_that("theta-correction gives error if founders are inbred", {
  x = singleton(1)
  founderInbreeding(x, 1) = 0.5
  m = marker(x, afreq = c(A = 0.3, B = 0.7), "1" = "A/B")

  expect_error(likelihood(x, m, theta = 0.05), "Theta correction cannot be used in pedigrees with inbred founders")
})
