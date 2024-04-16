test_that("IMVCT() works", {
  n=100
  x=rnorm(n)
  y=2*x+rnorm(n)
  expect_lt(IMVCT(x,y,K=5,type = "linear"), 0.1)
})

