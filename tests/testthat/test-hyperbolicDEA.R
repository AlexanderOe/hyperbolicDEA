library(Benchmarking)
library(hyperbolicDEA)


test_that("VRS", {
  X <- c(2,3,3,8,9)
  Y <- c(5,5,6,7,10)
  eff_comparison <- Benchmarking::dea(X, Y, RTS="vrs",ORIENTATION="graph")$eff
  eff <- hyperbolicDEA(X, Y, RTS = "vrs")$eff

  expect_equal(round(eff,3), round(eff_comparison,3))
})

test_that("CRS", {
  X <- c(2,3,3,8,9)
  Y <- c(5,5,6,7,10)
  eff_comparison <- Benchmarking::dea(X, Y, RTS="crs",ORIENTATION="graph")$eff
  eff <- hyperbolicDEA(X, Y, RTS = "crs")$eff

  expect_equal(round(eff,3), round(eff_comparison,3))
})

test_that("Scaling", {
  X1 <- c(1,2,4,7,6,7)*100000000000
  Y1 <- c(1,3,2,5,4,6)*100000000000
  X2 <- c(1,3,4,5,6,4)/1000000000
  Y2 <- c(1,2,4,4,4,5)/1000000000

  X <- cbind.data.frame(X1,X2)
  Y <- cbind.data.frame(Y1,Y2)

  eff_comparison <- Benchmarking::dea(X, Y, RTS="vrs",ORIENTATION="graph")$eff
  eff <- hyperbolicDEA(X, Y, RTS = "vrs")$eff

  expect_equal(round(eff,3), round(eff_comparison,3))
})

test_that("Weight Restrictions", {
  X <- c(1,2,4,5,6,7)
  Y <- c(1,3,2,5,4,6)

  # Weight restrictions two reduction of input -> Three reduction of output
  # Similar to NIRS
  WR <- rbind(c(-3,-2))

  eff_WR <- hyperbolicDEA(X, Y, RTS = "vrs", WR = WR)$eff
  eff_NIRS <- hyperbolicDEA(X, Y, RTS = "nirs")$eff

  expect_equal(eff_WR, eff_NIRS)
})
