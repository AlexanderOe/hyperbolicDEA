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

test_that("ALPHA=1 output orientation", {
  X <- c(1,2,4,5,6,7)
  Y <- c(1,3,2,5,4,6)

  eff_alpha<- hyperbolicDEA(X, Y, RTS = "vrs", ALPHA = 1)$eff
  eff <- 1/Benchmarking::dea(X, Y, RTS = "vrs", ORIENTATION = "out")$eff

  expect_equal(round(eff_alpha,3), round(eff,3))
})

test_that("ALPHA=0 input orientation", {
  X <- c(1,2,4,5,6,7)
  Y <- c(1,3,2,5,4,6)

  eff_alpha<- hyperbolicDEA(X, Y, RTS = "vrs", ALPHA = 0)$eff
  eff <- Benchmarking::dea(X, Y, RTS = "vrs", ORIENTATION = "in")$eff

  expect_equal(round(eff_alpha,3), round(eff,3))
})

test_that("SLACK", {
  X <- c(1,1,2,3)
  Y <- c(1,2,4,3)

  effHyp<- hyperbolicDEA(X, Y, RTS = "vrs", ALPHA = 1, SLACK = T)
  eff <- Benchmarking::dea(X, Y, RTS = "vrs", ORIENTATION = "out",SLACK = T)

  logic_vec <- rowSums(round(effHyp$slack,3)) > 0
  logic_vec <- unname(logic_vec)

  expect_equal(logic_vec, eff$slack)
})

test_that("fdh", {
  X <- c(1,2,4,5,6,7)
  Y <- c(1,3,2,5,4,6)

  # Hyperbolic Orientation
  effHyp<- hyperbolicDEA(X, Y, RTS = "fdh", ALPHA = 0.5)
  eff <- Benchmarking::dea(X, Y, RTS = "fdh", ORIENTATION = "graph")
  expect_equal(round(effHyp$eff, 3), round(eff$eff, 3))

  # Output Orientation
  effHyp_out<- hyperbolicDEA(X, Y, RTS = "fdh", ALPHA = 1)
  eff_out <- Benchmarking::dea(X, Y, RTS = "fdh", ORIENTATION = "out")
  expect_equal(round(effHyp_out$eff, 3), round(1/eff_out$eff, 3))

  # Output Orientation
  effHyp_in<- hyperbolicDEA(X, Y, RTS = "fdh", ALPHA = 0)
  eff_in <- Benchmarking::dea(X, Y, RTS = "fdh", ORIENTATION = "in")
  expect_equal(round(effHyp_in$eff, 3), round(eff_in$eff, 3))

  # Multidimensional
  x <- c(1,1,1,1)
  y <- matrix(c(1,1,3,4,
                1,3,2,1), ncol = 2)
  est_hyp <- hyperbolicDEA(x, y, RTS = "fdh", ALPHA = 1)
  est <- Benchmarking::dea(x, y, RTS = "fdh", ORIENTATION = "out")

  expect_equal(round(est_hyp$eff, 3), round(1/est$eff, 3))


  # Testing Lambdas
  X <- matrix(c(1,2,4,5,7,7), ncol=1)
  Y <- as.matrix(c(1,3,2,5,4,6))

  effHyp <- hyperbolicDEA(X, Y, RTS = "fdh", ALPHA = 1)
  effBO <- Benchmarking::dea(X, Y, RTS = "fdh", ORIENTATION = "out")
  expect_equal(all.equal(effHyp$lambdas, effBO$lambda, check.attributes = FALSE), TRUE)

})


test_that("multiple weight restricitons", {

  X1 <- c(4,4,5,6,7)
  X2 <- c(2,4,3,7,2)
  Y <- c(1,1,1,1,1)

  WR <- matrix(c(1/2,1), nrow = 1)
  WR_hyp <- matrix(c(0,-1,2,
                     0,1,-1),
                   ncol = 3, nrow = 2, byrow = TRUE)

  BO_crs_WR <- Benchmarking::dea.dual(cbind(X1, X2), Y, RTS = "crs", ORIENTATION = "in", DUAL = WR)

  hyp_crs_WR <- hyperbolicDEA(cbind(X1, X2), Y, RTS = "crs", ALPHA = 0, WR = WR_hyp)
  expect_equal(round(BO_crs_WR$eff, 3), round(hyp_crs_WR$eff, 3))


})

test_that("costDEA", {
  
  X <- matrix(c(1,2,3,3,2,1,2,2), ncol = 2)
  Y <- matrix(c(1,1,1,1), ncol = 1)
  
  input_prices <- matrix(c(2,1,2,1,2,1,1,2), ncol =  2, byrow = TRUE)
  
  min_cost <- costDEA(X,Y,input_prices)
  BO_cost <- Benchmarking::cost.opt(X,Y,input_prices)
  
  expect_equal(all.equal(as.matrix(min_cost$opt_value), 
                         as.matrix(BO_cost$xopt), check.attributes = FALSE), TRUE)
  
  expect_equal(all.equal(as.matrix(min_cost$lambdas), 
                         as.matrix(BO_cost$lambda), check.attributes = FALSE), TRUE)
  
})

test_that("lprofitDEA", {
  
  X <- matrix(c(1,2,3,3,2,1,2,2), ncol = 2)
  Y <- matrix(c(1,1,1,1), ncol = 1)
  
  input_prices <- matrix(c(2,1,2,1,2,1,1,2), ncol =  2, byrow = TRUE)
  output_prices <- matrix(c(1,1,1,1), ncol = 1)
  
  max_profit <- lprofitDEA(X,Y,input_prices, output_prices)
  BO_revenue <- Benchmarking::profit.opt(X,Y,input_prices, output_prices)
  
  expect_equal(all.equal(as.matrix(max_profit$opt_value), 
                         as.matrix(cbind(BO_revenue$xopt, BO_revenue$yopt)), 
                                   check.attributes = FALSE), TRUE)
  
  expect_equal(all.equal(as.matrix(max_profit$lambdas), 
                         as.matrix(BO_revenue$lambda), check.attributes = FALSE), TRUE)
  
})

test_that("nlprofitDEA", {
  
  X <- matrix(c(1,2,3,3,2,1,2,2), ncol = 2)
  Y <- matrix(c(1,1,1,1), ncol = 1)
  
  input_prices <- matrix(c(2,1,2,1,2,1,1,2), ncol =  2, byrow = TRUE)
  output_prices <- matrix(c(1,1,1,1), ncol = 1)
  
  max_profit <- nlprofitDEA(X,Y,input_prices, output_prices)
  BO_revenue <- Benchmarking::profit.opt(X,Y,input_prices, output_prices)
  
  expect_equal(all.equal(as.matrix(max_profit$opt_value), 
                         as.matrix(cbind(BO_revenue$xopt, BO_revenue$yopt)), 
                         check.attributes = FALSE), TRUE)
  
  expect_equal(all.equal(as.matrix(max_profit$lambdas), 
                         as.matrix(BO_revenue$lambda), check.attributes = FALSE), TRUE)
  
})


