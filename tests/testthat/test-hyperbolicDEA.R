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
  
  in_WR_hyp <- hyperbolicDEA(X, Y, RTS = "vrs", WR = WR, ALPHA = 0)
  in_WR_lin <- wrDEA(X, Y, RTS = "vrs", ORIENTATION = "in", WR = WR)
  
  expect_equal(in_WR_hyp$eff, in_WR_lin$eff)
  
  expect_equal(all.equal(round(as.matrix(in_WR_hyp$lambdas),3), 
                         round(as.matrix(in_WR_lin$lambdas),3), check.attributes = FALSE), TRUE)

  expect_equal(all.equal(round(as.matrix(in_WR_hyp$mus),3), 
                         round(as.matrix(in_WR_lin$mus),3), check.attributes = FALSE), TRUE)
  
  
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

  effHyp<- hyperbolicDEA(X, Y, RTS = "vrs", ALPHA = 1, SLACK = TRUE)
  eff <- Benchmarking::dea(X, Y, RTS = "vrs", ORIENTATION = "out",SLACK = TRUE)

  logic_vec <- rowSums(round(effHyp$slack,3)) > 0
  logic_vec <- unname(logic_vec)

  expect_equal(logic_vec, eff$slack)
  

  X <- matrix(c(1,1,2,4,1.5,4,
                2,4,1,1,4,1.5), ncol = 2)
  Y <- c(1,1,1,1,1,1)
  Y <- as.matrix(Y)
  
  wr_dea <- wrDEA(X,Y, RTS = "crs", ORIENTATION = "in", SLACK = TRUE)
  hyp_dea <- hyperbolicDEA(X,Y, RTS = "crs", ALPHA = 0, SLACK = TRUE)
  BO_dea <- Benchmarking::dea(X,Y, RTS = "crs", ORIENTATION = "in", SLACK = TRUE)
  
  slack_sum <- rowSums(round(wr_dea$slack,3)) > 0
  slack_sum <- unname(slack_sum)
  
  expect_equal(all.equal(round(wr_dea$slack,3), round(hyp_dea$slack,3), check.attributes = FALSE), TRUE)
  expect_equal(slack_sum, BO_dea$slack)
  
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
  
  min_cost <- costDEA(X,Y,input_prices, RTS = "crs")
  BO_cost <- Benchmarking::cost.opt(X,Y,input_prices, RTS = "crs")
  
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
  
  max_lprofit <- lprofitDEA(X,Y,input_prices, output_prices, RTS = "crs")
  BO_profit <- Benchmarking::profit.opt(X,Y,input_prices, output_prices, RTS = "crs")
  
  expect_equal(all.equal(as.matrix(max_lprofit$opt_value), 
                         as.matrix(cbind(BO_profit$xopt, BO_profit$yopt)), 
                                   check.attributes = FALSE), TRUE)
  
  expect_equal(all.equal(as.matrix(max_lprofit$lambdas), 
                         as.matrix(BO_profit$lambda), check.attributes = FALSE), TRUE)
  
})

test_that("nlprofitDEA", {
  
  X <- matrix(c(1,2,3,3,2,1,2,2), ncol = 2)
  Y <- matrix(c(10,10,10,10), ncol = 1)
  
  input_prices <- matrix(c(2,1,2,1,2,1,1,2), ncol =  2, byrow = TRUE)
  output_prices <- matrix(c(2,2,2,2), ncol = 1)
  
  max_nlprofit <- nlprofitDEA(X,Y,input_prices, output_prices, RTS = "vrs")
  BO_profit2 <- Benchmarking::profit.opt(X,Y,input_prices, output_prices, RTS = "vrs")
  
  expect_equal(all.equal(as.matrix(max_nlprofit$opt_value), 
                         as.matrix(cbind(BO_profit2$xopt, BO_profit2$yopt)), 
                         check.attributes = FALSE), TRUE)
  
  expect_equal(all.equal(as.matrix(max_nlprofit$lambdas), 
                         as.matrix(BO_profit2$lambda), check.attributes = FALSE), TRUE)
  
})

test_that("wrDEA general test", {
  
  X <- matrix(c(1, 2, 3, 3, 6, 7, 8, 2, 1, 2), ncol = 2)
  Y <- matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3), ncol = 2)
  
  for (RTS in c("vrs", "crs", "fdh")) {
    for (ORIENTATION in c("in", "out")) {
      AO_dea <- wrDEA(X, Y, RTS = RTS, ORIENTATION = ORIENTATION)
      BO_dea <- Benchmarking::dea(X,Y,RTS = RTS, ORIENTATION = ORIENTATION)
      expect_equal(AO_dea$eff, BO_dea$eff)
    }
  }
  
  AO_dea_ndrs <- wrDEA(X, Y, RTS = "ndrs", ORIENTATION = "in")
  BO_dea_irs <- Benchmarking::dea(X,Y,RTS = "irs", ORIENTATION = "in")
  expect_equal(AO_dea_ndrs$eff, BO_dea_irs$eff)
  
  AO_dea_nirs <- wrDEA(X, Y, RTS = "nirs", ORIENTATION = "in")
  BO_dea_drs <- Benchmarking::dea(X,Y,RTS = "drs", ORIENTATION = "in")
  expect_equal(AO_dea_nirs$eff, BO_dea_drs$eff)
  
  AO_dea_ndrs <- wrDEA(X, Y, RTS = "ndrs", ORIENTATION = "out")
  BO_dea_irs <- Benchmarking::dea(X,Y,RTS = "irs", ORIENTATION = "out")
  expect_equal(AO_dea_ndrs$eff, BO_dea_irs$eff)
  
  AO_dea_nirs <- wrDEA(X, Y, RTS = "nirs", ORIENTATION = "out")
  BO_dea_drs <- Benchmarking::dea(X,Y,RTS = "drs", ORIENTATION = "out")
  expect_equal(AO_dea_nirs$eff, BO_dea_drs$eff)
  
  AO_dea_supereff <- wrDEA(X, Y, RTS = "vrs", ORIENTATION = "in", SUPEREFF = TRUE)
  BO_dea_supereff <- Benchmarking::sdea(X,Y,RTS = "vrs", ORIENTATION = "in")
  
  # Bogetoft package does not work robustly for super-efficiency -> sometimes error in the test
  # expect_equal(AO_dea_supereff$eff[AO_dea_supereff$eff > 0 & !is.infinite(AO_dea_supereff$eff)], 
  #             BO_dea_supereff$eff[BO_dea_supereff$eff > 0 & !is.infinite(BO_dea_supereff$eff)])
  
  
})

test_that("XREF YREF test", {
  
  X <- matrix(c(1, 2, 3, 3, 6, 7, 8, 2, 1, 2), ncol = 2)
  Y <- matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3), ncol = 2)
  
  XREF <- matrix(c(6, 7, 5, 5, 8, 4, 8, 5), ncol = 2)
  YREF <- matrix(c(1, 4, 5, 3, 10, 4, 2, 2), ncol = 2)
  
  AO_dea <- wrDEA(X, Y, XREF = XREF, YREF = YREF, ORIENTATION = "in", RTS = "vrs")
  BO_dea <- Benchmarking::dea(X,Y, XREF = XREF, YREF = YREF, ORIENTATION = "in", RTS = "vrs")
  AO_hyp <- hyperbolicDEA(X, Y, XREF = XREF, YREF = YREF, ALPHA = 0, RTS = "vrs")
  
  expect_equal(AO_dea$eff, BO_dea$eff)
  expect_equal(AO_hyp$eff, AO_dea$eff)
  

})






