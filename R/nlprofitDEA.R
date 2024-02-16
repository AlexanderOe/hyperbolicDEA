#' @title Non-linear profit DEA model
#' @description Non-linear profit DEA model optimizing the ratio of cost over revenue
#' 
#' @seealso Julia package BenchmarkingEconomicEfficiency.jl and the function deaprofitability()
#' for non-linear profit DEA model
#'
#' @param X Matrix or dataframe with DMUS as rows and inputs as columns
#' @param Y Matrix or dataframe with DMUs as rows and outputs as columns
#' @param pX Matrix or dataframe with prices for each DMU and input.
#' Therefore it mus have the same dimensions as X.
#' @param pY Matrix or dataframe with prices for each DMU and output.
#' Therefore it mus have the same dimensions as Y.
#' @param RTS Character string indicating the returns-to-scale, e.g. "crs", "vrs"
#' @return A list object containing optimal inputs and outputs, lambdas indicating
#' the peers for optimal allocation and profitability score as the ratio of revenue
#' over cost for optimal and observed allocation.
#' @examples
#' X <- matrix(c(1,2,3,3,2,1,2,2), ncol = 2)
#' Y <- matrix(c(1,1,1,1), ncol = 1)
#'
#' pX <- matrix(c(2,1,2,1,2,1,1,2), ncol =  2, byrow = TRUE)
#' pY <- matrix(c(1,1,1,1), ncol = 1)
#'
#' max_prof_nolin<- nlprofitDEA(X,Y,pX,pY)
#'
#' @import nloptr
#' @export


nlprofitDEA <- function(X, Y, pX, pY, RTS = "vrs"){

  # Check arguments given by user
  if (!is.matrix(X) && !is.data.frame(X) && !is.numeric(X)){
    stop("X must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(Y) && !is.data.frame(Y) && !is.numeric(Y)){
    stop("Y must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(pX) && !is.data.frame(pX) && !is.numeric(pX)){
    stop("pX must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(pY) && !is.data.frame(pY) && !is.numeric(pY)){
    stop("pY must be a numeric vector, matrix or dataframe")
  }
  if (!identical(dim(X), dim(pX))) {
    stop("Dimensions of pX and X are not the same.")
  }
  if (!identical(dim(Y), dim(pY))) {
    stop("Dimensions of pY and Y are not the same.")
  }

  RTS <- tolower(RTS)
  possible_rts <- c("crs", "vrs")
  if (!(RTS %in% possible_rts)){
    stop("Scale of returns not implemented:", RTS)
  }

  # uniform data structure
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  pX <- as.matrix(pX)
  pY <- as.matrix(pY)


  # storage for results
  lambdas <- data.frame()
  opt_value <- data.frame()

  # Solve for all DMUs
  for (i in 1:nrow(X)) {

    # Definition of the objective function. We take cost over revenue to
    # have a minimization problem instead of a maximization problem
    eval_f <- function(controls){
      cost <- sum(controls[(nrow(X)+1):(nrow(X)+ncol(X))]*pX[i,])
      revenue <- sum(controls[(nrow(X)+ncol(X)+1):(nrow(X)+ncol(X)+ncol(Y))]*pY[i,])
      return(cost/revenue)
    }

    # Definition of the gradient of the objective function
    eval_grad_f <- function(controls){
      grad <- c(rep(0, nrow(X)))
      cost <-  sum(controls[(nrow(X)+1):(nrow(X)+ncol(X))]*pX[i,])
      revenue <- sum(controls[(nrow(X)+ncol(X)+1):(nrow(X)+ncol(X)+ncol(Y))]*pY[i,])
      # Partial derivatives with respect to the inputs
      for (m in 1:ncol(X)){
        grad <- c(grad, pX[i,m]/revenue)
      }
      # Partial derivatives with respect to the outputs
      for (n in 1:ncol(Y)){
        grad <- c(grad, -pY[i,n]*cost/revenue^2)
      }
      return(grad)
    }

    # Definition of the inequality constraints (strored in a vector)
    eval_g_ineq <- function(controls){
      constr <- c()
      for (j in c(1:ncol(X))){
        constraint_x <- controls[1:nrow(X)]%*%X[,j] - controls[nrow(X)+j]
        constr <- c(constr, constraint_x)
      }

      for (k in c(1:ncol(Y))){
        constraint_y <- -(controls[1:nrow(X)])%*%Y[,k] + controls[nrow(X)+ncol(X)+k]
        constr <- c(constr, constraint_y)
      }

      return(constr)
    }

    # Definition of the jacobian of the inequality constraints
    # stored in a matrix
    eval_jac_g_ineq <- function(controls){
      jacobian_X <- matrix(0, ncol = length(controls), nrow = ncol(X))
      jacobian_Y <- matrix(0, ncol = length(controls), nrow = ncol(Y))

      for (j in seq_len(ncol(X))) {
        jacobian_X[j, 1:nrow(X)] <- X[, j]
        jacobian_X[j, nrow(X) + j] <- -1
      }

      for (k in seq_len(ncol(Y))) {
        jacobian_Y[k, 1:nrow(X)] <- -Y[, k]
        jacobian_Y[k, nrow(X) + ncol(X) + k] <- 1
      }

      jacobian <- rbind(jacobian_X, jacobian_Y)
      return(jacobian)
    }

    # Equality constraint and jacobian for VRS assumption
    if (RTS == "vrs"){
      eval_g_eq <- function(controls){
        constr <- c(sum(controls[1:nrow(X)]) - 1)
        return(constr)
      }
      eval_jac_g_eq <- function(controls){
        jacobian <- c(rep(1, nrow(X)), rep(0, ncol(X)+ncol(Y)))
        return(jacobian)
      }
    } else{
      eval_g_eq <- NULL
      eval_jac_g_eq <- NULL
    }

    # Define the bounds for the controls in VRS lambdas constraint
    # to be between 0 and 1
    lb <- c(rep(0, (nrow(X)+ncol(X)+ncol(Y))))
    ub <- c(rep(1, nrow(X)), rep(Inf, (ncol(X)+ncol(Y))))

    # Start searching at best ratio with given prices of the DMU
    # Transposing for vector multiplication
    cost <- rowSums(t(t(X)*pX[i,]))
    revenue <- rowSums(t(t(Y)*pY[i,]))

    # Calculate the ratios
    ratios <- revenue/cost

    # Find the index of the observation with the highest ratio
    index <- which.max(ratios)
    start_values <- rep(0, nrow(X))
    start_values[index] <- 1
    controls0 <- c(start_values, X[index,], Y[index,])

    # Number of inequality constraints
    no_constr_ineq <- ncol(X)+ncol(Y)

    # Define the options for the optimization algorithm
    opts <-list("algorithm"="NLOPT_LD_SLSQP",
                "xtol_rel"=1.0e-10,
                "maxeval"=1000000,
                "tol_constraints_ineq"=rep(1.0e-10,no_constr_ineq))
    res <-nloptr(x0=controls0,
                 eval_f=eval_f,
                 eval_grad_f = eval_grad_f,
                 lb=lb,
                 ub=ub,
                 eval_g_ineq =eval_g_ineq,
                 eval_jac_g_ineq = eval_jac_g_ineq,
                 eval_g_eq = eval_g_eq,
                 eval_jac_g_eq = eval_jac_g_eq,
                 opts=opts)

    # Store the results
    lambdas <- rbind(lambdas, res$solution[1:nrow(X)])
    opt_value <- rbind(opt_value, res$solution[(nrow(X)+1):(nrow(X)+ncol(X)+ncol(Y))])

  }

  # Add column names
  colnames(lambdas) <- c(paste("Lambda",1:nrow(X),sep=""))
  colnames(opt_value) <- c(paste("X",1:ncol(X),sep=""), paste("Y",1:ncol(Y),sep=""))
  
  # Profitability efficiency
  profit_eff <- (rowSums(Y*pY)/rowSums(X*pX)) / (rowSums(pY*opt_value[,(ncol(X)+1):(ncol(X)+ncol(Y))])/rowSums(pX*opt_value[,1:ncol(X)]))


  # Return the results
  return(list(lambdas = lambdas, opt_value = opt_value, profit_eff = profit_eff))

}
