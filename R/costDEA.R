#' @title Cost DEA model
#' @description Cost DEA model optimizing the input allocation with given prices.
#' It returns the estimated lambdas as well as the optimal values for inputs.
#' 
#' @seealso [Benchmakring::cost.opt] for a similar function
#'
#' @param X Matrix or dataframe with DMUS as rows and inputs as columns
#' @param Y Matrix or dataframe with DMUs as rows and outputs as columns
#' @param pX Matrix or dataframe with prices for each DMU and input.
#' Therefore it mus have the same dimensions as X.
#' @param RTS Character string indicating the returns-to-scale, e.g. "crs", "vrs"
#' @return A list object containing cost optimal inputs and lambdas indiceting
#' the the peer for optimal input allocation. Additionally, it returns the cost 
#' efficiency as the ratio of the optimal cost to the observed cost.
#' @examples
#' X <- matrix(c(1,2,3,3,2,1,2,2), ncol = 2)
#' Y <- matrix(c(1,1,1,1), ncol = 1)
#'
#' pX <- matrix(c(2,1,2,1,2,1,1,2), ncol =  2, byrow = TRUE)
#' 
#'
#' cost_eff_input <- costDEA(X,Y,pX)
#'
#' @import lpSolveAPI
#' @export

costDEA <- function(X, Y, pX, RTS = "vrs") {
  
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
  if (!identical(dim(X), dim(pX))) {
    stop("Dimensions of pX and X are not the same.")
  }
  
  RTS <- tolower(RTS)
  possible_rts <- c("crs", "vrs")
  if (!(RTS %in% possible_rts)){
    stop("Scale of returns not implemented:", RTS)
  }
  
  # Change data structure to uniformly be matrices
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  pX <- as.matrix(pX)
  
  # Change data structure to add columns
  in_out_data <- rbind(t(X),t(Y))
  
  # storage for results
  lambdas <- data.frame()
  opt_value <- data.frame()
  
  # Solve for all DMUs
  for (i in 1:ncol(in_out_data)) {
    
    # Create the linear programming problem
    # Constraints are nrow and columns are ncol plus
    # the number of rows for optimizing inputs
    cost_lp <- make.lp(nrow = nrow(in_out_data), ncol = ncol(in_out_data)+ncol(X))
    
    # Change to minimization problem
    lp.control(cost_lp, sense="min")
    
    # Set the coefficients of the linear programming problem (frontier)
    for (column in 1:ncol(in_out_data)) {
      set.column(cost_lp, column, in_out_data[,column])
    }
    
    # Set decision variables for the objective function
    # optimizing each respective input
    for (j in 1:ncol(X)){
      column <- c(rep(0, nrow(in_out_data)))
      column[j] <- -1
      set.column(cost_lp, ncol(in_out_data)+j, column)
    }
    
    # Set the objective function
    set.objfn(cost_lp, c(rep(0,ncol(in_out_data)),pX[i,]))
    
    # Set constraint-types
    set.constr.type(cost_lp, c(rep("<=", ncol(X)), rep(">=", ncol(Y))))
    
    # Set the lower bounds of the decision variables
    set.bounds(cost_lp, lower = c(rep(0,(ncol(in_out_data)+ncol(X)))))
    
    # Set the right-hand side of the constraints
    set.rhs(cost_lp, c(rep(0, ncol(X)), Y[i,]))
    
    # Set the sum of the lambdas to 1
    if (RTS == "vrs") {
      add.constraint(cost_lp, c(rep(1, ncol(in_out_data))), 
                     indices = c(1:ncol(in_out_data)), "=", 1)
      set.bounds(cost_lp, upper = c(rep(1,(ncol(in_out_data)))), columns = c(1:ncol(in_out_data)))
    }
    
    # Solve the linear programming problem
    solve(cost_lp)
    
    # Get the optimal solution and store results in a data frame
    variables <- get.variables(cost_lp)
    lambdas <- rbind(lambdas, variables[1:ncol(in_out_data)])
    opt_value <- rbind(opt_value, variables[(ncol(in_out_data)+1):(ncol(in_out_data)+ncol(X))])
    
  }
  
  colnames(lambdas) <- c(paste("Lambda",1:ncol(in_out_data),sep=""))
  colnames(opt_value) <- c(paste("X",1:ncol(X),sep=""))
  
  # Estimate cost efficiency
  cost_eff <- rowSums(opt_value * pX) / rowSums(X * pX)
  
  
  # Return the results
  return(list(lambdas = lambdas, opt_value = opt_value, cost_eff = cost_eff))
}