#' @title Linear profit DEA model
#' @description Linear profit DEA model optimizing the difference of cost to revenue.
#' It returns the estimated lambdas as well as the optimal values for inputs and outputs.
#'
#' @param X Matrix or dataframe with DMUS as rows and inputs as columns
#' @param Y Matrix or dataframe with DMUs as rows and outputs as columns
#' @param input_prices Matrix or dataframe with prices for each DMU and input.
#' Therefore it mus have the same dimensions as X.
#' @param output_prices Matrix or dataframe with prices for each DMU and output.
#' Therefore it mus have the same dimensions as Y.
#' @param RTS Character string indicating the returns-to-scale, e.g. "crs", "vrs"
#' @return A list object containing optimal inputs and outputs and lambdas
#' @examples
#' X <- matrix(c(1,2,3,3,2,1,2,2), ncol = 2)
#' Y <- matrix(c(1,1,1,1), ncol = 1)
#'
#' input_prices <- matrix(c(2,1,2,1,2,1,1,2), ncol =  2, byrow = T)
#' output_prices <- matrix(c(1,1,1,1), ncol = 1)
#'
#' max_prof_lin<- lprofitDEA(X,Y,input_prices,output_prices)
#'
#' @import lpSolveAPI
#' @export


lprofitDEA <- function(X, Y, input_prices, output_prices, RTS = "vrs") {
  
  # Check arguments given by user 
  if (!is.matrix(X) && !is.data.frame(X) && !is.numeric(X)){
    stop("X must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(Y) && !is.data.frame(Y) && !is.numeric(Y)){
    stop("Y must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(input_prices) && !is.data.frame(input_prices) && !is.numeric(input_prices)){
    stop("input_prices must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(output_prices) && !is.data.frame(output_prices) && !is.numeric(output_prices)){
    stop("output_prices must be a numeric vector, matrix or dataframe")
  }
  if (!identical(dim(X), dim(input_prices))) {
    stop("Dimensions of input_prices and X are not the same.")
  }
  if (!identical(dim(Y), dim(output_prices))) {
    stop("Dimensions of output_prices and Y are not the same.")
  }
  
  RTS <- tolower(RTS)
  possible_rts <- c("crs", "vrs")
  if (!(RTS %in% possible_rts)){
    stop("Scale of returns not implemented:", RTS)
  }
  
  
  # Change data structure to uniformly be matrices
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  input_prices <- as.matrix(input_prices)
  output_prices <- as.matrix(output_prices)
  
  # Change data structure to add columns
  in_out_data <- rbind(t(X),t(Y))
  
  # storage for results
  lambdas <- data.frame()
  opt_value <- data.frame()
  
  
  # Solve for all DMUs
  for (i in 1:ncol(in_out_data)) {
    
    # Create the linear programming problem
    # Constraints are nrow and columns are ncol plus
    # the number of rows for optimizing inputs and outputs
    profit_lp <- make.lp(nrow = nrow(in_out_data), ncol = ncol(in_out_data)+nrow(in_out_data))
    
    # Change to maximization problem
    lp.control(profit_lp, sense="max")
    
    # Set the coefficients of the linear programming problem (frontier)
    for (column in 1:ncol(in_out_data)) {
      set.column(profit_lp, column, in_out_data[,column])
    }
    
    # Set decision variables for the objective function
    # optimizing inputs and outputs of the observed DMU
    for (j in 1:nrow(in_out_data)){
      column <- c(rep(0, nrow(in_out_data)))
      column[j] <- -1
      set.column(profit_lp, ncol(in_out_data)+j, column)
    }
    
    # Set the objective function minus for costs and plus for revenues
    # with the respective prices
    set.objfn(profit_lp, c(rep(0,ncol(in_out_data)),-input_prices[i,],output_prices[i,]))
    
    # Set constraint-types
    set.constr.type(profit_lp, c(rep("<=", ncol(X)), rep(">=", ncol(Y))))
    
    # Set the lower bounds of the decision variables
    set.bounds(profit_lp, lower = c(rep(0,(ncol(in_out_data)+nrow(in_out_data)))))
    
    # Set the right-hand side of the constraints
    set.rhs(profit_lp, c(rep(0, nrow(in_out_data))))
    
    # Set the sum of the lambdas to 1
    if (RTS == "vrs") {
      add.constraint(profit_lp, c(rep(1, ncol(in_out_data))), 
                     indices = c(1:ncol(in_out_data)), "=", 1)
      set.bounds(profit_lp, upper = c(rep(1,(ncol(in_out_data)))), columns = c(1:ncol(in_out_data)))
    }
    
    # Solve the linear programming problem
    solve(profit_lp)
    
    # Get the optimal solution and store results in a data frame
    variables <- get.variables(profit_lp)
    lambdas <- rbind(lambdas, variables[1:ncol(in_out_data)])
    opt_value <- rbind(opt_value, variables[(ncol(in_out_data)+1):(ncol(in_out_data)+nrow(in_out_data))])
  }
  
  # Add column names
  colnames(lambdas) <- c(paste("Lambda",1:nrow(X),sep=""))
  colnames(opt_value) <- c(paste("X",1:ncol(X),sep=""), paste("Y",1:ncol(Y),sep=""))
  
  # Return the results
  return(list(lambdas = lambdas, opt_value = opt_value))
}

