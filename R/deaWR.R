deaWR<- function(X, Y, WR = NULL, ORIENTATION = "out", RTS = "vrs") {
  
  # Check arguments given by user 
  if (!is.matrix(X) && !is.data.frame(X) && !is.numeric(X)){
    stop("X must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(Y) && !is.data.frame(Y) && !is.numeric(Y)){
    stop("Y must be a numeric vector, matrix or dataframe")
  }

  # Change data structure to uniformly be matrices
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  if (!is.null(WR)){
    if (!is.matrix(WR) && !is.data.frame(WR) && !is.numeric(WR)){
      stop("WR must be a numeric matrix or dataframe")
    }
    if (ncol(WR) != ncol(X) + ncol(Y)){
      stop("WR must be a matrix of weight restrictions in standard form,
           ncol(WR) = ncol(Y) + ncol(X)")
    }
    WR <- as.matrix(WR)
  }

  RTS <- tolower(RTS)
  possible_rts <- c("crs", "vrs")
  if (!(RTS %in% possible_rts)){
    stop("Scale of returns not implemented:", RTS)
  }

  
  # storage for results
  lambdas <- data.frame()
  eff <- c()
  
  if (!is.null(WR)){
    mu <- data.frame()
    
    # Change data structure to add columns
    # first outputs then inputs to fit the WR standard form
    in_out_data <- cbind(rbind(t(Y), t(X)),t(WR))
    
  } else {
    # Change data structure to add columns
    in_out_data <- rbind(t(Y), t(X))
  }
  
  
  # Solve for all DMUs
  for (i in 1:nrow(X)) {
    
    # Create the linear programming problem
    # Constraints are nrow and columns are ncol plus 1 for efficiency score
    dea_model <- make.lp(nrow = nrow(in_out_data), ncol = ncol(in_out_data)+1)
    
    # Change to maximization problem
    if (ORIENTATION == "in") {
      lp.control(dea_model, sense="min")
    } else {
      lp.control(dea_model, sense="max")
    }  
    
    # Set the coefficients of the linear programming problem (frontier)
    for (column in 1:ncol(in_out_data)) {
      set.column(dea_model, column, in_out_data[,column])
    }
    
    # Set decision variables for the objective function last row
    if (ORIENTATION == "in") {
      objective <- c(rep(0, ncol(Y)), -in_out_data[(ncol(Y)+1):(ncol(X)+ncol(Y)),i])
    } else {
      objective <- c(-in_out_data[1:ncol(Y),i],rep(0, ncol(X)))
    }
    
    set.column(dea_model, ncol(in_out_data)+1, objective)
    
    # Set the objective function just one for the last decision variable
    # which is the efficiency score
    set.objfn(dea_model, c(rep(0,ncol(in_out_data)),1))
    
    # Set constraint-types
    set.constr.type(dea_model, c(rep(">=", ncol(Y)), rep("<=", ncol(X))))
    
    # Set the lower bounds of the decision variables
    set.bounds(dea_model, lower = c(rep(0,(ncol(in_out_data)+1))))
    
    # Set the right-hand side of the constraints
    if (ORIENTATION == "in") {
      set.rhs(dea_model, c(in_out_data[1:ncol(Y),i], rep(0, ncol(X))))
    } else {
      set.rhs(dea_model, c(rep(0, ncol(Y)), in_out_data[(ncol(Y)+1):(ncol(X)+ncol(Y)),i]))
    }  
    
    # Set the sum of the lambdas to 1
    if (RTS == "vrs") {
      add.constraint(dea_model, c(rep(1, nrow(X))), 
                     indices = c(1:nrow(X)), "=", 1)
      set.bounds(dea_model, upper = c(rep(1,(nrow(X)))), columns = c(1:nrow(X)))
    }
    
    # Solve the linear programming problem
    solve(dea_model)
    
    # Get the optimal solution and store results in a data frame
    variables <- get.variables(dea_model)
    if (!is.null(WR)){
      lambdas <- rbind(lambdas, variables[1:nrow(X)])
      mu <- rbind(mu, variables[(nrow(X)+1):(nrow(X)+nrow(WR))])
      eff <- c(eff, variables[(nrow(X)+nrow(WR)+1)])
    } else {
      lambdas <- rbind(lambdas, variables[1:ncol(in_out_data)])
      eff <- c(eff, variables[(ncol(in_out_data)+1)])

    }
  }
  
  # Add column names
  colnames(lambdas) <- c(paste("Lambda",1:nrow(X),sep=""))
  
  if (!is.null(WR)){
    colnames(mu) <- c(paste("WR",1:nrow(WR),sep=""))
  }
  
  # Return the results
  if (!is.null(WR)){
    return(list(lambdas = lambdas, mu = mu, eff = eff))
  } else{
    return(list(lambdas = lambdas, eff = eff))
  }
}


