slack <- function(X, Y, RTS, ALPHA, EFF, XREF, YREF, WR = NULL){
  
  # if there is no specific XREF YREF, then use the same as X and Y
  # this is automatically given by the functions that call slack
  
  if (!is.null(WR)){
    # first outputs then inputs to fit the WR standard form
    in_out_data <- cbind(rbind(t(YREF), t(XREF)),t(WR))
    
  } else {
    # Change data structure to add columns
    in_out_data <- rbind(t(YREF), t(XREF))
  }
  
  
  for (i in 1:nrow(X)) {
    # Create the linear programming problem
    # Constraints are nrow and columns are ncol plus slack for each input and output
    slack_model <- make.lp(nrow = nrow(in_out_data), ncol = ncol(in_out_data)+nrow(in_out_data))
    
    # Slack is a maximization problem
    lp.control(slack_model, sense="max")
    
    # Set the coefficients of the linear programming problem (frontier)
    for (column in 1:ncol(in_out_data)) {
      set.column(slack_model, column, in_out_data[,column])
    }
    
    # Add decision variable for output slack in separate columns
    for (j in 1:ncol(Y)){
      column <- c(rep(0, nrow(in_out_data)))
      column[j+ncol(X)] <- -1
      set.column(slack_model, nrow(XREF)+j, column)
    }
    
    # Add decision variable for input slack in separate columns
    for (j in 1:ncol(X)){
      column <- c(rep(0, nrow(in_out_data)))
      column[j+ncol(Y)] <- 1
      set.column(slack_model, nrow(XREF)+ncol(Y)+j, column)
    }
    
    # set objective function only slack decision variables
    set.objfn(slack_model, c(rep(0,nrow(XREF)),rep(1,ncol(X)+ncol(Y))))
    
    # Set constraint-types
    set.constr.type(slack_model, c(rep("=", nrow(in_out_data))))
    
    # Set the lower bounds of the decision variables
    set.bounds(slack_model, lower = c(rep(0,(ncol(in_out_data)+nrow(in_out_data)))))
    
    # Set RTS assumption
    if (RTS == "vrs") {
      add.constraint(slack_model, c(rep(1, nrow(XREF))), 
                     indices = c(1:nrow(XREF)), "=", 1)
      set.bounds(slack_model, upper = c(rep(1,(nrow(XREF)))), columns = c(1:nrow(XREF)))
    }
    if (RTS == "ndrs") {
      add.constraint(slack_model, c(rep(1, nrow(XREF))), 
                     indices = c(1:nrow(XREF)), ">=", 1)
    }
    if (RTS == "nirs") {
      add.constraint(slack_model, c(rep(1, nrow(XREF))), 
                     indices = c(1:nrow(XREF)), "<=", 1)
    }
    if (RTS == "fdh") {
      add.constraint(slack_model, c(rep(1, nrow(XREF))), 
                     indices = c(1:nrow(XREF)), "=", 1)
      set.bounds(slack_model, upper = c(rep(1,(nrow(XREF)))), columns = c(1:nrow(XREF)))
      set.type(slack_model, columns = c(1:nrow(XREF)), type = "binary")
    }
    
    # Solve the linear programming problem
    solve(slack_model)
    
    # Get the optimal solution and store results in a data frame
    variables <- get.variables(slack_model)
    if (!is.null(WR)){
      mu <- rbind(mu, variables[(nrow(XREF)+1):(nrow(XREF)+nrow(WR))])
    } else {
      mu <- NULL
    }
    slack <- rbind(slack, variables[(ncol(in_out_data)+1):(ncol(in_out_data)+nrow(in_out_data))])
    lambdas <- rbind(lambdas, variables[1:nrow(XREF)])
  }
  
  # Return the results
  return(list(mu = mu, slack = slack, lambdas = lambdas))
}