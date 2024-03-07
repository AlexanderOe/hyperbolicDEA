#' @title Estimation of DEA efficiency scores with linear input or output orientation
#'
#' @description Linear DEA estimation including the possibility of trade-off weight restrictions,
#' external referencing, and super-efficiency scores. 
#'
#' @param X Matrix or dataframe with DMUs as rows and inputs as columns
#' @param Y Matrix or dataframe with DMUs as rows and outputs as columns
#' @param ORIENTATION Character string indicating the orientation of the DEA model, e.g. "in", "out"
#' @param RTS Character string indicating the returns-to-scale, e.g. "crs", "vrs", "ndrs", "nirs", "fdh"
#' @param WR Matrix with one row per homogeneous linear weight restriction in standard form. The columns are 
#' ncol(WR) = ncol(Y) + ncol(X). Hence the first ncol(Y) columns are the restrictions on outputs and the last ncol(X) columns are the 
#' restrictions on inputs. 
#' @param XREF Matrix or dataframe with firms defining the technology as rows and inputs as columns
#' @param YREF Matrix or dataframe with firms defining the technology as rows and outputs as columns
#' @param SUPEREFF Boolean variable indicating whether super-efficiencies shall be estimated
#' 
#'
#' @return A list object containing the following information:
#' \item{eff}{Are the estimated efficiency scores for the DMUs under observation stored 
#' in a vector with the length nrow(X).}
#' \item{lambdas}{Estimated values for the composition of the respective Benchmarks.
#' The lambdas are stored in a matrix with the dimensions nrow(X) x nrow(X), where
#' the row is the DMU under observation and the columns the peers used for the Benchmark.}
#' \item{mu}{If WR != NULL, the estimated decision variables for the imposed weight restrictions
#'  are stored in a matrix with the dimensions nrow(X) x nrow(WR), where the rows are the DMUs and 
#'  columns the weight restrictions. If the values are positive, the WR is binding for the respective DMU.}
#' 
#'
#' @examples
#' X <- c(1,1,2,4,1.5,2,4,3)
#' Y <- c(1,2,4,4,0.5,2.5,3.5,4)
#' deaWR(X,Y,RTS="vrs", SUPEREFF = FALSE)
#'
#' @import lpSolveAPI
#'
#' @export

wrDEA<- function(X, Y, ORIENTATION = "out", RTS = "vrs", WR = NULL,
                 XREF = NULL, YREF = NULL, SUPEREFF = FALSE) {
  
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
  
  if (is.null(XREF)&&is.null(YREF)){
    
    XREF <- X
    YREF <- Y
    
  } else{
    if (!is.matrix(XREF) && !is.data.frame(XREF) && !is.numeric(XREF)){
      stop("XREF must be a numeric vector, matrix or dataframe")
    }
    if (!is.matrix(XREF)){
      XREF <- as.matrix(XREF)
    }
    if (!is.matrix(YREF) && !is.data.frame(YREF) && !is.numeric(YREF)){
      stop("YREF must be a numeric vector, matrix or dataframe")
    }
    if (!is.matrix(YREF)){
      YREF <- as.matrix(YREF)
    }
    if ((ncol(as.matrix(YREF))+ncol(as.matrix(XREF))) != (ncol(as.matrix(X)) + ncol(as.matrix(Y)))){
      stop("XREF and YREF must be the same input-output combination:
           ncol(XREF) = ncol(X); ncol(YREF) = ncol(Y)")
    }
  }
  
  RTS <- tolower(RTS)
  possible_rts <- c("crs", "vrs", "ndrs", "nirs", "fdh")
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
    in_out_data <- cbind(rbind(t(YREF), t(XREF)),t(WR))
    
  } else {
    # Change data structure to add columns
    in_out_data <- rbind(t(YREF), t(XREF))
  }
  if (SUPEREFF){
     supereff_dat <- in_out_data
  }
  
  
  # Solve for all DMUs
  for (i in 1:nrow(X)) {
    
    # Change data for superefficiency
    if (SUPEREFF){
      in_out_data <- supereff_dat
      in_out_data[,i] <- 0
    }
    
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
      objective <- c(rep(0, ncol(Y)), -X[i,])
    } else {
      objective <- c(-Y[i,], rep(0, ncol(X)))
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
      set.rhs(dea_model, c(Y[i,], rep(0, ncol(X))))
    } else {
      set.rhs(dea_model, c(rep(0, ncol(Y)), X[i,]))
    }  
    
    # Set RTS assumption
    if (RTS == "vrs") {
      add.constraint(dea_model, c(rep(1, nrow(XREF))), 
                     indices = c(1:nrow(XREF)), "=", 1)
      set.bounds(dea_model, upper = c(rep(1,(nrow(XREF)))), columns = c(1:nrow(XREF)))
    }
    if (RTS == "ndrs") {
      add.constraint(dea_model, c(rep(1, nrow(XREF))), 
                     indices = c(1:nrow(XREF)), ">=", 1)
    }
    if (RTS == "nirs") {
      add.constraint(dea_model, c(rep(1, nrow(XREF))), 
                     indices = c(1:nrow(XREF)), "<=", 1)
    }
    if (RTS == "fdh") {
      add.constraint(dea_model, c(rep(1, nrow(XREF))), 
                     indices = c(1:nrow(XREF)), "=", 1)
      set.bounds(dea_model, upper = c(rep(1,(nrow(XREF)))), columns = c(1:nrow(XREF)))
      set.type(dea_model, columns = c(1:nrow(XREF)), type = "binary")
    }
    
    # Solve the linear programming problem
    solve(dea_model)
    
    # Get the optimal solution and store results in a data frame
    variables <- get.variables(dea_model)
    if (!is.null(WR)){
      mu <- rbind(mu, variables[(nrow(XREF)+1):(nrow(XREF)+nrow(WR))])
      eff <- c(eff, variables[(nrow(XREF)+nrow(WR)+1)])
    } else {
      mu <- NULL
      eff <- c(eff, variables[(ncol(in_out_data)+1)])
    }
    lambdas <- rbind(lambdas, variables[1:nrow(XREF)])
  }
  
  # Add column names
  colnames(lambdas) <- c(paste("L",1:nrow(XREF),sep=""))
  
  if (!is.null(WR)){
    colnames(mu) <- c(paste("MU",1:nrow(WR),sep=""))
  }
  
  # Return the results
  return(list(lambda = lambdas, mus = mu, eff = eff))

}


