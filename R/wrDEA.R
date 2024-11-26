#' @title Estimation of DEA efficiency scores with linear input or output orientation and trade-off weight restrictions
#'
#' @description Linear DEA estimation including the possibility of trade-off weight restrictions,
#' external referencing, and super-efficiency scores. Furthermore, in a second stage slacks can be estimated.
#' The function returns efficiency scores and adjusted lambdas according to the imposed weight restrictions and slack estimation. Additionally, 
#' mus are returned if weight restrictions are imposed that highlight binding restrictions for DMUs and the absolute slack values 
#' if slack-based estimation is applied.
#' 
#' @author Alexander Öttl
#'
#' @param X Vector, matrix or dataframe with DMUs as rows and inputs as columns
#' @param Y Vector, matrix or dataframe with DMUs as rows and outputs as columns
#' @param ORIENTATION Character string indicating the orientation of the DEA model, e.g. "in", "out"
#' @param RTS Character string indicating the returns-to-scale, e.g. "crs", "vrs", "ndrs", "nirs", "fdh"
#' @param WR Matrix with one row per homogeneous linear weight restriction in standard form. The columns are 
#' ncol(WR) = ncol(Y) + ncol(X). Hence the first ncol(Y) columns are the restrictions on outputs and the last ncol(X) columns are the 
#' restrictions on inputs. 
#' @param XREF Vector, matrix or dataframe with firms defining the technology as rows and inputs as columns
#' @param YREF Vector, matrix or dataframe with firms defining the technology as rows and outputs as columns
#' @param SUPEREFF Boolean variable indicating whether super-efficiencies shall be estimated
#' @param SLACK Boolean variable indicating whether slack-based estimation should be applied
#' @param ECONOMIC Boolean variable indicating whether economic efficiency (cost or revenue) is estimated using trade-off weight restrictions. 
#' Therefore mus can be positive and negative when set to TRUE and you have a multi-stage optimization where
#' the first stage is the efficiency estimation and the second stage is the minimization of the mus. You have to
#' use single trade-off, hence a 1 to 1 substitution of either inputs or outputs given by the focus on
#' either cost or revenue efficiency. The data must be value data. 
#'
#' @return A list object containing the following information:
#' \item{eff}{Are the estimated efficiency scores for the DMUs under observation stored 
#' in a vector with the length nrow(X).}
#' \item{lambdas}{Estimated values for the composition of the respective Benchmarks.
#' The lambdas are stored in a matrix with the dimensions nrow(X) x nrow(X), where
#' the row is the DMU under observation and the columns the peers used for the Benchmark. NOTE:
#' Lambdas are automatically slack optimized.}
#' \item{mus}{If WR != NULL, the estimated decision variables for the imposed weight restrictions
#'  are stored in a matrix with the dimensions nrow(X) x nrow(WR), where the rows are the DMUs and 
#'  columns the weight restrictions. If the values are positive, the WR is binding for the respective DMU.}
#' \item{slack}{If SLACK = TRUE, the slacks are estimated and stored in a matrix with the dimensions
#' nrow(X) x (ncol(X) + ncol(Y)). Showing the Slack of each DMU (row) for each input and output
#' (column).}
#'
#' @examples
#' X <- c(1,1,2,4,1.5,2,4,3)
#' Y <- c(1,2,4,4,0.5,2.5,3.5,4)
#' 
#' # Two weight restrictions in standard form first on output then input.
#' # The first WR shows the trade-off that inputs can be reduced by one unit
#' # which reduces outputs by four units. The second WR shows that outputs can 
#' # be increased by one unit when inputs are increased by four units.
#' 
#' WR <- matrix(c(-4,-1,1,4), nrow = 2, byrow = TRUE)
#' 
#' wrDEA(X, Y, ORIENTATION = "in", RTS="vrs", WR = WR)
#' 
#' # For an estimation just focusing on one DMU one can for example use 
#' # XREF and YREF to define the technology and then estimate the efficiency for 
#' # the DMU under observation (here DMU 1). Let's additionally estimate the slacks. 
#' 
#' wrDEA(X[1], Y[1], ORIENTATION = "in", RTS="vrs", XREF = X, YREF = Y, SLACK = TRUE, WR = WR)
#'
#' @import lpSolveAPI
#'
#' @export
#' 

wrDEA <- function(X, Y, ORIENTATION = "out", RTS = "vrs", WR = NULL,
                 XREF = NULL, YREF = NULL, SUPEREFF = FALSE, SLACK = FALSE,
                 ECONOMIC = FALSE) {
  
  # Check arguments given by user 
  check_arguments(X = X, Y = Y, ORIENTATION = ORIENTATION, RTS = RTS, 
                  WR = WR, XREF = XREF, YREF = YREF)

  # Change data structure to uniformly be matrices
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  
  if (!is.null(WR)){
    if (!is.matrix(WR) && !is.data.frame(WR)){
      WR <- t(as.matrix(WR))
    }
  }
  
  if (is.null(XREF)&&is.null(YREF)){
    
    XREF <- X
    YREF <- Y
    
  } else{
    if (!is.matrix(XREF)){
      XREF <- as.matrix(XREF)
    }
    if (!is.matrix(YREF)){
      YREF <- as.matrix(YREF)
    }
  }
  
  ORIENTATION <- tolower(ORIENTATION)
  RTS <- tolower(RTS)

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
    
    # Set decision variables for the objective function last column
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
    if (ECONOMIC) {
      set.bounds(dea_model, lower = c(rep(0,nrow(XREF)),
                                      rep(-Inf, nrow(WR)), 0))
    } else {
      set.bounds(dea_model, lower = c(rep(0,(ncol(in_out_data)+1))))
    }
    
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
    
    # Estimate and store the mus
    if (!is.null(WR)){
      if (ECONOMIC) {
        # Here we have a Multi-Objective Optimization as we try
        # to keep the  absolute values of mu's as small as possible relevant for
        # cost- and revenue efficiency analysis to obtain optimal quantities.
        # First we set the efficiency as fixed from the previous optimization
        add.constraint(dea_model, c(rep(0,ncol(in_out_data)),1), "=", variables[(ncol(in_out_data)+1)])
        
        # Add auxiliary variables for each mu and two constraints to estimate
        # absolute values of mu's (see https://lpsolve.sourceforge.net/5.5/absolute.htm)
        for (j in 1:nrow(WR)) {
          # Add one auxiliary variable y_j for each mu_j
          add.column(dea_model, rep(0, nrow(dea_model)))
        }
        for (j in 1:nrow(WR)) {
          # Constraint 1: mu_j + y_j >= 0
          # Add a row if needed and set coefficients
          add.constraint(dea_model, c(1, 1), 
                         indices = c((nrow(XREF) + j), (ncol(in_out_data)+ 1 + j)), type = ">=", rhs = 0)
          
          # Constraint 2: -mu_j + y_j >= 0
          add.constraint(dea_model, c(-1, 1), 
                         indices = c((nrow(XREF) + j), (ncol(in_out_data)+ 1 + j)), type = ">=", rhs = 0)
        }
        # Set objective function to minimize the sum of all auxiliary variables
        # with negative to adjust for maximization problem in output orientation
        if (ORIENTATION == "in") {
          set.objfn(dea_model, c(rep(0,ncol(in_out_data)), 0,
                                 rep(1, nrow(WR))))
        } else {
          set.objfn(dea_model, c(rep(0,ncol(in_out_data)), 0,
                                 rep(-1, nrow(WR))))
        }
        # Solve second stage optimization
        solve(dea_model)
      }
      # Store the mu's
      variables <- get.variables(dea_model)
      mu <- rbind(mu, variables[(nrow(XREF)+1):(nrow(XREF)+nrow(WR))])
    } else {
      mu <- NULL
    }
    
    # Store the efficiency score
    if (ORIENTATION == "in") {
      eff <- c(eff, variables[(ncol(in_out_data)+1)])
    } else {
      eff <- c(eff, 1/variables[(ncol(in_out_data)+1)])
    }
    
    # Store the lambdas
    lambdas <- rbind(lambdas, variables[1:nrow(XREF)])
  }
  
  if (SLACK){
    # use slack function from other file 
    if (ORIENTATION == "in") {
      slack_est <- slack(X = X, Y = Y, XREF = XREF, YREF = YREF, RTS = RTS, EFF = eff, ALPHA = 0, WR = WR)
    } else {
      slack_est <- slack(X = X, Y = Y, XREF = XREF, YREF = YREF, RTS = RTS, EFF = eff, ALPHA = 1, WR = WR)
    }
    lambdas <- slack_est$lambda
    slack_results <- slack_est$slack
  } else {
    slack_results <- NULL
  }
  
  # Add column names
  colnames(lambdas) <- c(paste("L",1:nrow(XREF),sep=""))
  rownames(lambdas) <- paste("DMU", 1:nrow(X))
  
  if (!is.null(WR)){
    colnames(mu) <- c(paste("MU",1:nrow(WR),sep=""))
  }
  
  # Return the results
  return(list(eff = eff, lambdas = lambdas, mus = mu, slack = slack_results))

}


