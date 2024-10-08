#' @title Hyperbolic estimation of DEA efficiency scores
#'
#' @description Hyperbolic DEA implementation including weight restrictions,
#' non-discretionary variables, gerenralized distance function, external referencing, estimation of slacks and
#' super-efficiency scores. The mathematical and theoretical foundations of the code are presented in the paper
#' "Data Envelopment Analysis and Hyperbolic Efficiency Measures: Extending Applications and Possiblities
#' for Between-Group Comparisons" (2023) by Alexander Öttl, Mette Asmild, and Daniel Gulde.
#'
#' @param X Vector, matrix or dataframe with DMUs as rows and inputs as columns
#' @param Y Vecotr, matrix or dataframe with DMUs as rows and outputs as columns
#' @param RTS Character string indicating the returns-to-scale, e.g. "crs", "vrs", "ndrs", "nirs", "fdh"
#' @param WR Matrix with one row per homogeneous linear weight restriction in standard form. The columns are 
#' ncol(WR) = ncol(Y) + ncol(X). Hence the first ncol(Y) columns are the restrictions on outputs and the last ncol(X) columns are the 
#' restrictions on inputs. 
#' @param SLACK Boolean variable indicating whether an additional estimation of slacks shall be performed when set to 'TRUE'.  
#' Be aware that SLACK estimation can change the lambda values.
#' @param ACCURACY Accuracy value for non-linear programm solver.
#' @param XREF Vector, matrix or dataframe with firms defining the technology as rows and inputs as columns
#' @param YREF Vector, matrix or dataframe with firms defining the technology as rows and outputs as columns
#' @param SUPEREFF Boolean variable indicating whether super-efficiencies shall be estimated
#' @param NONDISC_IN Vector containing column indices of the input matrix that are non-discretionary variables e.g. c(1,3) so the first and the third input are non-discretionary
#' @param NONDISC_OUT Vector containing column indices of the output matrix that are non-discretionary variables e.g. c(1,3) so the first and the third output are non-discretionary
#' @param PARALLEL Integer of amount of cores that should be used for parallel computing (Check availability of computer)
#' @param ALPHA ALPHA can be chosen between [0,1]. It indicates the relative weights given to the distance function to
#' both outputs and inputs when approaching the frontier. More weight on the input orientation is set by alpha < 0.5. Here,
#' the input efficiency score is estimated in the package. To receive the corresponding output efficiency score, estimate: e^((1-alpha)/alpha).
#' Vice versa for an output weighted model alpha > 0.5. The output efficiency is given and the input efficiency can
#' be recovered with: e^(alpha/(1-alpha))
#' 
#'
#' @return A list object containing the following information:
#' \item{eff}{Are the estimated efficiency scores for the DMUs under observation stored 
#' in a vector with the length nrow(X).}
#' \item{lambdas}{Estimated values for the composition of the respective Benchmarks.
#' The lambdas are stored in a matrix with the dimensions nrow(X) x nrow(X), where
#' the row is the DMU under observation and the columns the peers used for the Benchmark.}
#' \item{mus}{If WR != NULL, the estimated decision variables for the imposed weight restrictions
#'  are stored in a matrix with the dimensions nrow(X) x nrow(WR), where the rows are the DMUs and 
#'  columns the weight restrictions. If the values are positive, the WR is binding for the respective DMU.}
#' \item{slack}{If SLACK = TRUE, the slacks are estimated and stored in a matrix with the dimensions
#' nrow(X) x (ncol(X) + ncol(Y)). Showing the Slack of each DMU (row) for each input and output
#' (column).}
#' 
#'
#' @examples
#' X <- c(1,1,2,4,1.5,2,4,3)
#' Y <- c(1,2,4,4,0.5,2.5,3.5,4)
#' # we now impose linked weght restrictions. We assume outputs decrease by 
#' # four units when inputs are reduced by one. And we assume that outputs can
#' # can be increased by one when inputs are increased by four 
#'
#' WR <- matrix(c(-4,-1,1,4), nrow = 2, byrow = TRUE)
#' hyperbolicDEA(X,Y,RTS="vrs", WR = WR)
#'
#' # Another example having the same data but just estimate the results for DMU 1
#' # using XREF YREF and and a higher focus on inputs adjusting the ALPHA towards 0.
#' # Additionally, slacks are estimated.
#'
#' hyperbolicDEA(X[1],Y[1],RTS="vrs", XREF = X, YREF = Y, WR = WR, ALPHA = 0.1, SLACK = TRUE)
#'
#' @import dplyr
#' @import nloptr
#' @import lpSolveAPI
#' @import foreach
#' @import doParallel
#' @import Benchmarking
#'
#' @export

hyperbolicDEA <- function(X, Y, RTS = "vrs", WR = NULL, SLACK=FALSE,
                           ACCURACY = 1.0e-10, XREF = NULL, YREF = NULL,
                           SUPEREFF = FALSE, NONDISC_IN = NULL, NONDISC_OUT = NULL,
                           PARALLEL = 1, ALPHA = 0.5){

  # Check arguments given by user
  check_arguments(X = X, Y = Y, RTS = RTS, WR = WR, XREF = XREF, 
                  YREF = YREF, NONDISC_IN = NONDISC_IN, NONDISC_OUT = NONDISC_OUT,
                  ALPHA = ALPHA)
  
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
  if (!is.null(XREF)&&!is.null(YREF)){
    if (!is.matrix(XREF)){
      XREF <- as.matrix(XREF)
    }
    if (!is.matrix(YREF)){
      YREF <- as.matrix(YREF)
    }
  }  

  RTS <- tolower(RTS)
  

  # Variable for if condition in SLACK estimation and WR,X,Y,XREF,YREF for post scaling
  XREF_YREF <- FALSE
  WR_slack <- WR
  X_slack <- X
  Y_slack <- Y
  ifelse(!is.null(XREF), XREF_slack <- XREF, XREF_slack <- X)
  ifelse(!is.null(YREF), YREF_slack <- YREF, YREF_slack <- Y)

  # scaling adjustments
  # and referring X and Y to XREF and YREF as well as matrix definition
  WR_scaler <- NULL # used to adjust mus again to original scale
  if (is.null(XREF)&&is.null(YREF)){
    if (!is.null(WR)){
      # Scaling weights to observed values - taking mean of all variables that 
      # are in the WR and adjusting with mean of absolute values of WR
      WR_scaler <- mean(cbind(Y,X)[,which(colSums(WR != 0) > 0)])/mean(abs(WR[WR != 0]))
      WR <- t(t(WR*WR_scaler)/c(apply(Y, 2, median),apply(X, 2, median)))
    }
    X <- t(t(X)/apply(X, 2, median))
    Y <- t(t(Y)/apply(Y, 2, median))

    XREF <- X
    YREF <- Y

  } else{
    # equal scaling of XREF YREF and Y and X
    if (!is.null(WR)){
      WR_scaler <- mean(cbind(YREF,XREF)[,which(colSums(WR != 0) > 0)])/mean(abs(WR[WR != 0]))
      WR <- matrix(t(t(WR * WR_scaler)/c(apply(YREF, 2, median),apply(XREF, 2, median))), nrow = nrow(WR))
    }
    X <- t(t(X)/apply(XREF, 2, median))
    Y <- t(t(Y)/apply(YREF, 2, median))
    XREF <- t(t(XREF)/apply(XREF, 2, median))
    YREF <- t(t(YREF)/apply(YREF, 2, median))
  }

  # adjustments to NONDISC variables
  DISC_IN <- 1:ncol(XREF)
  DISC_OUT <- 1:ncol(YREF)
  if (!is.null(NONDISC_IN)){
    DISC_IN <- DISC_IN[-NONDISC_IN]
  }
  if (!is.null(NONDISC_OUT)){
    DISC_OUT <- DISC_OUT[-NONDISC_OUT]
  }

  # Setting up result matrices and vectors
  eff <- c()
  theta <- c()

  lambdas <- matrix(NA, nrow = nrow(X), ncol = nrow(XREF))

  if (!is.null(WR)){
    mus <- matrix(NA, nrow = nrow(X), ncol = nrow(WR))
  } else{
    mus <- NULL
  }



  #register multicore estimation if adjusted
  if (PARALLEL > 1){
    registerDoParallel(PARALLEL)
  }

  # start of the main loop for non-linear solver (suppress warning of
  # parallel backend not registerd if core =1 )
  suppressWarnings(
  result_list <- foreach (i = 1:nrow(X), .packages = "nloptr") %dopar% {

    # Change XREF YREF to check superefficiency
    if (SUPEREFF){
      XREF <- X
      YREF <- Y
      XREF <- as.matrix(XREF[-i,])
      YREF <- as.matrix(YREF[-i,])
    }

    # controls is a vector containing all lambdas, Eff score,
    # and all the mu's (one mu per row of WR)
    # Objective function is first decision variable after lambdas
    eval_f <- function(controls){
      return(controls[nrow(XREF)+1])
    }

    # Setting up optimization without weight restrictions
    if (is.null(WR)){

      # Gradient of the objective function
      eval_grad_f <- function(controls){
        return(c(rep(0, nrow(XREF)), 1))
      }

      # Inequality constraints as vector for non-disc and disc variables
      # as well as ndrs and nirs if specified
      eval_g_ineq <- function(controls){
        constr <- c()
        for (j in c(1:ncol(XREF))[DISC_IN]){
          constraint_x <- controls[1:nrow(XREF)]%*%XREF[,j] - controls[nrow(XREF)+1]^(1-ALPHA)*X[i,j]
          constr <- c(constr, constraint_x)
        }
        for (j in c(1:ncol(XREF))[NONDISC_IN]){
          constraint_x <- controls[1:nrow(XREF)]%*%XREF[,j] - X[i,j]
          constr <- c(constr, constraint_x)
        }
        for (k in c(1:ncol(YREF))[DISC_OUT]){
          constraint_y <- -(controls[1:nrow(XREF)])%*%YREF[,k] + (1/controls[nrow(XREF)+1]^ALPHA)*Y[i,k]
          constr <- c(constr, constraint_y)
        }
        for (k in c(1:ncol(YREF))[NONDISC_OUT]){
          constraint_y <- -(controls[1:nrow(XREF)])%*%YREF[,k] + Y[i,k]
          constr <- c(constr, constraint_y)
        }
        if (RTS == "ndrs"){
          constraint_lambdas <- 1 - sum(controls[1:nrow(XREF)])
          constr <- c(constr, constraint_lambdas)
        }
        if (RTS == "nirs"){
          constraint_lambdas <- sum(controls[1:nrow(XREF)]) - 1
          constr <- c(constr, constraint_lambdas)
        }
        return(constr)
      }

      # Jacobian of the inequality constraints stored in a matrix for
      # each constraint
      eval_jac_g_ineq <- function(controls){
        jacobian_X <- c()
        jacobian_Y <- c()
        for (j in c(1:ncol(XREF))[DISC_IN]){
          jacobian_X_j <- c(XREF[,j],-((1-ALPHA)*controls[nrow(XREF)+1]^(-ALPHA)*X[i,j]))
          jacobian_X <- rbind(jacobian_X, jacobian_X_j)
        }
        for (j in c(1:ncol(XREF))[NONDISC_IN]){
          jacobian_X_j <- c(XREF[,j], 0)
          jacobian_X <- rbind(jacobian_X, jacobian_X_j)
        }
        for (k in c(1:ncol(YREF))[DISC_OUT]){
          jacobian_Y_k <- c(-YREF[,k], -(ALPHA*controls[nrow(XREF)+1]^(-ALPHA-1))*Y[i,k])
          jacobian_Y <- rbind(jacobian_Y, jacobian_Y_k)
        }
        for (k in c(1:ncol(YREF))[NONDISC_OUT]){
          jacobian_Y_k <- c(-YREF[,k], 0)
          jacobian_Y <- rbind(jacobian_Y, jacobian_Y_k)
        }
        jacobian <- rbind(jacobian_X,jacobian_Y)

        if (RTS == "ndrs"){
          jacobian_ndrs <- c(rep(-1, nrow(XREF)), 0)
          jacobian <- rbind(jacobian, jacobian_ndrs)
        }
        if (RTS == "nirs"){
          jacobian_nirs <- c(rep(1, nrow(XREF)), 0)
          jacobian <- rbind(jacobian, jacobian_nirs)
        }

        return(jacobian)
      }

      # Equality constraints and jacobian as vector for vrs
      # otherwise NULL
      if (RTS == "vrs"){
        eval_g_eq <- function(controls){
          constr <- c(sum(controls[1:nrow(XREF)]) - 1)
          return(constr)
        }

        eval_jac_g_eq <- function(controls){
          jacobian <- c(rep(1, nrow(XREF)), 0)
          return(jacobian)
        }
      } else{
        eval_g_eq <- NULL
        eval_jac_g_eq <- NULL
      }

    # Setting up optimization with weight restrictions
    } else{
      # Gradient of the objective function
      eval_grad_f <- function(controls){
        return(c(rep(0, nrow(XREF)), 1, c(rep(0, nrow(WR)))))
      }

      # Inequality constraints as vector for non-disc and disc variables
      # with WR as well as ndrs and nirs if specified
      eval_g_ineq <- function(controls){
        constr <- c()
        for (j in c(1:ncol(XREF))[DISC_IN]){
          constraint_x <- controls[1:nrow(XREF)]%*%XREF[,j] - controls[nrow(XREF)+1]^(1-ALPHA)*X[i,j] + controls[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))]%*%WR[,ncol(Y)+j]
          constr <- c(constr, constraint_x)
        }
        for (j in c(1:ncol(XREF))[NONDISC_IN]){
          constraint_x <- controls[1:nrow(XREF)]%*%XREF[,j] - X[i,j] + controls[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))]%*%WR[,ncol(Y)+j]
          constr <- c(constr, constraint_x)
        }
        for (k in c(1:ncol(YREF))[DISC_OUT]){
          constraint_y <- -(controls[1:nrow(XREF)]%*%YREF[,k]) + (1/controls[nrow(XREF)+1]^ALPHA)*Y[i,k] - (controls[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))])%*%WR[,k]
          constr <- c(constr, constraint_y)
        }
        for (k in c(1:ncol(YREF))[NONDISC_OUT]){
          constraint_y <- -(controls[1:nrow(XREF)]%*%YREF[,k]) + Y[i,k] - (controls[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))])%*%WR[,k]
          constr <- c(constr, constraint_y)
        }
        if (RTS == "ndrs"){
          constraint_lambdas <- 1 - sum(controls[1:nrow(XREF)])
          constr <- c(constr, constraint_lambdas)
        }
        if (RTS == "nirs"){
          constraint_lambdas <- sum(controls[1:nrow(XREF)]) - 1
          constr <- c(constr, constraint_lambdas)
        }
        return(constr)
      }

      # Jacobian of the inequality constraints as matrix
      eval_jac_g_ineq <- function(controls){
        jacobian_X <- c()
        jacobian_Y <- c()
        for (j in c(1:ncol(XREF))[DISC_IN]){
          jacobian_X_j <- c(XREF[,j],-((1-ALPHA)*controls[nrow(XREF)+1]^(-ALPHA)*X[i,j]), WR[,ncol(Y)+j])
          jacobian_X <- rbind(jacobian_X, jacobian_X_j)
        }
        for (j in c(1:ncol(XREF))[NONDISC_IN]){
          jacobian_X_j <- c(XREF[,j], 0, WR[,ncol(Y)+j])
          jacobian_X <- rbind(jacobian_X, jacobian_X_j)
        }
        for (k in c(1:ncol(YREF))[DISC_OUT]){
          jacobian_Y_k <- c(-YREF[,k], -(ALPHA*controls[nrow(XREF)+1]^(-ALPHA-1))*Y[i,k], -WR[,k])
          jacobian_Y <- rbind(jacobian_Y, jacobian_Y_k)
        }
        for (k in c(1:ncol(YREF))[NONDISC_OUT]){
          jacobian_Y_k <- c(-YREF[,k], 0, -WR[,k])
          jacobian_Y <- rbind(jacobian_Y, jacobian_Y_k)
        }

        jacobian <- rbind(jacobian_X,jacobian_Y)

        if (RTS == "ndrs"){
          jacobian_ndrs <- c(rep(-1, nrow(XREF)), 0, rep(0, nrow(WR)))
          jacobian <- rbind(jacobian, jacobian_ndrs)
        }
        if (RTS == "nirs"){
          jacobian_nirs <- c(rep(1, nrow(XREF)), 0, rep(0, nrow(WR)))
          jacobian <- rbind(jacobian, jacobian_nirs)
        }

        return(jacobian)
      }

      # Equality constraints and jacobian as vector for VRS
      if (RTS == "vrs"){
        eval_g_eq <- function(controls){
          constr <- c(sum(controls[1:nrow(XREF)]) - 1)
          return(constr)
        }

        eval_jac_g_eq <- function(controls){
          jacobian <- c(rep(1, nrow(XREF)), 0, rep(0, nrow(WR)))
          return(jacobian)
        }
      } else{
        eval_g_eq <- NULL
        eval_jac_g_eq <- NULL
      }
    }

    # Initial values in controls0 to start optimization
    # and bounds for the optimization
    if (is.null(WR)){
      controls0 <- c(rep(0, nrow(XREF)), 1)
      if (!XREF_YREF){
        controls0[i] <- 1
      }
      lb <- c(rep(0, nrow(XREF)), 0)
      ub <- c(rep(Inf, nrow(XREF)), Inf)
    } else{
      controls0 <- c(rep(0, nrow(XREF)), 1, c(rep(0, nrow(WR))))
      if (!XREF_YREF){
        controls0[i] <- 1
      }
      lb <- c(rep(0, nrow(XREF)), 0, c(rep(0, nrow(WR))))
      ub <- c(rep(Inf, nrow(XREF)), Inf, c(rep(Inf, nrow(WR))))
    }


    # setting number of inequality constraints
    if (RTS == "ndrs" || RTS == "nirs"){
      no_constr_ineq <- ncol(X)+ncol(Y)+1
    } else{
      no_constr_ineq <- ncol(X)+ncol(Y)
    }

    # Set up non-linear optimizer
    opts <-list("algorithm"="NLOPT_LD_SLSQP",
                "xtol_rel"=ACCURACY,
                "maxeval"=1000000,
                "tol_constraints_ineq"=rep(ACCURACY,no_constr_ineq))
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

    # res will be stored in result_list for each i given by the foreach loop

  }
  ) # end of suppress warnings

  # Deregister parallel core
  if (PARALLEL > 1){
    registerDoSEQ()
  }

  # Extract results for weight restrictions
  for (i in 1:nrow(X)){
    if (!is.null(WR)){
      mus[i,] <- result_list[[i]]$solution[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))]
    }

    # fdh estimation
    # solved manually by calculating the rhs of the fdh problem for all possible solutions
    # (iterating over all lambdas being 1 once),
    # calculating the corresponding "hypothetical" eff scores, picking the binding one
    # (highest eff score) for each possible solution and subsequently picking the best (minimum)
    # eff score and the corresponding lambda vector

    if (RTS == "fdh"){
      if (is.null(WR)){
        peer_list <- c()
        eff_list <- c()
        for (j in 1:nrow(XREF)){
          lambdas_fdh <- c(rep(0, nrow(XREF)))
          lambdas_fdh[j] <- 1
          rhs <- c(lambdas_fdh%*%XREF, lambdas_fdh%*%YREF)
          # The following line has a ifelse condition to select only DMUs that
          # are worse of for comparison. Hence, for input efficiency we are only
          # interested in DMUs that have similar or lower output levels. For output
          # efficiency we are only interested in DMUs that have similar or higher
          # input levels. If the conditions is not fullfilled we set the score to 1
          # with efficiency by default. Then we only use the max efficiency score
          # of the DMUs. Hence the output or input the DMU is relatively best at.
          if (ALPHA != 1) {
            eff_X <- ifelse(any(rhs[c((ncol(XREF)+1):(ncol(XREF)+ncol(YREF)))[DISC_OUT]]<Y[i, DISC_OUT]),1,max((rhs[c(1:ncol(XREF))[DISC_IN]]/X[i, DISC_IN])^(1/(1-ALPHA)), na.rm = T))
          } else{
            eff_X <- 0
          }

          if (ALPHA != 0) {
            eff_Y <- ifelse(any(rhs[c(1:ncol(XREF))[DISC_IN]]>X[i, DISC_IN]), 1, max((1/(rhs[c((ncol(XREF)+1):(ncol(XREF)+ncol(YREF)))[DISC_OUT]])*Y[i, DISC_OUT])^(1/ALPHA), na.rm = T))
          } else{
            eff_Y <- 0
          }
          peer_list <- c(peer_list, j)
          # Depending on which efficiency is better output or input efficiency is
          # selected. This depends also on the respective alpha value
          eff_list <- c(eff_list, ifelse(eff_Y>=eff_X, eff_Y^(ALPHA), eff_X^(1-ALPHA)))
          # Theta is stored for slack estimations.
          theta_list <- c(eff_list, ifelse(eff_Y>=eff_X, eff_Y, eff_X))
        }

        # Extract results for fdh optimization
        possible_eff <- data.frame(peer_list, eff_list)
        eff_fdh <- min(possible_eff$eff_list)
        peer <- which(possible_eff[, "eff_list"] == eff_fdh)[1]
        lambdas_fdh <- c(rep(0, nrow(XREF)))
        lambdas_fdh[peer] <- 1
        eff <- c(eff, eff_fdh)
        lambdas[i,] <- lambdas_fdh
        theta <- c(theta, min(theta_list))
      } else{
        stop("FDH cannot be combined with weight restrictions")
      }
    } else{
      # Storage of results if not fdh
      if (SUPEREFF){
        eff <- c(eff, ifelse(ALPHA >= 0.5, result_list[[i]]$solution[nrow(XREF)]^ALPHA,
                             result_list[[i]]$solution[nrow(XREF)]^(1 - ALPHA)))
        lambda <- result_list[[i]]$solution[1:(nrow(XREF)-1)]
        lambda <- append(lambda, NA, after = i-1)
        lambdas[i,] <- lambda
      } else{
        # Storage of results if not fdh or super-efficiency
        eff <- c(eff, ifelse(ALPHA >= 0.5, result_list[[i]]$solution[nrow(XREF)+1]^ALPHA,
                             result_list[[i]]$solution[nrow(XREF)+1]^(1 - ALPHA)))
        lambdas[i,] <- result_list[[i]]$solution[1:nrow(XREF)]
        # For Slack estimation Theta
        theta <- c(theta, result_list[[i]]$solution[nrow(XREF)+1])
      }
    }
  }


  ####### Slack estimation if True ######
  if (SLACK){

    slack_est <- slack(X = X_slack, Y = Y_slack, XREF = XREF_slack, 
                       YREF = YREF_slack, RTS = RTS, 
                       EFF = theta, ALPHA = ALPHA, WR = WR_slack, 
                       NONDISC_IN = NONDISC_IN, NONDISC_OUT = NONDISC_OUT)
    
    lambdas <- slack_est$lambda
    slack_results <- slack_est$slack
  } else{
    slack_results <- NULL
  }
  
  # Renaming results
  if (!is.null(WR)){
    mus <- mus*WR_scaler
    colnames(mus) <- paste("MU", 1:nrow(WR), sep = "")
  }
  
  colnames(lambdas) <- paste("L", 1:nrow(XREF), sep = "")
  rownames(lambdas) <- paste("DMU", 1:nrow(X))

  return(list(eff = eff, lambdas = lambdas, mus = mus, slack = slack_results))
}
