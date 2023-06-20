#' @title hyperbolicDEA
#'
#' @description Hyperbolic DEA implementation including weight restrictions,
#' exogenous variables, external referencing, estimation of slacks and
#' superefficiency scores
#'
#' @param X Matrix or dataframe with DMUS as rows and inputs as columns
#' @param Y Matrix or dataframe with DMUs as rows and outputs as columns
#' @param RTS Character string indicating the returns-to-scale, e.g. "crs", "vrs", "ndrs", "nirs", "fdh"
#' @param WR Matrix with one row per homogenous linear weight restriction in standard form, ncol(WR) = ncol(X) + ncol(Y)
#' @param SLACK Variable indicating whether an additional estimation of slacks shall be performed
#' @param ACCURACY Accuracy value for non-linear programmer
#' @param XREF Matrix or dataframe with firms defining the technology as rows and inputs as columns
#' @param YREF Matrix or dataframe with firms defining the technology as rows and outputs as columns
#' @param SUPEREFF Variable indicating whether superefficiencies shall be estimated
#' @param NONDISC_IN Vector containing indices of the input matrix that are non-discretionary variables
#' @param NONDISC_OUT Vector containing indices of the output matrix that are non-discretionary variables
#' @param PARALLEL Integer of amount of cores that should be used for estimation (Check availability of computer)
#'
#' @return A list object containing efficiency scores, lambdas, and potentially slacks and
#' binding parameters in the weight restrictions (mus)
#' @examples
#' X <- c(1,1,2,4,1.5,2,4,3)
#' Y <- c(1,2,4,4,0.5,2.5,3.5,4)
#' hyperbolicDEA(X,Y,RTS="vrs", SUPEREFF = F)
#' @export hyperbolicDEA
#' @import dplyr
#' @import nloptr
#' @import lpSolveAPI
#' @import foreach
#' @import doParallel

hyperbolicDEA <- function(X, Y, RTS = "vrs", WR = NULL, SLACK=F,
                           ACCURACY = 1.0e-10, XREF = NULL, YREF = NULL,
                           SUPEREFF = F, NONDISC_IN = NULL, NONDISC_OUT = NULL,
                           PARALLEL = 1){

  if (!is.matrix(X) && !is.data.frame(X) && !is.numeric(X)){
    stop("X must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }
  if (!is.matrix(Y) && !is.data.frame(Y) && !is.numeric(Y)){
    stop("Y must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  if (!is.null(WR)){
    if (ncol(WR) != ncol(X) + ncol(Y)){
      stop("WR must be a matrix of weight restrictions in standard form,
           ncol(WR) = ncol(Y) + ncol(X)")
    }
  }

  possible_rts <- c("crs", "vrs", "ndrs", "nirs", "fdh")

  # Variable for if condition in SLACK estimation
  XREF_YREF <- F

  # scaling adjustments
  # and referring X and Y to XREF and YREF as well as matrix definition
  if (is.null(XREF)&&is.null(YREF)){

    scaled_values <- scale(rbind(cbind(Y,X),WR), center = F)
    Y <- as.matrix(scaled_values[1:nrow(X),1:ncol(Y)])
    X <- as.matrix(scaled_values[1:nrow(X),(ncol(Y)+1):(ncol(X)+ncol(Y))])
    if (!is.null(WR)){
      WR <- matrix(scaled_values[(nrow(X)+1):(nrow(X)+nrow(WR)),1:(ncol(X)+ncol(Y))],
                   ncol = (ncol(X)+ncol(Y)))
    }

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

    # equal scaling of XREF YREF and Y and X
    XREF_YREF <- T
    scaled_values <- scale(rbind(cbind(rbind(Y,YREF),rbind(X,XREF)), WR) , center = F)
    Y <- as.matrix(scaled_values[1:nrow(X),1:ncol(Y)])
    X <-  as.matrix(scaled_values[1:nrow(X),(ncol(Y)+1):(ncol(X)+ncol(Y))])
    YREF <- as.matrix(scaled_values[(nrow(X)+1):(nrow(X)+nrow(XREF)),1:ncol(Y)])
    XREF <- as.matrix(scaled_values[(nrow(X)+1):(nrow(X)+nrow(XREF)),(ncol(Y)+1):(ncol(X)+ncol(Y))])
    if (!is.null(WR)){
      WR <- matrix(scaled_values[(nrow(X)+nrow(XREF)+1):(nrow(X)+nrow(XREF)+nrow(WR)),
                                 1:(ncol(X)+ncol(Y))], ncol = (ncol(X)+ncol(Y)))
    }
  }

  if (!(RTS %in% possible_rts)){
    stop("Unknown scale of returns:", RTS)
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

  lambdas <- matrix(NA, nrow = nrow(X), ncol = nrow(XREF))
  colnames(lambdas) <- paste("L", 1:nrow(XREF), sep = "")
  rownames(lambdas) <- paste("DMU", 1:nrow(X))

  if (is.null(WR)){
    mus <- NULL
  } else{
    mus <- matrix(NA, nrow = nrow(X), ncol = nrow(WR))
    colnames(mus) <- paste("MU", 1:nrow(WR), sep = "")
  }

  results <- list(eff = eff, lambdas = lambdas, mus = mus)

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

    # controls is a vector containing all lambdas, G (Eff score),
    # and all the mu's (one mu per row of WR)
    eval_f <- function(controls){
      return(controls[nrow(XREF)+1])
    }

    if (is.null(WR)){

      eval_grad_f <- function(controls){
        return(c(rep(0, nrow(XREF)), 1))
      }

      eval_g_ineq <- function(controls){
        constr <- c()
        for (j in c(1:ncol(XREF))[DISC_IN]){
          constraint_x <- controls[1:nrow(XREF)]%*%XREF[,j] - controls[nrow(XREF)+1]*X[i,j]
          constr <- c(constr, constraint_x)
        }
        for (j in c(1:ncol(XREF))[NONDISC_IN]){
          constraint_x <- controls[1:nrow(XREF)]%*%XREF[,j] - X[i,j]
          constr <- c(constr, constraint_x)
        }
        for (k in c(1:ncol(YREF))[DISC_OUT]){
          constraint_y <- -(controls[1:nrow(XREF)])%*%YREF[,k] + (1/controls[nrow(XREF)+1])*Y[i,k]
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

      eval_jac_g_ineq <- function(controls){
        jacobian_X <- c()
        jacobian_Y <- c()
        for (j in c(1:ncol(XREF))[DISC_IN]){
          jacobian_X_j <- c(XREF[,j],-X[i,j])
          jacobian_X <- rbind(jacobian_X, jacobian_X_j)
        }
        for (j in c(1:ncol(XREF))[NONDISC_IN]){
          jacobian_X_j <- c(XREF[,j], 0)
          jacobian_X <- rbind(jacobian_X, jacobian_X_j)
        }
        for (k in c(1:ncol(YREF))[DISC_OUT]){
          jacobian_Y_k <- c(-YREF[,k], -(controls[nrow(XREF)+1]^(-2))*Y[i,k])
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

    } else{
      eval_grad_f <- function(controls){
        return(c(rep(0, nrow(XREF)), 1, c(rep(0, nrow(WR)))))
      }

      eval_g_ineq <- function(controls){
        constr <- c()
        for (j in c(1:ncol(XREF))[DISC_IN]){
          constraint_x <- controls[1:nrow(XREF)]%*%XREF[,j] - controls[nrow(XREF)+1]*X[i,j] + controls[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))]%*%WR[,ncol(Y)+j]
          constr <- c(constr, constraint_x)
        }
        for (j in c(1:ncol(XREF))[NONDISC_IN]){
          constraint_x <- controls[1:nrow(XREF)]%*%XREF[,j] - X[i,j] + controls[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))]%*%WR[,ncol(Y)+j]
          constr <- c(constr, constraint_x)
        }
        for (k in c(1:ncol(YREF))[DISC_OUT]){
          constraint_y <- -(controls[1:nrow(XREF)]%*%YREF[,k]) + (1/controls[nrow(XREF)+1])*Y[i,k] - (controls[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))])%*%WR[,k]
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

      eval_jac_g_ineq <- function(controls){
        jacobian_X <- c()
        jacobian_Y <- c()
        for (j in c(1:ncol(XREF))[DISC_IN]){
          jacobian_X_j <- c(XREF[,j],-X[i,j], WR[,ncol(Y)+j])
          jacobian_X <- rbind(jacobian_X, jacobian_X_j)
        }
        for (j in c(1:ncol(XREF))[NONDISC_IN]){
          jacobian_X_j <- c(XREF[,j], 0, WR[,ncol(Y)+j])
          jacobian_X <- rbind(jacobian_X, jacobian_X_j)
        }
        for (k in c(1:ncol(YREF))[DISC_OUT]){
          jacobian_Y_k <- c(-YREF[,k], -(controls[nrow(XREF)+1]^(-2))*Y[i,k], -WR[,k])
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


    if (RTS == "ndrs" || RTS == "nirs"){
      no_constr_ineq <- ncol(X)+ncol(Y)+1
    } else{
      no_constr_ineq <- ncol(X)+ncol(Y)
    }


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

  }
  ) # end of suppress warnings

  # Deregister parallel core
  if (PARALLEL > 1){
    registerDoSEQ()
  }

  for (i in 1:nrow(X)){
    if (!is.null(WR)){
      results$mus[i,] <- result_list[[i]]$solution[(nrow(XREF)+2):(nrow(XREF)+1+nrow(WR))]
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
          lambdas <- c(rep(0, nrow(XREF)))
          lambdas[j] <- 1
          rhs <- c(lambdas%*%XREF, lambdas%*%YREF)
          eff_X <- rhs[c(1:ncol(XREF))[DISC_IN]]/X[i,DISC_IN]
          eff_Y <- 1/(rhs[c((ncol(XREF)+1):(ncol(XREF)+ncol(YREF)))[DISC_OUT]])*Y[i, DISC_OUT]
          peer_list <- c(peer_list, j)
          eff_list <- c(eff_list, max(c(eff_X, eff_Y)))
        }
        possible_eff <- data.frame(peer_list, eff_list)
        eff_fdh <- min(possible_eff$eff_list)
        peer <- possible_eff[possible_eff$eff_list == eff_fdh, "peer_list"]
        lambdas <- c(rep(0, nrow(XREF)))
        lambdas[peer] <- 1
        eff <- c(eff, eff_fdh)
        results$lambdas[i,] <- lambdas
      } else{
        stop("FDH cannot be combined with weight restrictions")
      }
    } else{
      if (SUPEREFF){
        eff <- c(eff, result_list[[i]]$solution[nrow(XREF)])
        lambda <- result_list[[i]]$solution[1:(nrow(XREF)-1)]
        lambda <- append(lambda, NA, after = i-1)
        results$lambdas[i,] <- lambda
      } else{
        eff <- c(eff, result_list[[i]]$solution[nrow(XREF)+1])
        results$lambdas[i,] <- result_list[[i]]$solution[1:nrow(XREF)]
      }
    }
  }

  results$eff <- eff

  ####### Slack estimation if True ######
  if (SLACK){

    # unscale for slack estimation
    if (XREF_YREF){
      non_scaled_values <- t(apply(scaled_values, 1, function(r)r*attr(scaled_values,'scaled:scale')))
      X <- as.matrix(non_scaled_values[1:nrow(X),(ncol(X)+1):(ncol(X)+ncol(Y))])
      Y <-  as.matrix(non_scaled_values[1:nrow(X),1:ncol(X)])
      XREF <- as.matrix(non_scaled_values[(nrow(X)+1):(nrow(X)+nrow(XREF)),(ncol(X)+1):(ncol(X)+ncol(Y))])
      YREF <- as.matrix(non_scaled_values[(nrow(X)+1):(nrow(X)+nrow(XREF)),1:ncol(X)])

    } else{
      # YREF and XREF are the same as X and Y if not specified
      non_scaled_values <- t(apply(scaled_values, 1, function(r)r*attr(scaled_values,'scaled:scale')))
      X <- as.matrix(non_scaled_values[,(ncol(X)+1):(ncol(X)+ncol(Y))])
      Y <-  as.matrix(non_scaled_values[,1:ncol(X)])
      XREF <- as.matrix(non_scaled_values[,(ncol(X)+1):(ncol(X)+ncol(Y))])
      YREF <-  as.matrix(non_scaled_values[,1:ncol(X)])
    }

    results_slack <- c()

    for (i in 1:nrow(X)){
      # set up the problem ( number of constraints, number of decision variables)
      lprec <- make.lp((ncol(X)+ncol(Y)+1),(nrow(XREF)+ncol(X)+ncol(Y)))

      # Set up first part with lambdas and inputs outputs as well as 1 for
      # the lambda constraint in VRS (sum lambdas = 1)
      for (k in 1:nrow(XREF)){
        column <- set.column(lprec, k, c(XREF[k,],YREF[k,],1))

        if (RTS == "fdh"){
          set.type(lprec, k, "binary")
        }
      }

      # Add decision variable for input slack in separate columns
      for (j in 1:ncol(X)){
        column <- c(rep(0, ncol(X)+ncol(Y)+1))
        column[j] <- 1
        set.column(lprec, nrow(XREF)+j, column)
      }

      # Add decision variable for output slack in separate columns
      for (j in 1:ncol(Y)){
        column <- c(rep(0, ncol(X)+ncol(Y)+1))
        column[j+ncol(X)] <- -1
        set.column(lprec, nrow(XREF)+ncol(X)+j, column)
      }

      # Change to maximazation problem
      lp.control(lprec, sense="max")

      # set objective function only slack decision variables
      set.objfn(lprec, c(rep(0,nrow(XREF)),rep(1,ncol(X)+ncol(Y))))

      # set the constraint type
      set.constr.type(lprec, c(rep("=", ncol(X)+ncol(Y)+1)))

      # set bounds 0 to 1 for lambdas and 0 to inf for slack
      set.bounds(lprec, lower = c(rep(0,nrow(XREF)+ncol(X)+ncol(Y))), columns = c(1:(nrow(XREF)+ncol(X)+ncol(Y))))
      set.bounds(lprec, upper = c(rep(1,nrow(XREF))), columns = c(1:nrow(XREF)))
      set.bounds(lprec, upper = c(rep(Inf,ncol(X)+ncol(Y))), columns = c((nrow(XREF)+1):(nrow(XREF)+ncol(X)+ncol(Y))))


      # Solve the problem
      set.rhs(lprec, c(X[i,]*results$eff[i],Y[i,]/results$eff[i],1))
      solve(lprec)
      slacks <- get.variables(lprec)
      results_slack <- rbind(results_slack,slacks)
    }

    # slack
    results$slack <- results_slack[,(nrow(XREF)+1):(nrow(XREF)+ncol(X)+ncol(Y))]
    rownames(results$slack) <- paste("DMU",1:nrow(results$slack))
    colnames(results$slack) <- c(paste("Sx", 1:ncol(X), sep = ""),
                                 paste("Sy", 1:ncol(Y), sep = ""))

    # New peers due to slack in lambdas
    results$lambdas <- results_slack[,1:nrow(XREF)]
    rownames(results$lambdas) <- paste("DMU",1:nrow(results$lambdas))
    colnames(results$lambdas) <- paste("L", 1:ncol(results$lambdas), sep = "")
  }

  return(results)
}
