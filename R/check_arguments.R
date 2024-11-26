.onAttach <- function(libname, pkg) {
  packageStartupMessage("Welcome to hyperbolicDEA")
  packageStartupMessage("Version: ", utils::packageVersion(pkg))
  packageStartupMessage("Visit https://github.com/AlexanderOe/hyperbolicDEA for more information or opening issues. You can also contact me via alexander@ifro.ku.dk.")
  packageStartupMessage(" ")
  packageStartupMessage("When using the package, please cite the following paper: \u00d6ttl, A., Asmild, M., & Gulde, D. (2023). Data Envelopment Analysis and hyperbolic efficiency measures: Extending applications and possibilities for between-group comparisons.")
}


check_arguments <- function(X, Y, XREF = NULL, YREF = NULL, 
                            WR = NULL, RTS = "vrs", 
                            NONDISC_IN = NULL, NONDISC_OUT = NULL, 
                            ALPHA = NULL, ORIENTATION = NULL){
  
  # General checks for all functions in the package
  if (!is.matrix(X) && !is.data.frame(X) && !is.numeric(X)){
    stop("X must be a numeric vector, matrix or dataframe")
  }
  if (!is.matrix(Y) && !is.data.frame(Y) && !is.numeric(Y)){
    stop("Y must be a numeric vector, matrix or dataframe")
  } 
  
  if (!is.null(WR)){
    if (!is.matrix(WR) && !is.data.frame(WR)){
      WR <- t(as.matrix(WR))
    } 
    if (!is.null(XREF)&&!is.null(YREF)){
      if (ncol(WR) != ncol(as.matrix(XREF)) + ncol(as.matrix(YREF))){
        stop("WR must be a matrix of weight restrictions in standard form,
           ncol(WR) = ncol(Y) + ncol(X)")
      }
    } else {
      if (ncol(WR) != ncol(as.matrix(X)) + ncol(as.matrix(Y))){
        stop("WR must be a matrix of weight restrictions in standard form,
           ncol(WR) = ncol(Y) + ncol(X)")
      }
    }
  }
  if (!is.null(XREF)&&!is.null(YREF)){
    if (!is.matrix(XREF) && !is.data.frame(XREF) && !is.numeric(XREF)){
      stop("XREF must be a numeric vector, matrix or dataframe")
    } 
    if (!is.matrix(YREF) && !is.data.frame(YREF) && !is.numeric(YREF)){
      stop("YREF must be a numeric vector, matrix or dataframe")
    } 
    if (!is.numeric(X) || !is.numeric(Y)){
      if ((ncol(as.matrix(YREF))+ncol(as.matrix(XREF))) != (ncol(as.matrix(X)) + ncol(as.matrix(Y)))){
        stop("XREF and YREF must be the same input-output combination:
           ncol(XREF) = ncol(X); ncol(YREF) = ncol(Y)")
      }
    }
  }  
  if (anyNA(X) || anyNA(Y) || anyNA(XREF) || anyNA(YREF) || anyNA(WR)) {
    stop("The optimizer cannot handle missing or NaN values.")
  }
  
  possible_rts <- c("crs", "vrs", "ndrs", "nirs", "fdh")
  RTS <- tolower(RTS)
  if (!(RTS %in% possible_rts)){
    stop("Unknown scale of returns:", RTS)
  }
  if (!is.null(ALPHA)){
    if (ALPHA < 0 || ALPHA > 1){
      warning("ALPHA should be a number between 0 and 1 otherwise inputs can be increased
            or outputs can be decreased")
    }
  }
  if (!is.null(ORIENTATION)){
    ORIENTATION <- tolower(ORIENTATION)
    if (!(ORIENTATION %in% c("in", "out"))){
      stop("ORIENTATION must be either 'in' or 'out'")
    }
  }
  if (!is.null(NONDISC_IN)){
    if (!(all(NONDISC_IN %in% 1:ncol(as.matrix(X))))){
      stop("NONDISC_IN must be an available index of the input matrix")
    } 
  }
  if (!is.null(NONDISC_OUT)){
    if (!(all(NONDISC_OUT %in% 1:ncol(as.matrix(Y))))){
      stop("NONDISC_OUT must be an available index of the output matrix")
    } 
  }
  
}