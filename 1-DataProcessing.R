# Loading functions
# source("0-Functions.R")
# Loading battery charges and initial power values from Matlab
for (rr.limit in c(1, 5, 7)) {
  
  charge <- R.matlab::readMat(paste0("Data/Data", rr.limit, "/charge", rr.limit, ".mat"))
  power <- R.matlab::readMat(paste0("Data/Data", rr.limit, "/power", rr.limit, ".mat"))
  # Adding together original battery charges and initial power values
  wind.data <- sapply(1:4, function(i) {
    
    sapply(1:100, function(j){
      
      # Getting charges
      C <- charge[[i]][[j]][[1]]
      # Checking charges existence
      if(is.null(C)) return(NULL)
      # Getting power
      P <- t(power[[i]][[j]][[1]])
      n.paths <- nrow(P) # number of paths
      # Taking the absolute value of the charges and adding initial and terminal 0s
      C <- cbind(rep(0, n.paths), t(abs(C))[1:n.paths, , drop = FALSE], rep(0, n.paths))
      # Creating time line. 
      # It's assumed (A1) that the first and last time step is half of a unit (hours)
      t.line <- c(0, (1:j)-0.5, j)
      # Creating the initial and final powers according to assumption (A1), 
      # as well as the power barrier, and Bessel bridges
      # Transitioning from negative to whatever...
      if (i <= 2) {
        IP <- pmax(pmin(P[, 1] + C[, 2] + 0.02 * rr.limit * t.line[2], 2), 0.02 * rr.limit * t.line[j+2])
        FP = pmax(P[, j] + C[, j+1] - 0.02 * rr.limit * (t.line[j+2] - t.line[j+1]), 0)
        P <- cbind(IP, P, FP)
        UB <- P + C # upper power barrier
        if (any(P < 0)) {
          print(paste0("Negative power! Generating a max absolute error of ", round(max(-P[P < 0]), 4),
                       " (",  round(-max(P[P < 0])/rr.limit/0.02*100, 3), "% of the ramp-rate limit)"))
        }
        if (any(P > 2)) {
          print(paste0("Greater than 2 power! Generating a max absolute error of ", round(max(P[P > 2] - 2), 4),
                       " (",  round(max(P[P > 2] - 2)/rr.limit/0.02*100, 3), "% of the ramp-rate limit)"))
        }
        P[P < 0] <- 0; P[P > 2] <- 2
        UB[UB < 0] <- 0; UB[UB > 2] <- 2
        # Creating Bessel bridges
        Br <- UB - P
        Br[P == 0 | P == 2] <- NaN 
        Br[, c(1, j+2)] <- 0
        # Returning list with the data
        return(list(
          charge = C, power = P, boundary = UB, bridges = Br, time = t.line
        ))
        # ... from positive to whatever
      } else {
        IP <- pmin(pmax(P[, 1] - C[, 2] - 0.02 * rr.limit * t.line[2], 0), 2 - 0.02 * rr.limit * t.line[j+2])
        FP = pmin(P[, j] - C[, j+1] + 0.02 * rr.limit * (t.line[j+2] - t.line[j+1]), 2)
        P <- cbind(IP, P, FP)
        LB <- P - C # lower power barrier
        if (any(P < 0)) {
          print(paste0("Negative power! Generating a max absolute error of ", round(max(-P[P < 0]), 4),
                       " (",  round(max(-P[P < 0])/rr.limit/0.02*100, 2), "% of the ramp-rate limit)"))
        }
        if (any(P > 2)) {
          print(paste0("Greater than 2 power! Generating a max absolute error of ", round(max(P[P > 2] - 2), 4),
                       " (",  round(max(P[P > 2] - 2)/rr.limit/0.02*100, 2), "% of the ramp-rate limit = ", rr.limit, ")"))
        }
        P[P < 0] <- 0; P[P > 2] <- 2
        # Creating Bessel bridges
        Br <-  P - LB
        Br[P == 0 | P == 2] <- NaN
        Br[, c(1, j+2)] <- 0
        # Returning list with the data
        return(list(
          charge = C, power = P, boundary = LB, bridge = Br, time = t.line
        ))
      }
      
    })
  }, simplify = FALSE)
  # Reordering and renaming the data in a more intuitive way
  wind.data <- list("-1->0" = wind.data[[2]], "-1->+1" = wind.data[[1]],
                    "+1->-1" = wind.data[[3]], "+1->0" = wind.data[[4]])
  # Saving data
  save(wind.data, file = paste0("Data/Data", rr.limit, "/wind.data", rr.limit, ".RData"))
  
}
