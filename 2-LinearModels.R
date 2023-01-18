# Loading libraries
library(ggplot2)
library(olsrr)
# Loading functions
source("0-Functions.R")
# Running over different ram rate limits
for (rr.limit in c(1, 5, 7)) {
  
  # Loading wind data  
  load(file = paste0("Data/Data", rr.limit, "/wind.data", rr.limit, ".RData"))
  
  # Visualizing and generating data
  lin.Fit <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
  vol.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
  height.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
  t.height.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
  T.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
  power.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
  
  for (from_to in c("-1->0", "-1->+1", "+1->-1", "+1->0")) {
    
    ### DATA PROCESSING
    
    for (i in 2:100) { 
      
      if (!is.null(wind.data[[from_to]][[i]])) {
        
        # Allocating data
        C <- wind.data[[from_to]][[i]]$charge
        P <- wind.data[[from_to]][[i]]$power
        Br <- wind.data[[from_to]][[i]]$bridge
        t.line <- wind.data[[from_to]][[i]]$time
        t.final <- t.line[i+2]
        t.line <- t.line/t.final
        
        ### Parameters estimation
        n.paths <- nrow(C) # sample size
        # Saving total time
        T.total[[from_to]] <- c(T.total[[from_to]], rep(t.final, n.paths)) # saving total heights
        # Saving initial powers
        power.total[[from_to]] <- c(power.total[[from_to]], P[, 1]) # saving total heights
        # Saving maximum heights
        height <- apply(C, 1, max)
        height.total[[from_to]] <- c(height.total[[from_to]], height) # saving total heights
        # Saving the times in which heights are attained
        t.height <- t.line[sapply(1:n.paths, 
                                  function(j) min(which(C[j, ] == height[j])))]
        t.height.total[[from_to]] <- c(t.height.total[[from_to]], t.height) # saving total heights
        
        ### Creating triangles and computing errors
        triangle <- matrix(nrow = n.paths, ncol = i+2)
        errors <- matrix(nrow = n.paths, ncol = i+2)
        for (j in 1:n.paths) {
          triangle[j, ] <- (t.line <= t.height[j]) * (t.line * height[j] / t.height[j]) + 
            (t.line > t.height[j]) * ((1 - t.line) * height[j] / (1 - t.height[j]))
          errors[j, ] <- C[j, ] - triangle[j, ]
        }
        errors.na <- errors
        errors.na[is.na(Br)] <- NaN
        
        # Estimating volatilities for each error path
        vol <- rep(0, n.paths)
        for (j in 1:n.paths) {
          vol[j]  <- vol.mle(t.line = t.line[!is.na(errors.na[j, ])],
                             path = errors[j, !is.na(errors.na[j, ])],
                             t.cut = t.height[j])
        }
        vol.total[[from_to]] <- c(vol.total[[from_to]], vol)
        
        # removing NaNs
        no_nan.idx <- !is.na(vol.total[[from_to]])
        vol.total[[from_to]] <- vol.total[[from_to]][no_nan.idx]
        power.total[[from_to]] <- power.total[[from_to]][no_nan.idx]
        T.total[[from_to]] <- T.total[[from_to]][no_nan.idx]
        height.total[[from_to]] <- height.total[[from_to]][no_nan.idx]
        t.height.total[[from_to]] <- t.height.total[[from_to]][no_nan.idx]
        
        ### Printing information
        # cat("ramp rate = ", rr.limit, ", transition = ", from_to, ", time length = ", t.final, ", number of paths = ", n.paths, "\n", sep = "")
        
      }
      
    }
    
    ### VISUALIZATION
    
    df <- data.frame(vol = vol.total[[from_to]], t.final = T.total[[from_to]],
                     t.height = t.height.total[[from_to]], 
                     height = height.total[[from_to]],
                     power = power.total[[from_to]])
    
    frm <- vol ~ (t.final + t.height + height + power)^2
    linFit <- lm(frm, data = df)
    linFit_trafo <- trafo::trafo_lm(linFit)
    # Removing influential points
    cooksD <- cooks.distance(linFit_trafo$trafo_mod)
    influential <- which((cooksD > (3 * mean(cooksD, na.rm = TRUE))))
    # Fitting linear model without influential points
    df <- df[-influential,]
    linFit <- lm(frm, data = df)
    linFit_trafo <- trafo::trafo_lm(linFit)
    # visualizing transformed linear fit
    layout(matrix(c(1, 2, 3, 4), nrow = 2))
    plot(linFit_trafo$trafo_mod)
    print(summary(linFit_trafo$trafo_mod))
    
    lin.Fit[[from_to]] <- linFit_trafo
    
  }
  
  # Saving linear models
  save(lin.Fit, file = paste0("Data/Data", rr.limit, "/linear.models", rr.limit, ".RData"))
  
}
