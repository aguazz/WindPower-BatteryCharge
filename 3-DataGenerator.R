# Loading functions
source("0-Functions.R")

# Loading libraries
library(ks)

# Running over different ram rate limits
for (rr.limit in c(1, 5, 7)) {

# Loading data  
load(file = paste0("Data/Data", rr.limit, "/linear.models", rr.limit, ".RData"))
load(file = paste0("Data/Data", rr.limit, "/wind.data", rr.limit, ".RData"))

# Initializing simulated data
data.sim <- list("-1->0" = lapply(1:100, function(x) list()), 
                 "-1->+1" = lapply(1:100, function(x) list()),
                 "+1->-1" = lapply(1:100, function(x) list()),
                 "+1->0" = lapply(1:100, function(x) list()))

# # Visualizing and generating data
# vol.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
# height.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
# t.height.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())
# power.total <- list("-1->0" = c(), "-1->+1" = c(), "+1->-1" = c(), "+1->0" = c())

for (from_to in c("-1->0", "-1->+1", "+1->-1", "+1->0")) {
# for (from_to in c("+1->-1", "+1->0")) {
  for (i in 1:100) { 

    if (!is.null(wind.data[[from_to]][[i]])) {
      
      # Allocating data
      C <- wind.data[[from_to]][[i]]$charge
      P <- wind.data[[from_to]][[i]]$power
      B <- wind.data[[from_to]][[i]]$boundary
      t.line <- wind.data[[from_to]][[i]]$time
      t.final <- t.line[i+2]
      t.line <- t.line/t.final
      Br <- wind.data[[from_to]][[i]]$bridge
      
      ### Parameters estimation
      n.paths <- nrow(Br) # sample size
      # Saving maximum heights
      height <- apply(C, 1, max)
      # height.total[[from_to]] <- c(height.total[[from_to]], height) # saving total heights
      # Saving the times in which heights are attained
      t.height <- t.line[sapply(1:n.paths, 
                                function(j) min(which(C[j, ] == height[j])))]
      # t.height.total[[from_to]] <- c(t.height.total[[from_to]], t.height) # saving total heights
      
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
      # Errors magnitude
      error.var <- rowSums(errors.na^2, na.rm = TRUE)
      
      # # Estimating volatilities for each error path
      # vol <- rep(0, n.paths)
      # for (j in 1:n.paths) {
      #   vol[j]  <- vol.mle(t.line = t.line[!is.na(errors.na[j, ])], 
      #                      path = errors[j, !is.na(errors.na[j, ])],
      #                      t.cut = t.height[j])
      # }
      # vol.total[[from_to]] <- c(vol.total[[from_to]], vol)
      
      ### Computing kernel density of (t.height, height)
      # if((length(unique(height)) == 1) && (length(unique(t.height)) > 1)) {
      #   bw.h <- 0.15
      #   bw.t <- 0
      #   q1 <- 0.25; q2 <- 0.75
      #   while(bw.t == 0){
      #     bw.t <- bandwidth.nrd_(t.height, q1, q2)
      #     q1 <- max(q1 * 0.9, 0) 
      #     q2 <- min(q2 * 1.1, 1)
      #   }
      # } else if((length(unique(height)) > 1) && (length(unique(t.height)) == 1)) {
      #   bw.t <- 0.15
      #   bw.h <- 0
      #   q1 <- 0.25; q2 <- 0.75
      #   while(bw.h == 0){
      #     bw.h <- bandwidth.nrd_(height, q1, q2)
      #     q1 <- max(q1 * 0.9, 0) 
      #     q2 <- min(q2 * 1.1, 1)
      #   }
      # } else if((length(unique(height)) == 1) && (length(unique(t.height)) == 1)) {
      #   bw.h <- 0.15
      #   bw.t <- 0.15
      # } else {
      # bw.h <- 0; bw.t <- 0
      # q1 <- 0.25; q2 <- 0.75
      # while((bw.h == 0) | (bw.t == 0)){
      #   bw.h <- bandwidth.nrd_(height, q1, q2)
      #   bw.t <- bandwidth.nrd_(t.height, q1, q2)
      #   q1 <- max(q1 * 0.9, 0) 
      #   q2 <- min(q2 * 1.1, 1)
      # }
      # }
      
      ### Printing information
      cat("ramp rate = ", rr.limit, ", transition = ", from_to, ", time length = ", t.final, ", number of paths = ", n.paths, "\n", sep = "")
      
      # Computing power support
      power.upper <- ifelse(from_to %in% c("-1->0", "-1->+1"), 
                            2, 2 - 0.02 * rr.limit * t.final)
      power.lower <- ifelse(from_to %in% c("-1->0", "-1->+1"), 
                            0.02 * rr.limit * t.final, 0)
      
      # Adjusting power
      # power.total[[from_to]] <- c(power.total[[from_to]], P[, 1]) # saving total heights
      if(from_to %in% c("-1->0", "-1->+1")){
        power.adjusted <- pmin(P[, 1] + runif(n.paths, 0, 0.01), power.upper)
      } else {
        power.adjusted <- pmax(P[, 1] - runif(n.paths, 0, 0.01), power.lower)
      }
      
      M <- 200 # length of the time and space partition
      if(length(unique(height)) <= 3 || length(unique(t.height)) <= 3){
        
        H <- matrix(c(1e-3, 1e-4, 1e-4,
                      1e-4, 1e-3, 1e-4,
                      1e-4, 1e-4, 1e-3), ncol = 3)
        k.den <- ks::kde.boundary(x = cbind(t.height, height, power.adjusted), 
                                  xmin = c(0.01, 0, power.lower),
                                  xmax = c(0.99, 2, power.upper),
                                  gridsize = c(M, M, M), H = H)

      } else {
        
        k.den <- ks::kde.boundary(x = cbind(t.height, height, power.adjusted), 
                                  xmin = c(0.01, 0, power.lower), 
                                  xmax = c(0.99, 2, power.upper),
                                  gridsize = c(M, M, M)) #, H = diag(c(bw.t, bw.h)))  
        
      }
      
      k.den <- list(w = k.den$estimate, x = k.den$eval.points[[1]], 
                    y = k.den$eval.points[[2]], z = k.den$eval.points[[3]])
      # k.den <- MASS::kde2d(x = t.height, y = height, lims = c(0.01, 0.99, 0.001, 1.999), 
      #                      n = c(M, M), h = c(bw.t, bw.h))
      # Adjusting kernel density
      # k.den$z <- k.den$z / (sum(outer(diff(k.den$x), diff(k.den$y), FUN = "*") * 
      #                             (k.den$z[-1, -1] + 
      #                             k.den$z[-length(k.den$x), -length(k.den$x)])/2))
      
      ### Computing kernel density and of initial powers
      
      # Estimating power density
      
      # if(n.paths > 1){
      #   # kd.power <- ks::kde.boundary(x = power.adjusted, 
      #   #                              xmin = power.lower, xmax = power.upper,
      #   #                              boundary.supp = c(power.lower, power.upper),
      #   #                              gridsize = M)
      #   kd.power <- kdensity::kdensity(x = P[, 1], kernel = "gamma", adjust = 5,
      #                                  support = c(power.lower, power.upper))
      # } else {
      #   # kd.power <- ks::kde.boundary(x = power.adjusted, xmin = power.lower, xmax = power.upper,
      #   #                              boundary.supp = c(power.lower, power.upper),
      #   #                              gridsize = M, h = 0.005)
      #   kd.power <- kdensity::kdensity(x = P[, 1], bw = 0.005, kernel = "gamma",
      #                                  support = c(power.lower, power.upper))
      # }
      # # kd.power <- list(x = kd.power$eval.points, y = kd.power$estimate)
      # power.line <- seq(power.lower, power.upper, l = M)
      # kd.power <- list(x = power.line, y = kd.power(power.line))
      # kd.power$y <- kd.power$y / sum(diff(kd.power$x)*(kd.power$y[-1]+kd.power$y[-M])/2)
      # dist.power <- cumsum(diff(kd.power$x)*(kd.power$y[-1]+kd.power$y[-M])/2)
      
      ### Simulating initial powers
      # power.sim <- inverse.dist.sim(m, kd.power$x, dist.power)
      
      ### Simulating t.heights, heights and initial powers
      m <-  max(3*n.paths, 100) # simulation size
        # Adjusting 2d density to be bounded ramp rate slope
        # k.den.power <- adjust.2d.kernel(pdf.xy = k.den$z, x = k.den$x, y = k.den$y,
        #                                 l = ifelse(from_to %in% c("-1->0", "-1->+1"), 
        #                                            power.sim[j], 2-power.sim[j]),
        #                                 r = rr.limit, t.final = t.final)
        # k.den.power <- list(z = k.den.power, x = k.den$x, y = k.den$y)
        # # BEGIN MINI-PLOT
        # # Plot's settings
        # par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        #     mgp = c(2, 0.6, 0), tcl = -0.4)
        # layout(matrix(c(1), nrow = 1, byrow = FALSE), heights = c(1))
        # # kde of t.heights + heights
        # contour(k.den, xlab = "t.height", ylab = "heights",
        #         levels  =  c(0.1, 0.9),
        #         ylim = c(0, max(height)+0.2))
        # par(new=TRUE)
        # contour(k.den.power, xlab = "t.height", ylab = "heights",
        #         levels  =  c(0.1, 0.9), col = "red",
        #         ylim = c(0, max(height)+0.2))
        # points(t.height, height, lwd = 2)
        # lines(k.den.power$x, ifelse(from_to %in% c("-1->0", "-1->+1"),
        #                             power.sim[j], 2-power.sim[j]) - 0.02 * rr.limit * k.den.power$x * t.final,
        #       lty = 2, lwd = 2, col = "red")
        # # Make system wait to visualize the images
        # Sys.sleep(1)
        # # END MINI-PLOT
        
        # Generating height and t.height
        sim <- reject.sample.3d(n = m, pdf.xyz = k.den$w, 
                                       x = k.den$x, y = k.den$y,
                                       z = k.den$z)
        t.height.sim <- sim[, 1]
        height.sim <- sim[, 2]
        power.sim <- sim[, 3]
      
      ### Generating Brownian bridges
      # Continuous version of the time interval
      N <- 200
      t.line.full <- seq(0, 1, l = N)
      # Setting BB matrix
      BB <- matrix(nrow = m, ncol = N)
      # Predicting volatility
      vol.sim <- c()
      linModel <- lin.Fit[[from_to]]$trafo_mod
      lambda <- lin.Fit[[from_to]]$lambdahat
      for (j in 1:m) {
        
        vol.sim <- c(vol.sim, predict(linModel, 
                                      newdata = data.frame(t.final = t.final,
                                                           t.height = t.height.sim[j],
                                                           height = height.sim[j],
                                                           power = power.sim[j])))
        vol.sim[j] <- (lambda*vol.sim[j] + 1)^(1/lambda)
          
      }
      # Simulating BB paths
      for (j in 1:m) {
        
        t1 <- t.line.full[t.line.full <= t.height.sim[j]]
        t2 <- t.line.full[t.line.full > t.height.sim[j]]
        t2 <- c(t1[length(t1)], t2)
        
        dW1 <- vol.sim[j] * rnorm(length(t1)-1) * sqrt(diff(t1))
        W1 <-  c(0, cumsum(dW1))
        BB1 <- W1 - t1/t1[length(t1)] * W1[length(W1)]
        
        dW2 <- vol.sim[j] * rnorm(length(t2)-1) * sqrt(diff(t2))
        W2 <-  cumsum(dW2)
        BB2 <- W2 - (t2[-1] - t2[1])/(1-t2[1]) * W2[length(W2)]
        
        BB[j, ] <- c(BB1, BB2)
        
      }
      
      ### Creating simulated triangles
      triangle.sim <- matrix(nrow = m, ncol = N)
      for (j in 1:m) {
        triangle.sim[j, ] <- (t.line.full <= t.height.sim[j]) * (t.line.full * height.sim[j] / t.height.sim[j]) + 
                         (t.line.full > t.height.sim[j]) * ((1 - t.line.full) * height.sim[j] / (1 - t.height.sim[j]))
      }
      
      ### Adding noise to the triangles to create the final charges
      ### and computing discertized charges to match same t.line
      C.sim <- abs(triangle.sim + BB)
      C.sim.discrete <- matrix(nrow = m, ncol = i)
      
      ### Bounding each charge path by the power limit
      for (j in 1:m) {
          
        C.sim[j, ] <- pmin(C.sim[j, ], 
                           ifelse(from_to %in% c("-1->0", "-1->+1"), power.sim[j], 2 - power.sim[j]) - 
                           0.02 * rr.limit * t.line.full * t.final)
        
        C.sim.discrete[j, ] <- abs(splinefun(t.line.full, C.sim[j, ])(t.line[2:(i+1)]))
          
      }
      
      ### Storing simulated discretized  charges
      data.sim[[from_to]][[i]] <- C.sim.discrete
      
      ### Goodness-of-fit metrics
      
      ### BEGIN PLOT
      # Plot's settings
      par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
          mgp = c(2, 0.6, 0), tcl = -0.4) 
      layout(matrix(c(1, 2, 3), nrow = 1, byrow = FALSE), heights = c(1))
      
      # Battery charges
      matplot(t.line , t(C), type = "l", lwd = 2, xlab = "Time", ylab = "Charges")
      lines(c(0, 1), c(0, 0), lty = 3, lwd = 1)
      
      # Simulated discretized Charges
      matplot(t.line , t(cbind(rep(0, m), C.sim.discrete, rep(0, m))), type = "l", lwd = 2, xlab = "Time", ylab = "Simulated charges")
      lines(c(0, 1), c(0, 0), lty = 3, lwd = 1)
      
      # Simulated Charges
      matplot(t.line.full , t(C.sim), type = "l", lwd = 2, xlab = "Time", ylab = "Simulated charges")
      lines(c(0, 1), c(0, 0), lty = 3, lwd = 1)
      
      # # errors NaN
      # matplot(t.line , t(errors.na), type = "l", lwd = 2, xlab = "Time", ylab = "Errors")
      # lines(c(0, 1), c(0, 0), lty = 3, lwd = 1)
      # 
      # # errors brides
      # matplot(t.line.full , t(BB), type = "l", lwd = 2, xlab = "Time", ylab = "Simulated Errors")
      # lines(c(0, 1), c(0, 0), lty = 3, lwd = 1)
      
      # # kernel density of power
      # plot(kd.power$x, kd.power$y, type = "l", lwd = 2,
      #      xlab = "Initial power", ylab = "KDensity", xaxt="n")
      # axis(side = 1, at = c(power.lower, power.sim, power.upper), 
      #      labels = c(power.lower, rep("", m), power.upper))
      # lines(c(power.lower, power.upper), c(0, 0), lty = 3, lwd = 1)
      # lines(c(power.lower, power.lower), c(0, max(kd.power$y)), lty = 3, lwd = 1)
      
      # # power vs height
      # library(RColorBrewer)
      # cols <- brewer.pal(4, "Blues")
      # pal <- colorRampPalette(c("blue", "red"))
      # pal <- colorRampPalette(cols)
      # color.pal <- findInterval(t.height.total[[from_to]], 
      #                      sort(t.height.total[[from_to]]))
      # 
      # plot(power.total[[from_to]], height.total[[from_to]], lwd = 2, 
      #      xlab = "Initial power", ylab = "Height", 
      #      col = pal(length(color.pal))[color.pal])
      
      # # kde of t.heights + heights
      # contour(k.den, xlab = "t.height", ylab = "heights",
      #         levels  =  c(0.1, 0.9),
      #         ylim = c(0, max(height)+0.2))
      # points(t.height, height, lwd = 2)
      
      # Make system wait to visualize the images
      #Sys.sleep(5)
      ### END PLOT
      
    }
      
  }
  
}

save(data.sim, file = paste0("Data/Data", rr.limit, "/wind.data.sim", rr.limit, ".RData"))

}
