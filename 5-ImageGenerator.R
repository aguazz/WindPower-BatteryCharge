# Loading libraries
library(latex2exp)
library(scales)
library(ks)
# Loading data
load(file = "Data/Data1/linear.models1.RData")
load(file = "Data/Data1/wind.data.sim1.RData")
# Loading functions
source(file = "0-Functions.R")
# Settings
color <- c("firebrick3", "deepskyblue3", "darkolivegreen4")
FigGen <- FALSE
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Figure: Battery charges as the difference between power and corrected power ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Setting ramp-rate limit
rr.limit <- 1
# Loading wind data  
load(file = paste0("Data/Data", rr.limit, "/wind.data", rr.limit, ".RData"))
# Setting transition and Soujorn Time
for (from_to in c("-1->0", "+1->0")) {
i <- 13
# Getting data
C <- wind.data[[from_to]][[i]]$charge[1:3,]
P <- wind.data[[from_to]][[i]]$power[1:3,]
B <- wind.data[[from_to]][[i]]$boundary[1:3,]
t.line <- wind.data[[from_to]][[i]]$time
t.final <- t.line[i+2]
# Printing information
cat("ramp rate = ", rr.limit, ", transition = ", from_to, 
    ", time length = ", t.final, ", number of paths = ", nrow(B), "\n", sep = "")
### BEGIN PLOT
  # Plot's settings   
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/power-boundary_", from_to, "-", i, ".pdf"), 
        width = 8, height = 4)
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }
  
  # Power + corrected power
  matplot(t.line , t(P), type = "l", lty = 1, lwd = 2, xlab = "", ylab = "", 
          ylim = range(B, P), col = color)
  matlines(t.line , t(B), lty = 2, lwd = 2, col = color)

  if (FigGen) dev.off()

}
### END PLOT
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Figure: Real vs synthetic data. RampRate = 1 ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Fixing seed
set.seed(71)
# Setting ramp-rate limit
rr.limit <- 1
# Loading wind data  
load(file = paste0("Data/Data", rr.limit, "/wind.data", rr.limit, ".RData"))
# Setting transition and Soujorn Time
from_to <- "-1->0"
i <- 13
# Getting data
C <- wind.data[[from_to]][[i]]$charge
P <- wind.data[[from_to]][[i]]$power
t.line <- wind.data[[from_to]][[i]]$time
t.final <- t.line[i+2]
t.line <- t.line/t.final
Br <- wind.data[[from_to]][[i]]$bridge
# Parameters estimation
n.paths <- nrow(Br) # sample size
height <- apply(C, 1, max) # maximum heights
t.height <- t.line[sapply(1:n.paths, 
                          function(j) min(which(C[j, ] == height[j])))] # times in which heights are attained
# Printing information
cat("ramp rate = ", rr.limit, ", transition = ", from_to, 
    ", time length = ", t.final, ", number of paths = ", n.paths, "\n", sep = "")
# Creating triangles
triangle <- matrix(nrow = n.paths, ncol = i+2)
for (j in 1:n.paths) {
  triangle[j, ] <- (t.line <= t.height[j]) * (t.line * height[j] / t.height[j]) + 
    (t.line > t.height[j]) * ((1 - t.line) * height[j] / (1 - t.height[j]))
}
# Computing errors
E <- C - triangle
### BEGIN PLOT
# Plot's settings  
if(FigGen) { # Setting for saving the PDF image
  pdf(paste0("Img/charge_", from_to, "-", i, ".pdf"), 
      width = , height = )
  par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
      mgp = c(2, 1.6, 0), tcl = -0.8) 
} else {
  par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
      mgp = c(2, 0.6, 0), tcl = -0.4) 
}

# Battery charges
matplot(t.line*t.final , t(C[1:3, ]), type = "l", lty = 1, lwd = 2, xlab = "", ylab = "", col = color)
lines(c(0, t.final), c(0, 0), lty = 3, lwd = 1)

if (FigGen) dev.off()

###
# 2
###

if(FigGen) { # Setting for saving the PDF image
  pdf(paste0("Img/triangle_", from_to, "-", i, ".pdf"), 
      width = , height = )
  par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
      mgp = c(2, 1.6, 0), tcl = -0.8) 
  # par(mar=c(3, 3, 1, 1), cex.axis = 1.5, cex.lab = 2,
  #     mgp = c(2.5, 1.6, 0), fin = c(6, 6.5), tcl = -0.8) # Setting for printing in RStudio
} else {
  par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
      mgp = c(2, 0.6, 0), tcl = -0.4) 
}

# Simulated discretized Charges
matplot(t.line*t.final , t(triangle[1:3,]), type = "l", lty = 1, lwd = 2, xlab = "", ylab = "", col = color)
lines(c(0, t.final), c(0, 0), lty = 3, lwd = 1)

if (FigGen) dev.off()

###
# 3
###

if(FigGen) { # Setting for saving the PDF image
  pdf(paste0("Img/error_", from_to, "-", i, ".pdf"), 
      width = , height = )
  par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
      mgp = c(2, 1.6, 0), tcl = -0.8) 
} else {
  par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
      mgp = c(2, 0.6, 0), tcl = -0.4) 
}

# Simulated Charges
matplot(t.line*t.final , t(E[1:3,]), type = "l", lty = 1, lwd = 2, xlab = "", ylab = "", col = color)
lines(c(0, t.final), c(0, 0), lty = 3, lwd = 1)

if (FigGen) dev.off()
### END PLOT
# Computing power support
power.upper <- ifelse(from_to %in% c("-1->0", "-1->+1"), 2, 2 - 0.02 * rr.limit * t.final)
power.lower <- ifelse(from_to %in% c("-1->0", "-1->+1"), 0.02 * rr.limit * t.final, 0)
# Adjusting power
if(from_to %in% c("-1->0", "-1->+1")){
  power.adjusted <- pmin(P[, 1] + runif(n.paths, 0, 0.01), power.upper)
} else {
  power.adjusted <- pmax(P[, 1] - runif(n.paths, 0, 0.01), power.lower)
}
M <- 200 # length of the time and space partition for the kd estimation
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
# Simulating t.heights, heights and initial powers
m <-  max(n.paths, 50) # simulation size
# Generating height and t.height
sim <- reject.sample.3d(n = m, pdf.xyz = k.den$w, x = k.den$x, y = k.den$y, z = k.den$z)
t.height.sim <- sim[, 1]
height.sim <- sim[, 2]
power.sim <- sim[, 3]
# Generating Brownian bridges
N <- 200 # continuous version of the time interval
t.line.full <- seq(0, 1, l = N)
BB <- matrix(nrow = m, ncol = N)  # setting BB matrix
vol.sim <- c() # Predicting volatility
linModel <- linFitList[[from_to]]$trafo_mod
lambda <- linFitList[[from_to]]$lambdahat
for (j in 1:m) {
  vol.sim <- c(vol.sim, predict(linModel, newdata = data.frame(x_ = height.sim[j])))
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
# Creating simulated triangles
triangle.sim <- matrix(nrow = m, ncol = N)
for (j in 1:m) {
  triangle.sim[j, ] <- (t.line.full <= t.height.sim[j]) * (t.line.full * height.sim[j] / t.height.sim[j]) + 
                       (t.line.full > t.height.sim[j]) * ((1 - t.line.full) * height.sim[j] / (1 - t.height.sim[j]))
}
# Adding noise to the triangles to create the final charges and computing discretized charges to match same t.line
C.sim <- abs(triangle.sim + BB)
C.sim.discrete <- matrix(nrow = m, ncol = i)
# Bounding each charge path by the power limit
for (j in 1:m) {
  C.sim[j, ] <- pmin(C.sim[j, ], 
                     ifelse(from_to %in% c("-1->0", "-1->+1"), power.sim[j], 2 - power.sim[j]) - 
                     0.02 * rr.limit * t.line.full * t.final)
  C.sim.discrete[j, ] <- abs(splinefun(t.line.full, C.sim[j, ])(t.line[2:(i+1)]))
}
### BEGIN PLOT
  # y-axis limits
  y.max <- max(C, C.sim)
  # Plot's settings  
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/real-data_", from_to, "-", i, ".pdf"), 
        width = , height = )
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
    # par(mar=c(3, 3, 1, 1), cex.axis = 1.5, cex.lab = 2,
    #     mgp = c(2.5, 1.6, 0), fin = c(6, 6.5), tcl = -0.8) # Setting for printing in RStudio
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }

  # Battery charges
  matplot(t.line , t(C), type = "l", lwd = 2, xlab = "", ylab = "", 
          ylim = c(0, y.max))
  lines(c(0, 1), c(0, 0), lty = 3, lwd = 1)
          
  if (FigGen) dev.off()
  
  ###
  # 2
  ###
  
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/simulated-data-discrete_", from_to, "-", i, ".pdf"), 
        width = , height = )
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
    # par(mar=c(3, 3, 1, 1), cex.axis = 1.5, cex.lab = 2,
    #     mgp = c(2.5, 1.6, 0), fin = c(6, 6.5), tcl = -0.8) # Setting for printing in RStudio
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }

  # Simulated discretized Charges
  matplot(t.line , t(cbind(rep(0, m), C.sim.discrete, rep(0, m))), 
          type = "l", lwd = 2, xlab = "", ylab = "", ylim = c(0, y.max))
  lines(c(0, 1), c(0, 0), lty = 3, lwd = 1)
  
  if (FigGen) dev.off()
  
  ###
  # 3
  ###
  
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/simulated-data-continuous_", from_to, "-", i, ".pdf"), 
        width = , height = )
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
    # par(mar=c(3, 3, 1, 1), cex.axis = 1.5, cex.lab = 2,
    #     mgp = c(2.5, 1.6, 0), fin = c(6, 6.5), tcl = -0.8) # Setting for printing in RStudio
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }
  
  # Simulated Charges
  matplot(t.line.full , t(C.sim), type = "l", lwd = 2, xlab = "", ylab = "", 
          ylim = c(0, y.max))
  lines(c(0, 1), c(0, 0), lty = 3, lwd = 1)
  
  if (FigGen) dev.off()
### END PLOT
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Figure: Linear models. RampRate = 1 ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
for (from_to in c("-1->0", "-1->+1", "+1->-1", "+1->0")) {
  ### BEGIN PLOT
  
  ###
  # 1: Fitted vs Residuals (Original)
  ###
  
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/fit-res_O_", from_to, ".pdf"), width = 8, height = 4)
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }
  
  fit <- lin.Fit[[from_to]]$orig_mod$fitted.values
  res <- lin.Fit[[from_to]]$orig_mod$residuals
  
  plot(fit, res, lwd = 1, xlab = "", ylab = "")
  lines(range(fit), c(0, 0), lty = 1, lwd = 3, col = color[1])
  
  if (FigGen) dev.off()
 
  ###
  # 2: Q-Q plot (Original)
  ###
  
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/q-q_O_", from_to, ".pdf"), width = 8, height = 4)
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }
  
  res <- lin.Fit[[from_to]]$orig_mod$residuals
  
  qqnorm(res, lwd = 1, xlab = "", ylab = "", main = "")
  qqline(res, lty = 1, lwd = 3, col = color[1])
  
  if (FigGen) dev.off()
  
  ###
  # 3: Fitted vs Residuals (Transformed)
  ###
  
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/fit-res_T_", from_to, ".pdf"), width = 8, height = 4)
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }
  
  fit <- lin.Fit[[from_to]]$trafo_mod$fitted.values
  res <- lin.Fit[[from_to]]$trafo_mod$residuals
  
  plot(fit, res, lwd = 1, xlab = "", ylab = "")
  lines(range(fit), c(0, 0), lty = 1, lwd = 3, col = color[1])
  
  if (FigGen) dev.off()
  
  ###
  # 4: Q-Q plot (Transformed)
  ###
  
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/q-q_T_", from_to, ".pdf"), width = 8, height = 4)
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }
  
  res <- lin.Fit[[from_to]]$trafo_mod$residuals
  
  qqnorm(res, lwd = 1, xlab = "", ylab = "", main = "")
  qqline(res, lty = 1, lwd = 3, col = color[1])
  
  if (FigGen) dev.off()
  
  
  ### END PLOT 
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Figure: Volatility vs height. RampRate = 1 ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
for (from_to in c("-1->0", "-1->+1", "+1->-1", "+1->0")) {
  ### BEGIN PLOT
  
  ###
  # 1: Fitted vs Residuals (Original)
  ###
  
  if(FigGen) { # Setting for saving the PDF image
    pdf(paste0("Img/vol-vs-height_O_", from_to, ".pdf"), width = 8, height = 4)
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 1.6, 0), tcl = -0.8) 
  } else {
    par(mar = c(3, 3, 1, 1), cex.axis = 1, cex.lab = 1,
        mgp = c(2, 0.6, 0), tcl = -0.4) 
  }
  
  vol <- lin.Fit[[from_to]]$orig_mod$model$vol
  height <- lin.Fit[[from_to]]$orig_mod$model$height
  
  plot(height, vol, lwd = 1, xlab = "", ylab = "")
  
  if (FigGen) dev.off()
  
  ### END PLOT 
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Figure: R Squared. RampRate = 1 ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
for (r in c(1, 5, 7)) {
  # Loading data
  load(file = paste0("Data/Data", r, "/linear.models", r, ".RData"))
  
  for (from_to in c("-1->0", "-1->+1", "+1->-1", "+1->0")) {
    
     adj.r.squared <- summary(lin.Fit[[from_to]]$trafo_mod)$adj.r.squared
     cat("RampRate: ", r, " FromTo: ", from_to,
         " Adjusted R Squared: ", adj.r.squared, "\n")
    
  }
  
}
