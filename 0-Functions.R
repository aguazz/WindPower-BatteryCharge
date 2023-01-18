# Functions
reject.sample.3d <- function(n, pdf.xyz, x, y, z)
{
  smpl <- data.frame(x = numeric(n), y = numeric(n), z = numeric(n))
  max.val <- max(pdf.xyz)
  i <- 0
  
  while (i < n){
    
    x.val <- sample(x, size = 1)
    x.t <- which(x.val == x)
    
    y.val <- sample(y, size = 1)
    y.t <- which(y.val == y)
    
    z.val <- sample(z, size = 1)
    z.t <- which(z.val == z)
    
    if (runif(1) < pdf.xyz[x.t, y.t, z.t] / max.val){
      
      i <- i + 1
      smpl[i, ] <- c(x.val, y.val, z.val)
      
      cat("simulating path ", i, " out of ", n, "\r", sep = "")
      
    }
    
  }
  
  return(smpl)
  
}
bandwidth.nrd_ <- function(x, q1 = 0.25, q2 = 0.75)
{
  r <- quantile(x, c(q1, q2))
  h <- (r[2] - r[1])/1.34
  4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
}
inverse.dist.sim <- function(n, pdf.x, pdf.y){
  
  pdf.x[findInterval(runif(n), pdf.y)+1] 
  
}
adjust.2d.kernel <- function(pdf.xy, x, y, x.lim, y.lim, l, r, t.final){
  
  idx <- outer(x, y, FUN = function(x, y) y < l - 0.02 * r * x * t.final)
  
  c <- sum(
    (idx[-length(x), -length(y)] + idx[-1, -1])/2 *
      outer(diff(x), diff(y), FUN = "*") * 
      (pdf.xy[-1, -1] + pdf.xy[-length(x), -length(y)])/2
  )
  
  pdf.xy <- idx * pdf.xy / c 
  
  return(pdf.xy)
  
}
vol.mle <- function(t.line, path, t.cut){
  
  t1 <- t.line[t.line <= t.cut]
  path1 <- path[t.line <= t.cut]
  n1 <- length(t1)
  t1.f <- t.cut
  t2 <- t.line[t.line >= t.cut] - t.cut
  path2 <- path[t.line >= t.cut]
  n2 <- length(t2)
  t2.f <- t2[n2]
  
  s1 <- 0
  s2 <- 0
  n <- 0
  
  if(n1 > 2){
    
    for (i in 2:(n1-1)) {
      
      
      s1 <- s1 + (path1[i] - path1[i-1]*(t1.f - t1[i]) / (t1.f - t1[i-1]))^2 /
        ((t1[i] - t1[i-1]) * (t1.f - t1[i]) / (t1.f - t1[i-1]))
      
    }
    
    n <- n + (n1-2)
    
  } 
  
  if(n2 > 2){
    
    for (i in 2:(n2-1)) {
      
      s2 <- s2 + (path2[i] - path2[i-1]*(t2.f - t2[i]) / (t2.f - t2[i-1]))^2 /
        ((t2[i] - t2[i-1]) * (t2.f - t2[i]) / (t2.f - t2[i-1]))
      
    }
    
    n <- n + (n2-2)
    
  }
  
  if(n1 > 2 || n2 > 2){
    
    return(
      sqrt((s1 + s2)/n)
    )
    
  } else {
    
    return(NaN)
    
  }
  
}
