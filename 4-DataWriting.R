# Setting transition labels
transition <- c("negative-zero", "negative-positive", "positive-negative", "positive-zero")
# Running over different ram rate limits
for (rr.limit in c(1, 5, 7)) {
  
  # Loading data
  load(file = paste0("Data/Data", rr.limit, "/wind.data.sim", rr.limit, ".RData"))
  
  for (from_to in c("-1->0", "-1->+1", "+1->-1", "+1->0")) {
    
    for (i in 1:100) { 
      
      cat("Sojourn time", i, "\n")
      if (!is.null(data.sim[[from_to]][[i]])) {
        
        charges <- data.sim[[from_to]][[i]]
        # Saving charge as csv file
        write.csv(charges, file = paste0("Data/Data", rr.limit, "/CSV/", 
                                    transition[which(c("-1->0", "-1->+1", "+1->-1", "+1->0") == from_to)], "_", i))
        
      }
      
    }
    
  }

}