#Algebraic simulations generate the number of reference alleles in populations based on the binomial distribution.
calculate_next_gen_num <- function(current_num, N, a1, b1, m, p3) {
  p0 <- current_num / (2 * N)
  p_after_migration <- p0 * (1 - m) + m * p3
  p_after_mutation <- p_after_migration * (1 - a1 - b1) + b1
  rbinom(1, 2 * N, p_after_mutation)
}

#Store the number of reference alleles for each of the three populations in every generation.

F_simulation_modified <- function(N, K, a1, b1, m_ij, p0, p3, num_generations) {
  n_i <- matrix(nrow = K, ncol = num_generations)
  n_i[, 1] <- rep(2 * N * p0, K)
  n_i[3, ] <- rep(2 * N * p3, num_generations) 
  
  for (t in 1:(num_generations - 1)) {
    for (i in 1:2) {
      n_i[i, t + 1] <- calculate_next_gen_num(n_i[i, t], N, a1, b1, m_ij[i, 3], p3)
    }
  }
  
  p_i <- n_i / (2 * N)
  return(p_i)
}

# Under this combination of parameters, theoretically, an inflection point will occur in F2.
N <- 100  
a1 <- 0.0005  
b1 <- 0.0002  
m_ij <- matrix(0.015, nrow = 3, ncol = 3)
m_ij[1,2] <- m_ij[2,1] <- 0
diag(m_ij) <- 1 - 0.015 
m_ij[3,3] <- 1 - 0.015*2
p0 <- 0.2   
p3 <- 0.1121   
num_generations <- 2000  

repeat_simulation_2 <- function(N, num_simulations, K, a1, b1, m_ij, p0, p3, num_generations) {
  all_freqs <- matrix(0, nrow = num_generations, ncol = num_simulations)
  F2t <- matrix(0, nrow = num_generations, ncol = num_simulations) 
  for (sim in 1:num_simulations) {
    F_matrix <- F_simulation_modified(N, K, a1, b1, m_ij, p0, p3, num_generations)
    all_freqs[, sim] <- F_matrix[1, ] 
    F2t[, sim] <- (all_freqs[, sim] - p0) ^2
  }
  return(list(all_freqs = all_freqs, F2t = F2t))
}

set.seed(123)
num_simulations <- 1000  
num_generations <- 500  

results <- repeat_simulation_2(N, num_simulations, 3, a1, b1, m_ij, p0, p3, num_generations)
all_freqs_matrix <- results$all_freqs
F2t_matrix <- results$F2t

F2_values <- rowMeans(F2t_matrix)
plot(1:num_generations, F2_values, type = 'l', xlab = "Generation", ylab = "Mean F2", lwd = 2,
     main = "Mean F2 over Generations")

# LOESS “Mean F2 Over Generations”
generation_time <- 1:num_generations  
loess_fit_F2 <- loess(F2_values ~ generation_time, span = 0.28)  

pred_time_vals <- seq(min(generation_time), max(generation_time), length.out = 1000)
pred_F2_vals <- predict(loess_fit_F2, newdata = data.frame(generation_time = pred_time_vals))

plot(generation_time, F2_values, type = 'l', xlab = "Generation" , ylab = expression(italic("Mean F")[bold(2)]), lwd =2,
     main = expression(italic("Mean F")[bold(2)] ~ "over Generations"), col = "slateblue")
lines(pred_time_vals, pred_F2_vals, col = "brown1", lwd = 2) 

time_interval <- pred_time_vals >= 200 & pred_time_vals <= 280
preds_interval <- pred_F2_vals[time_interval]
t_vals_interval <- pred_time_vals[time_interval]

max_val <- max(preds_interval)
min_val <- min(preds_interval)
max_time <- t_vals_interval[which.max(preds_interval)]
min_time <- t_vals_interval[which.min(preds_interval)]

max_x <- max_time
max_y <- max_val
points(max_x, max_y, col = "tan1", pch = 20, cex = 1.5)
segments(max_x, par("usr")[3], max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
segments(par("usr")[1], max_y, max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
min_x <- min_time
min_y <- min_val
points(min_x, min_y, col = "limegreen", pch = 20, cex = 1.5)
segments(min_x, par("usr")[3], min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)
segments(par("usr")[1], min_y, min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)

legend("bottomright",                          
       legend = c("(208, 0.03975635)", "(255, 0.03974821)"), 
       col = c("tan1", "limegreen"),               
       pch = 20,                            
       pt.cex = 1.5,                          
       cex = 0.8,                         
       bty = "n",                            
       bg = 'white',                   
       box.lwd = 1,                     
       box.col = "grey50",              
       text.col = "black",                
       text.font = 2,                  
       horiz = FALSE)                        

# 打印最高点和最低点坐标
print(paste("Max value point: (", max_time, ", ", max_val, ")", sep = ""))
print(paste("Min value point: (", min_time, ", ", min_val, ")", sep = ""))
