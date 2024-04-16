#As same as In model part
calculate_next_gen_num <- function(current_num, N, a1, b1, m, p3) {
  p0 <- current_num / (2 * N)
  p_after_migration <- p0 * (1 - m) + m * p3
  p_after_mutation <- p_after_migration * (1 - a1 - b1) + b1
  
  rbinom(1, 2 * N, p_after_mutation)
}

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

#modify function for Fst scenario
repeat_simulation_2_modified <- function(N, num_simulations, K, a1, b1, m_ij, p0, p3, num_generations) {
  all_diffs_up <- matrix(0, nrow = num_generations, ncol = num_simulations)
  all_diffs_down <- matrix(0, nrow = num_generations, ncol = num_simulations)
  all_p1_freqs <- matrix(0, nrow = num_generations, ncol = num_simulations)  
  all_p2_freqs <- matrix(0, nrow = num_generations, ncol = num_simulations)  
  for (sim in 1:num_simulations) {
    F_matrix <- F_simulation_modified(N, K, a1, b1, m_ij, p0, p3, num_generations)
    p1_freqs <- F_matrix[1, ] 
    p2_freqs <- F_matrix[2, ] 
    all_p1_freqs[, sim] <- p1_freqs
    all_p2_freqs[, sim] <- p2_freqs
    all_diffs_up[, sim] <- (p1_freqs - p2_freqs)^2
    all_diffs_down[, sim] <- (p1_freqs + p2_freqs - 2*p1_freqs*p2_freqs)
  }
  mean_diffs <- rowMeans(all_diffs_up) / rowMeans(all_diffs_down)
  return(list(mean_diffs = mean_diffs, all_p1_freqs = all_p1_freqs, all_p2_freqs = all_p2_freqs))
}

N <- 100
a1 <- 0.0003  
b1 <- 0.0002  
m_ij <- matrix(0.005, nrow = 3, ncol = 3)
m_ij[1,2] <- m_ij[2,1] <- 0
diag(m_ij) <- 1 - 0.005 
m_ij[3,3] <- 1 - 0.005*2
p0 <- 0.4   
p3 <- 0.5   

num_generations <- 600
num_simulations <- 2000  
set.seed(2)
results <- repeat_simulation_2_modified(N, num_simulations, 3, a1, b1, m_ij, p0, p3, num_generations)
mean_diffs <- results$mean_diffs
all_p1_freqs <- results$all_p1_freqs
all_p2_freqs <- results$all_p2_freqs
generation_time2 <- 1:500
mean_diffs2 <- mean_diffs[1:500]
loess_fit2 <- loess(mean_diffs2 ~ generation_time2, span = 0.28)  
pred_time_vals2 <- seq(min(generation_time2), max(generation_time2), length.out = 1000)
pred_vals2 <- predict(loess_fit2, newdata = data.frame(generation_time2 = pred_time_vals2))

plot(generation_time2, mean_diffs2, type = 'l', xlab = "Generation", ylab = "Mean F2", lwd =2,
     main = "Mean Fst over Generations", col = "slateblue")
lines(pred_time_vals2, pred_vals2, col = "brown1", lwd = 2) 

time_interval2 <- pred_time_vals >= 200 & pred_time_vals <= 280
preds_interval2 <- pred_vals2[time_interval2]
t_vals_interval2 <- pred_time_vals2[time_interval2]

max_val2 <- max(preds_interval2)
min_val2 <- min(preds_interval2)
max_time2 <- t_vals_interval2[which.max(preds_interval2)]
min_time2 <- t_vals_interval2[which.min(preds_interval2)]
max_x2 <- max_time2
max_y2 <- max_val2
points(max_x2, max_y2, col = "tan1", pch = 20, cex = 1.5)
max_x2
max_y2
segments(max_x2, par("usr")[3], max_x2, max_y2, col = "tan1", lty = "dotted" , lwd = 2)
segments(par("usr")[1], max_y2, max_x2, max_y2, col = "tan1", lty = "dotted" , lwd = 2)

time_interval2 <- pred_time_vals2 >= 200 & pred_time_vals2 <= 400
preds_interval2 <- pred_vals2[time_interval2]
t_vals_interval2 <- pred_time_vals2[time_interval2]

min_val2 <- min(preds_interval2)
min_time2 <- t_vals_interval2[which.min(preds_interval2)]

min_x2 <- min_time2
min_y2 <- min_val2
min_x2
min_y2

points(min_x2, min_y2, col = "limegreen", pch = 20, cex = 1.5)
segments(min_x2, par("usr")[3], min_x2, min_y2, col = "limegreen", lty = "dotted" , lwd = 2)
segments(par("usr")[1], min_y2, min_x2, min_y2, col = "limegreen", lty = "dotted" , lwd = 2)

legend("bottomright",               
       legend = c("(233, 0.0000123)", "(310, 0.0000106)"), 
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















