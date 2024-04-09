
calculate_next_gen_num <- function(current_num, N, a1, b1, m, p3) {

  p0 <- current_num / (2 * N)
  p_after_migration <- p0 * (1 - m) + m * p3
  p_after_mutation <- p_after_migration * (1 - a1 - b1) + b1
  
  rbinom(1, 2 * N, p_after_mutation)
}


calculate_next_gen_num2 <- function(current_num, N, a1, b1, m1, p1, m2, p2) {
  p0 <- current_num / (2 * N)
  p_after_migration <- p0 * (1 - m1 - m2) + m1 * p1 + m2 * p2
  p_after_mutation <- p_after_migration * (1 - a1 - b1) + b1
  rbinom(1, 2 * N, p_after_mutation)
}

#----------All three populations are affected by linear evolutionary pressure, but what sets Population 3 apart is that its allele frequencies also change.

F_simulation_modified <- function(N, K, a1, b1, m_ij, p0, p3, num_generations) {
  n_i <- matrix(nrow = K, ncol = num_generations)
  n_i[1, 1] <- 2 * N * p0 
  n_i[2, 1] <- 2 * N * p0
  n_i[3, 1] <- 2 * N * p3  

  for (t in 1:(num_generations - 1)) {
    p1 = n_i[1 , t] / (2 * N)
    p2 = n_i[2 , t] / (2 * N)
    p3 = n_i[3, t] / (2 * N)
    n_i[1, t + 1] <- calculate_next_gen_num(n_i[1, t], N, a1, b1, m_ij[1, 3], p3)
    n_i[2, t + 1] <- calculate_next_gen_num(n_i[2, t], N, a1, b1, m_ij[2, 3], p3)
    n_i[3, t + 1] <- calculate_next_gen_num2(n_i[3, t], N, a1, b1, m_ij[3, 1], p1, m_ij[3, 2], p2)
  }
  p_i <- n_i / (2 * N)
  return(p_i)
}


repeat_simulation_2_modified <- function(N, num_simulations, K, a1, b1, m_ij, p0, p3, num_generations) {
  all_diffs_up <- matrix(0, nrow = num_generations, ncol = num_simulations)
  all_diffs_down <- matrix(0, nrow = num_generations, ncol = num_simulations)
  p1_diff <- matrix(0, nrow = num_generations, ncol = num_simulations)
  for (sim in 1:num_simulations) {
    F_matrix <- F_simulation_modified(N, K, a1, b1, m_ij, p0, p3, num_generations)
    
    p1_freqs <- F_matrix[1, ] / (2 * N)
    p2_freqs <- F_matrix[2, ] / (2 * N)

    p1_diff[, sim] <- (p1_freqs - p0)^2
    all_diffs_up[, sim] <- (p1_freqs - p2_freqs)^2
    all_diffs_down[, sim] <- (p1_freqs + p2_freqs - 2*p1_freqs*p2_freqs)
  }

  mean_diffs <- rowMeans(all_diffs_up) / rowMeans(all_diffs_down)
  mean_p1_diff <- rowMeans(p1_diff)
  return(list(mean_diffs = mean_diffs, mean_p1_diff = mean_p1_diff))
}



N <- 500  
a1 <- 0.0005 
b1 <- 0.0003  
m_ij <- matrix(0.005, nrow = 3, ncol = 3)
m_ij[1,2] <- m_ij[2,1] <- 0
diag(m_ij) <- 1 - 0.005 
m_ij[3,3] <- 1 - 0.005*2
p0 <- 0.3   
p3 <- 0.2   

num_generations <- 1000
num_simulations <- 1000 

set.seed(9)
results <- repeat_simulation_2_modified(N, num_simulations, 3, a1, b1, m_ij, p0, p3, num_generations)
generation_time <- 1:1000

mean_diffs <- results$mean_diffs
loess_fit_st <- loess(mean_diffs ~ generation_time, span = 0.3)  

pred_time_vals_st <- seq(min(generation_time), max(generation_time), length.out = 1000)
pred_vals_st <- predict(loess_fit_st, newdata = data.frame(generation_time = pred_time_vals_st))

F2_values <- results$mean_p1_diff
loess_fit_F2 <- loess(F2_values ~ generation_time, span = 0.3) 

pred_time_vals <- seq(min(generation_time), max(generation_time), length.out = 1000)
pred_F2_vals <- predict(loess_fit_F2, newdata = data.frame(generation_time = pred_time_vals))

par(mfrow=c(1, 2))

plot(generation_time, F2_values, type = 'l', xlab = expression(bold("time")), ylab = expression(italic("F")[bold(2)]), lwd =2,
     main = expression(italic("F")[bold(2)] ~ "over Time"), col = "slateblue")
lines(pred_time_vals, pred_F2_vals, col = "brown1", lwd = 2)  
mtext("A", side = 3, line = 1, at = par("usr")[1], adj = 0, cex = 1.5, font = 2)

# time_interval <- pred_time_vals >= 200 & pred_time_vals <= 280
# preds_interval <- pred_F2_vals[time_interval]
# t_vals_interval <- pred_time_vals[time_interval]

# max_val <- max(preds_interval)
# min_val <- min(preds_interval)
# max_time <- t_vals_interval[which.max(preds_interval)]
# min_time <- t_vals_interval[which.min(preds_interval)]
# max_x <- max_time
# max_y <- max_val
# points(max_x, max_y, col = "tan1", pch = 20, cex = 1.5)
# segments(max_x, par("usr")[3], max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
# segments(par("usr")[1], max_y, max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
# min_x <- min_time
# min_y <- min_val
# points(min_x, min_y, col = "limegreen", pch = 20, cex = 1.5)
# segments(min_x, par("usr")[3], min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)
# segments(par("usr")[1], min_y, min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)
# legend("bottomright",                      
#        legend = c("(208, 0.03975635)", "(255, 0.03974821)"), 
#        col = c("tan1", "limegreen"),               
#        pch = 20,                          
#        pt.cex = 1.5,                           
#        cex = 0.8,                          
#        bty = "n",                        
#        bg = 'white',                       
#        box.lwd = 1,                           
#        box.col = "grey50",                     
#        text.col = "black",                   
#        text.font = 2,                          
#        horiz = FALSE)                         
# print(paste("Max value point: (", max_time, ", ", max_val, ")", sep = ""))
# print(paste("Min value point: (", min_time, ", ", min_val, ")", sep = ""))

plot(generation_time, mean_diffs, type = 'l', xlab = expression(bold("time")), ylab = expression(italic("F")[bold(st)]), lwd =2,
     main = expression(italic("F")[bold(st)] ~ "over Time"), col = "slateblue")
lines(pred_time_vals_st, pred_vals_st, col = "brown1", lwd = 2)  
mtext("B", side = 3, line = 1, at = par("usr")[1], adj = 0, cex = 1.5, font = 2)

# time_interval2 <- pred_time_vals >= 200 & pred_time_vals <= 280
# preds_interval2 <- pred_vals2[time_interval2]
# t_vals_interval2 <- pred_time_vals2[time_interval2]
# max_val2 <- max(preds_interval2)
# min_val2 <- min(preds_interval2)
# max_time2 <- t_vals_interval2[which.max(preds_interval2)]
# min_time2 <- t_vals_interval2[which.min(preds_interval2)]
# max_x2 <- max_time2
# max_y2 <- max_val2
# points(max_x2, max_y2, col = "tan1", pch = 20, cex = 1.5)
# max_x2
# max_y2
# segments(max_x2, par("usr")[3], max_x2, max_y2, col = "tan1", lty = "dotted" , lwd = 2)
# segments(par("usr")[1], max_y2, max_x2, max_y2, col = "tan1", lty = "dotted" , lwd = 2)
# time_interval2 <- pred_time_vals2 >= 200 & pred_time_vals2 <= 400
# preds_interval2 <- pred_vals2[time_interval2]
# t_vals_interval2 <- pred_time_vals2[time_interval2]
# min_val2 <- min(preds_interval2)
# min_time2 <- t_vals_interval2[which.min(preds_interval2)]
# min_x2 <- min_time2
# min_y2 <- min_val2
# min_x2
# min_y2
# points(min_x2, min_y2, col = "limegreen", pch = 20, cex = 1.5)
# segments(min_x2, par("usr")[3], min_x2, min_y2, col = "limegreen", lty = "dotted" , lwd = 2)
# segments(par("usr")[1], min_y2, min_x2, min_y2, col = "limegreen", lty = "dotted" , lwd = 2)
# legend("bottomright",                         
#        legend = c("(233, 0.0000123)", "(310, 0.0000106)"),
#        col = c("tan1", "limegreen"),           
#        pch = 20,                         
#        pt.cex = 1.5,                       
#        cex = 0.8,                           
#        bty = "n",                        
#        bg = 'white',                        
#        box.lwd = 1,                   
#        box.col = "grey50",               
#        text.col = "black",                 
#        text.font = 2,                       
#        horiz = FALSE)                        
