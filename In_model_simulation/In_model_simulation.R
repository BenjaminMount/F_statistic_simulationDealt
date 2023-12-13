#First, if you only want to simulate a scenario with two populations, you can use the F2_simulation which is commented out part.

# F2_simulation <- function(N, a1, b1, m12, m21, p0, num_generations) {
#   p_i <- matrix(rep(p0, 2 * num_generations), nrow = 2)
#   
#   for (t in 1:(num_generations - 1)) {
#     Z_i <- rbinom(2, 2 * N, p_i[, t])
#     
#     p1_next <- (1 - m12) * ((1 - a1) * Z_i[1] / (2 * N) + b1 * (1 - Z_i[1] / (2 * N))) + 
#       m21 * ((1 - a1) * Z_i[2] / (2 * N) + b1 * (1 - Z_i[2] / (2 * N)))
#     p2_next <- m12 * ((1 - a1) * Z_i[1] / (2 * N) + b1 * (1 - Z_i[1] / (2 * N))) + 
#       (1 - m21) * ((1 - a1) * Z_i[2] / (2 * N) + b1 * (1 - Z_i[2] / (2 * N)))
#     
#     p_i[1, t + 1] <- p1_next
#     p_i[2, t + 1] <- p2_next
#   }
#   
#   return(p_i)
# }
# 
# N <- 50  
# a1 <- 0.01  
# b1 <- 0.02  
# m12 <- 0.05  
# m21 <- 0.05  
# p0 <- 0.60  
# num_generations <- 300  
# set.seed(1)
# F2_matrix <- F2_simulation(N, a1, b1, m12, m21, p0, num_generations)
# population_1_F2 <- as.data.frame(t(F2_matrix)[,1])
# population_1_F2$time <- 1:300 
# colnames(population_1_F2) <-c("P(A)","time")
# p0 <- population_1_F2[population_1_F2$time == 2, 1]
# population_1_F2$F2 <- (population_1_F2[, 1] - p0) ^ 2
# 
# loess_fit <- loess(F2 ~ time, data = population_1_F2, span = 0.1)
# t_vals <- seq(min(population_1_F2$time), max(population_1_F2$time), length.out = 1000)
# preds <- predict(loess_fit, newdata = data.frame(time = t_vals))
# 
# plot(population_1_F2$time, population_1_F2$F2,
#      col = "slateblue",
#      pch = 21,           
#      xlab = "time", 
#      ylab = "F-2", 
#      main = "F-2-statistic over time",
#      font.lab = 2
# )
# points(population_1_F2$time, population_1_F2$F2,
#        col = "slateblue",   
#        pch = 21,
#        bg = "slateblue",    
#        cex = 0.8)       
# lines(t_vals, preds, col = "brown1", lwd = 2)

#After expanding the function to the scenario about K populations, we can obtain a function as followed.
#Design the function
F_simulation <- function(N, K, a1, b1, m_ij, p0, num_generations) {
  if(dim(m_ij)[1] != K || dim(m_ij)[2] != K) {
    stop("The migration matrix m_ij must be a K x K matrix.")
  }
  p_i <- matrix(rep(p0, K * num_generations), nrow = K)
  
  for (t in 1:(num_generations - 1)) {
    Z_i <- rbinom(K, 2 * N, p_i[, t])
    
    for (i in 1:K) {
      migrant_freq <- 0
      for (j in 1:K) {
        migrant_freq <- migrant_freq + m_ij[i, j] * ((1 - a1) * Z_i[j] / (2 * N) + b1 * (1 - Z_i[j] / (2 * N)))
        #Both of migration and mutation are considered  
      }
      p_i[i, t + 1] <- migrant_freq
      #loop, Iterate
    }
  }
  return(p_i)
}


N <- 50  #each population size
a1 <- 0.01  #mutation rate (A to a)
b1 <- 0.02  #mutation rate (a to A)

m_ij <- matrix(0, nrow = 4, ncol = 4) #migration matrix
m_ij[1, 3] <- 0.03  
m_ij[1, 4] <- 0.02  
m_ij[2, 3] <- 0.02  
m_ij[2, 4] <- 0.03  
m_ij[3, 1] <- m_ij[3, 2] <- m_ij[3, 4] <- m_ij[4, 1] <- m_ij[4, 2] <- m_ij[4, 3] <- 0.05
m_ij[3, 3] <- m_ij[4, 4] <- 0.85
m_ij[1, 1] <- m_ij[2, 2] <- 0.95
p0 <- 0.6  
num_generations <- 300  

set.seed(1)
F_matrix <- F_simulation(N, 4, a1, b1, m_ij, p0, num_generations)
#write.csv(F_matrix, file = "path/F_matrix.csv", row.names = FALSE)

population_1_2 <- as.data.frame(t(F_matrix)[,1:2])
population_1_2$time <- 1:num_generations 
colnames(population_1_2) <-c("P(A1)","P(A2)","time")
#calculate the F-statistic
population_1_2$Population1_F2 <- (population_1_2[, 1] - p0) ^ 2
population_1_2$Population2_F2 <- (population_1_2[, 2] - p0) ^ 2
population_1_2$Fst <- (population_1_2[, 1] - population_1_2[, 2]) / (1 - population_1_2[, 2])
#LOESS fitting
loess_fit <- loess(Population1_F2 ~ time, data = population_1_2, span = 0.1)  
t_vals <- seq(min(population_1_2$time), max(population_1_2$time), length.out = 1000)
preds <- predict(loess_fit, newdata = data.frame(time = t_vals))

#Scatterplot
plot(population_1_2$time, population_1_2$Population1_F2,
     col = "slateblue",   
     pch = 21,            
     xlab = "time",
     ylab = expression(bold(F[2])),
     main = expression(italic(F[bold(2)]) ~ "over Time"),
     font.lab = 2
)
points(population_1_2$time, population_1_2$Population1_F2,
       col = "slateblue",  
       pch = 21,
       bg = "slateblue",    
       cex = 0.8)         
lines(t_vals, preds, col = "brown1", lwd = 2)

# loess_fit_2 <- loess(Population2_F2 ~ time, data = population_1_2, span = 0.1) 
# t_vals_2 <- seq(min(population_1_2$time), max(population_1_2$time), length.out = 1000)
# preds_2 <- predict(loess_fit_2, newdata = data.frame(time = t_vals))
# 
# plot(population_1_2$time, population_1_2$Population2_F2,
#      col = "slateblue",   
#      pch = 21,            
#      xlab = "time",
#      ylab = "F-2",
#      main = "F-2-statistic over time",
#      font.lab = 2
# )
# points(population_1_2$time, population_1_2$Population2_F2,
#        col = "slateblue",   
#        pch = 21,
#        bg = "slateblue",   
#        cex = 0.8)        
# lines(t_vals_2, preds_2, col = "brown1", lwd = 2)

loess_fit_st <- loess(Fst ~ time, data = population_1_2, span = 0.15) 
t_vals_st <- seq(min(population_1_2$time), max(population_1_2$time), length.out = 1000)
preds_st <- predict(loess_fit_st, newdata = data.frame(time = t_vals_st))

plot(population_1_2$time, population_1_2$Fst,
     col = "slateblue",   
     pch = 21,          
     xlab = "time", 
     ylab = expression(bold(F[st])),
     main = expression(italic(F[italic(bold(st))]) ~ "over Time"),
     font.lab = 2)
points(population_1_2$time, population_1_2$Fst,
       col = "slateblue",   
       pch = 21,
       bg = "slateblue",    
       cex = 0.8)        

lines(t_vals_st, preds_st, col = "brown1", lwd = 2)

#combine the two plots into the one
par(mfrow= c(1,2))
plot(population_1_2$time, population_1_2$Population1_F2,
     col = "slateblue",   
     pch = 21,            
     xlab = "time",
     ylab = expression(bold(F[2])),
     main = expression(italic(F[bold(2)]) ~ "over Time"),
     font.lab = 2
)
points(population_1_2$time, population_1_2$Population1_F2,
       col = "slateblue",   
       pch = 21,
       bg = "slateblue",   
       cex = 0.8)         
lines(t_vals, preds, col = "brown1", lwd = 2)
mtext("A", side = 3, line = 1, at = par("usr")[1], adj = 0, cex = 1.5, font = 2)

plot(population_1_2$time, population_1_2$Fst,
     col = "slateblue",   
     pch = 21,            
     xlab = "time", 
     ylab = expression(bold(F[st])),
     main = expression(italic(F[italic(bold(st))]) ~ "over Time"),
     font.lab = 2)
points(population_1_2$time, population_1_2$Fst,
       col = "slateblue",   
       pch = 21,
       bg = "slateblue",    
       cex = 0.8)         
lines(t_vals_st, preds_st, col = "brown1", lwd = 2)
mtext("B", side = 3, line = 1, at = par("usr")[1], adj = 0, cex = 1.5, font = 2)
