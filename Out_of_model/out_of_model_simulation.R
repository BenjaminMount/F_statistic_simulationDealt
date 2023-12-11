data <- read.table("path/F2_final.txt", header = F, sep = "\t")
data <- data[,-6]
colnames(data) <- c('time','P_A1','P_a1','P_A2','P_a2')
head(data)
data[,1] <- 1:4000
frequency <- data[,c(1,2)]
plot(frequency) #preview
initial_frequency <- frequency[frequency$time == 1, 2]
#calculate F_2
frequency$F2t <- (frequency[, 2] - initial_frequency) ^ 2
plot(frequency$F2t) #preview

# LOESS Fitting
loess_fit <- loess(F2t ~ time, data = frequency, span = 0.15)  
t_vals <- seq(min(frequency$time), max(frequency$time), length.out = 1000)
preds <- predict(loess_fit, newdata = data.frame(time = t_vals))

# Scatter Plot
plot(frequency$time, frequency$F2t,
     col = "slateblue",   
     pch = 21,           
     xlab = "time", 
     ylab = expression(bold(F[2])),
     main = expression(italic(F[bold(2)]) ~ "over Time"),
     font.lab = 2
     )
points(frequency$time, frequency$F2t,
       col = "slateblue",   
       pch = 21,
       bg = "slateblue",    
       cex = 0.8)         
lines(t_vals, preds, col = "brown1", lwd = 2)

time_interval <- t_vals >= 500 & t_vals <= 800
preds_interval <- preds[time_interval]
t_vals_interval <- t_vals[time_interval]

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
       inset = c(0.05, 0.05),                  
       legend = c("(501.4, 0.14)", "(745.6, 0.08)"), 
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
print(paste("Max value point: (", max_time, ", ", max_val, ")", sep = ""))
print(paste("Min value point: (", min_time, ", ", min_val, ")", sep = ""))

#Fst
data_st <- read.table("path/Fst_final.txt", header = F, sep = "\t")
data_st <- data_st[,-8]
colnames(data_st) <- c('time','P_A1','P_a1','P_A2','P_a2','P_A3','P_a3')
data_st[,1] <- 1:6000
frequency_for_st <- data_st[,c(1,2,4)]
# calculate F_st
frequency_for_st$Fst <- (frequency_for_st[, 2] - frequency_for_st[, 3]) / (1 - frequency_for_st[, 3])

# LOESS Fitting
loess_fit_st <- loess(Fst ~ time, data = frequency_for_st, span = 0.09)  # span 控制拟合的平滑程度
t_vals_st <- seq(min(frequency_for_st$time), max(frequency_for_st$time), length.out = 1000)
preds_st <- predict(loess_fit_st, newdata = data.frame(time = t_vals_st))

#Scatter Plot
plot(frequency_for_st$time, frequency_for_st$Fst,
     col = "slateblue",  
     pch = 21,            
     xlab = "time", 
     ylab = expression(bold(F[st])),
     main = expression(italic(F[italic(bold(st))]) ~ "over Time"),
     font.lab = 2)
points(frequency_for_st$time, frequency_for_st$Fst,
       col = "slateblue",   
       pch = 21,
       bg = "slateblue",  
       cex = 0.8)       
lines(t_vals_st, preds_st, col = "brown1", lwd = 2)

time_interval_st <- t_vals_st >= 300 & t_vals_st <= 800
preds_interval_st <- preds_st[time_interval_st]
t_vals_interval_st <- t_vals_st[time_interval_st]
max_val_st <- max(preds_interval_st)
min_val_st <- min(preds_interval_st)
max_time_st <- t_vals_interval_st[which.max(preds_interval_st)]
min_time_st <- t_vals_interval_st[which.min(preds_interval_st)]
max_x_st <- max_time_st
max_y_st <- max_val_st
points(max_x_st, max_y_st, col = "tan1", pch = 20, cex = 1.5)
segments(max_x_st, par("usr")[3], max_x_st, max_y_st, col = "tan1", lty = "dotted" , lwd = 2)
segments(par("usr")[1], max_y_st, max_x_st, max_y_st, col = "tan1", lty = "dotted" , lwd = 2)
min_x_st <- min_time_st
min_y_st <- min_val_st
points(min_x_st, min_y_st, col = "limegreen", pch = 20, cex = 1.5)
segments(min_x_st, par("usr")[3], min_x_st, min_y_st, col = "limegreen", lty = "dotted" , lwd = 2)
segments(par("usr")[1], min_y_st, min_x_st, min_y_st, col = "limegreen", lty = "dotted" , lwd = 2)

legend("bottomright",                         
       inset = c(0.05, 0.05),                 
       legend = c("(373.3, 0.12)", "(655.5, -1.32)"), 
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
print(paste("Max value point: (", max_time_st, ", ", max_val_st, ")", sep = ""))
print(paste("Min value point: (", min_time_st, ", ", min_val_st, ")", sep = ""))

#combine the two plot into the one
par(mfrow = c(1, 2))
plot(frequency$time, frequency$F2t,
     col = "slateblue",   
     pch = 21,           
     xlab = "time", 
     ylab = expression(bold(F[2])),
     main = expression(italic(F[bold(2)]) ~ "over Time"),
     font.lab = 2
)
points(frequency$time, frequency$F2t,
       col = "slateblue",   
       pch = 21,
       bg = "slateblue",   
       cex = 0.8)         
lines(t_vals, preds, col = "brown1", lwd = 2)
points(max_x, max_y, col = "tan1", pch = 20, cex = 1.5)
segments(max_x, par("usr")[3], max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
segments(par("usr")[1], max_y, max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
points(min_x, min_y, col = "limegreen", pch = 20, cex = 1.5)
segments(min_x, par("usr")[3], min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)
segments(par("usr")[1], min_y, min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)

legend("bottomright",                          
       inset = c(0.05, 0.05),                
       legend = c("(501.4, 0.14)", "(745.6, 0.08)"), 
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
mtext("A", side = 3, line = 1, at = par("usr")[1], adj = 0, cex = 1.5, font = 2)

plot(frequency_for_st$time, frequency_for_st$Fst,
     col = "slateblue",   
     pch = 21,            
     xlab = "time", 
     ylab = expression(bold(F[st])),
     main = expression(italic(F[italic(bold(st))]) ~ "over Time"),
     font.lab = 2)
points(frequency_for_st$time, frequency_for_st$Fst,
       col = "slateblue",   
       pch = 21,
       bg = "slateblue",    
       cex = 0.8)      
lines(t_vals_st, preds_st, col = "brown1", lwd = 2)

points(max_x_st, max_y_st, col = "tan1", pch = 20, cex = 1.5)
segments(max_x_st, par("usr")[3], max_x_st, max_y_st, col = "tan1", lty = "dotted" , lwd = 2)
segments(par("usr")[1], max_y_st, max_x_st, max_y_st, col = "tan1", lty = "dotted" , lwd = 2)
points(min_x_st, min_y_st, col = "limegreen", pch = 20, cex = 1.5)
segments(min_x_st, par("usr")[3], min_x_st, min_y_st, col = "limegreen", lty = "dotted" , lwd = 2)
segments(par("usr")[1], min_y_st, min_x_st, min_y_st, col = "limegreen", lty = "dotted" , lwd = 2)

legend("bottomright",                          
       inset = c(0.05, 0.05),                 
       legend = c("(373.3, 0.12)", "(655.5, -1.32)"), 
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
mtext("B", side = 3, line = 1, at = par("usr")[1], adj = 0, cex = 1.5, font = 2)
