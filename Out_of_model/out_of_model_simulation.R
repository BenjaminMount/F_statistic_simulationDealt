data <- read.table("path/F2_final.txt", header = F, sep = "\t")
data <- data[,-6]
colnames(data) <- c('time','P_A1','P_a1','P_A2','P_a2')
# 查看前几行数据
head(data)
data[,1] <- 1:4000
frequency <- data[,c(1,2)]
plot(frequency)
initial_frequency <- frequency[frequency$time == 1, 2]
# 计算F2t并添加到数据框
frequency$F2t <- (frequency[, 2] - initial_frequency) ^ 2
plot(frequency$F2t) #观察

# 应用Loess平滑拟合
loess_fit <- loess(F2t ~ time, data = frequency, span = 0.15)  # span 控制拟合的平滑程度
# 预测并绘制拟合曲线
t_vals <- seq(min(frequency$time), max(frequency$time), length.out = 1000)
preds <- predict(loess_fit, newdata = data.frame(time = t_vals))

# 绘制单个散点图
plot(frequency$time, frequency$F2t,
     col = "slateblue",   # 外环颜色
     pch = 21,            # 空心点类型
     xlab = "time", 
     ylab = expression(bold(F[2])),
     main = expression(italic(F[bold(2)]) ~ "over Time"),
     font.lab = 2
     )
# 添加内填充的点
points(frequency$time, frequency$F2t,
       col = "slateblue",   # 内填充颜色，这是一个浅蓝色
       pch = 21,
       bg = "slateblue",    # bg参数设置背景填充颜色
       cex = 0.8)         # cex参数使内部点比外部点小
lines(t_vals, preds, col = "brown1", lwd = 2)

# 找到时间区间 [500, 800] 内的最高点和最低点
time_interval <- t_vals >= 500 & t_vals <= 800
preds_interval <- preds[time_interval]
t_vals_interval <- t_vals[time_interval]

# 最高点和最低点
max_val <- max(preds_interval)
min_val <- min(preds_interval)
max_time <- t_vals_interval[which.max(preds_interval)]
min_time <- t_vals_interval[which.min(preds_interval)]
# 最大值点
max_x <- max_time
max_y <- max_val
points(max_x, max_y, col = "tan1", pch = 20, cex = 1.5)
# 从最大值点出发绘制垂直线到x轴
segments(max_x, par("usr")[3], max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
# 从最大值点出发绘制水平线到y轴
segments(par("usr")[1], max_y, max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
# 最小值点
min_x <- min_time
min_y <- min_val
points(min_x, min_y, col = "limegreen", pch = 20, cex = 1.5)
# 从最小值点出发绘制垂直线到x轴
segments(min_x, par("usr")[3], min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)
# 从最小值点出发绘制水平线到y轴
segments(par("usr")[1], min_y, min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)

# 添加图例
legend("bottomright",                          # 图例位置
       inset = c(0.05, 0.05),                  # 图例内缩，向图内移动
       legend = c("(501.4, 0.14)", "(745.6, 0.08)"), # 图例文本
       col = c("tan1", "limegreen"),               # 图例颜色
       pch = 20,                               # 图例点的类型
       pt.cex = 1.5,                           # 图例点的大小
       cex = 0.8,                              # 图例文本的大小
       bty = "n",                              # 不显示图例边框
       bg = 'white',                           # 图例背景颜色
       box.lwd = 1,                            # 图例边框宽度
       box.col = "grey50",                     # 图例边框颜色
       text.col = "black",                     # 文本颜色
       text.font = 2,                          # 文本字体
       horiz = FALSE)                          # 图例方向，FALSE为垂直，TRUE为水平

# 打印最高点和最低点坐标
print(paste("Max value point: (", max_time, ", ", max_val, ")", sep = ""))
print(paste("Min value point: (", min_time, ", ", min_val, ")", sep = ""))


#Fst
data_st <- read.table("path/Fst_final.txt", header = F, sep = "\t")
data_st <- data_st[,-8]
colnames(data_st) <- c('time','P_A1','P_a1','P_A2','P_a2','P_A3','P_a3')
data_st[,1] <- 1:6000
frequency_for_st <- data_st[,c(1,2,4)]
# 计算FSt并添加到数据框
frequency_for_st$Fst <- (frequency_for_st[, 2] - frequency_for_st[, 3]) / (1 - frequency_for_st[, 3])


# 应用Loess平滑拟合
loess_fit_st <- loess(Fst ~ time, data = frequency_for_st, span = 0.09)  # span 控制拟合的平滑程度
# 预测并绘制拟合曲线
t_vals_st <- seq(min(frequency_for_st$time), max(frequency_for_st$time), length.out = 1000)
preds_st <- predict(loess_fit_st, newdata = data.frame(time = t_vals_st))

plot(frequency_for_st$time, frequency_for_st$Fst,
     col = "slateblue",   # 外环颜色
     pch = 21,            # 空心点类型
     xlab = "time", 
     ylab = expression(bold(F[st])),
     main = expression(italic(F[italic(bold(st))]) ~ "over Time"),
     font.lab = 2)
# 添加内填充的点
points(frequency_for_st$time, frequency_for_st$Fst,
       col = "slateblue",   # 内填充颜色，这是一个浅蓝色
       pch = 21,
       bg = "slateblue",    # bg参数设置背景填充颜色
       cex = 0.8)         # cex参数使内部点比外部点小

lines(t_vals_st, preds_st, col = "brown1", lwd = 2)

# 找到时间区间 [500, 1100] 内的最高点和最低点
time_interval_st <- t_vals_st >= 300 & t_vals_st <= 800
preds_interval_st <- preds_st[time_interval_st]
t_vals_interval_st <- t_vals_st[time_interval_st]

# 最高点和最低点
max_val_st <- max(preds_interval_st)
min_val_st <- min(preds_interval_st)
max_time_st <- t_vals_interval_st[which.max(preds_interval_st)]
min_time_st <- t_vals_interval_st[which.min(preds_interval_st)]
# 最大值点
max_x_st <- max_time_st
max_y_st <- max_val_st
points(max_x_st, max_y_st, col = "tan1", pch = 20, cex = 1.5)
# 从最大值点出发绘制垂直线到x轴
segments(max_x_st, par("usr")[3], max_x_st, max_y_st, col = "tan1", lty = "dotted" , lwd = 2)
# 从最大值点出发绘制水平线到y轴
segments(par("usr")[1], max_y_st, max_x_st, max_y_st, col = "tan1", lty = "dotted" , lwd = 2)
# 最小值点
min_x_st <- min_time_st
min_y_st <- min_val_st
points(min_x_st, min_y_st, col = "limegreen", pch = 20, cex = 1.5)
# 从最小值点出发绘制垂直线到x轴
segments(min_x_st, par("usr")[3], min_x_st, min_y_st, col = "limegreen", lty = "dotted" , lwd = 2)
# 从最小值点出发绘制水平线到y轴
segments(par("usr")[1], min_y_st, min_x_st, min_y_st, col = "limegreen", lty = "dotted" , lwd = 2)

# 添加图例
legend("bottomright",                          # 图例位置
       inset = c(0.05, 0.05),                  # 图例内缩，向图内移动
       legend = c("(373.3, 0.12)", "(655.5, -1.32)"), # 图例文本
       col = c("tan1", "limegreen"),               # 图例颜色
       pch = 20,                               # 图例点的类型
       pt.cex = 1.5,                           # 图例点的大小
       cex = 0.8,                              # 图例文本的大小
       bty = "n",                              # 不显示图例边框
       bg = 'white',                           # 图例背景颜色
       box.lwd = 1,                            # 图例边框宽度
       box.col = "grey50",                     # 图例边框颜色
       text.col = "black",                     # 文本颜色
       text.font = 2,                          # 文本字体
       horiz = FALSE)                          # 图例方向，FALSE为垂直，TRUE为水平

# 打印最高点和最低点坐标
print(paste("Max value point: (", max_time_st, ", ", max_val_st, ")", sep = ""))
print(paste("Min value point: (", min_time_st, ", ", min_val_st, ")", sep = ""))



#绘制合并图
par(mfrow = c(1, 2))

plot(frequency$time, frequency$F2t,
     col = "slateblue",   # 外环颜色
     pch = 21,            # 空心点类型
     xlab = "time", 
     ylab = expression(bold(F[2])),
     main = expression(italic(F[bold(2)]) ~ "over Time"),
     font.lab = 2
)
# 添加内填充的点
points(frequency$time, frequency$F2t,
       col = "slateblue",   # 内填充颜色，这是一个浅蓝色
       pch = 21,
       bg = "slateblue",    # bg参数设置背景填充颜色
       cex = 0.8)         # cex参数使内部点比外部点小
lines(t_vals, preds, col = "brown1", lwd = 2)

# 找到时间区间 [500, 800] 内的最高点和最低点
time_interval <- t_vals >= 500 & t_vals <= 800
preds_interval <- preds[time_interval]
t_vals_interval <- t_vals[time_interval]

# 最高点和最低点
max_val <- max(preds_interval)
min_val <- min(preds_interval)
max_time <- t_vals_interval[which.max(preds_interval)]
min_time <- t_vals_interval[which.min(preds_interval)]
# 最大值点
max_x <- max_time
max_y <- max_val
points(max_x, max_y, col = "tan1", pch = 20, cex = 1.5)
# 从最大值点出发绘制垂直线到x轴
segments(max_x, par("usr")[3], max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
# 从最大值点出发绘制水平线到y轴
segments(par("usr")[1], max_y, max_x, max_y, col = "tan1", lty = "dotted" , lwd = 2)
# 最小值点
min_x <- min_time
min_y <- min_val
points(min_x, min_y, col = "limegreen", pch = 20, cex = 1.5)
# 从最小值点出发绘制垂直线到x轴
segments(min_x, par("usr")[3], min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)
# 从最小值点出发绘制水平线到y轴
segments(par("usr")[1], min_y, min_x, min_y, col = "limegreen", lty = "dotted" , lwd = 2)

# 添加图例
legend("bottomright",                          # 图例位置
       inset = c(0.05, 0.05),                  # 图例内缩，向图内移动
       legend = c("(501.4, 0.14)", "(745.6, 0.08)"), # 图例文本
       col = c("tan1", "limegreen"),               # 图例颜色
       pch = 20,                               # 图例点的类型
       pt.cex = 1.5,                           # 图例点的大小
       cex = 0.8,                              # 图例文本的大小
       bty = "n",                              # 不显示图例边框
       bg = 'white',                           # 图例背景颜色
       box.lwd = 1,                            # 图例边框宽度
       box.col = "grey50",                     # 图例边框颜色
       text.col = "black",                     # 文本颜色
       text.font = 2,                          # 文本字体
       horiz = FALSE)                          # 图例方向，FALSE为垂直，TRUE为水平

# 在第一个图形外部的左上角标记'A'
mtext("A", side = 3, line = 1, at = par("usr")[1], adj = 0, cex = 1.5, font = 2)




plot(frequency_for_st$time, frequency_for_st$Fst,
     col = "slateblue",   # 外环颜色
     pch = 21,            # 空心点类型
     xlab = "time", 
     ylab = expression(bold(F[st])),
     main = expression(italic(F[italic(bold(st))]) ~ "over Time"),
     font.lab = 2)
# 添加内填充的点
points(frequency_for_st$time, frequency_for_st$Fst,
       col = "slateblue",   # 内填充颜色，这是一个浅蓝色
       pch = 21,
       bg = "slateblue",    # bg参数设置背景填充颜色
       cex = 0.8)         # cex参数使内部点比外部点小

lines(t_vals_st, preds_st, col = "brown1", lwd = 2)

# 找到时间区间 [500, 1100] 内的最高点和最低点
time_interval_st <- t_vals_st >= 300 & t_vals_st <= 800
preds_interval_st <- preds_st[time_interval_st]
t_vals_interval_st <- t_vals_st[time_interval_st]

# 最高点和最低点
max_val_st <- max(preds_interval_st)
min_val_st <- min(preds_interval_st)
max_time_st <- t_vals_interval_st[which.max(preds_interval_st)]
min_time_st <- t_vals_interval_st[which.min(preds_interval_st)]
# 最大值点
max_x_st <- max_time_st
max_y_st <- max_val_st
points(max_x_st, max_y_st, col = "tan1", pch = 20, cex = 1.5)
# 从最大值点出发绘制垂直线到x轴
segments(max_x_st, par("usr")[3], max_x_st, max_y_st, col = "tan1", lty = "dotted" , lwd = 2)
# 从最大值点出发绘制水平线到y轴
segments(par("usr")[1], max_y_st, max_x_st, max_y_st, col = "tan1", lty = "dotted" , lwd = 2)
# 最小值点
min_x_st <- min_time_st
min_y_st <- min_val_st
points(min_x_st, min_y_st, col = "limegreen", pch = 20, cex = 1.5)
# 从最小值点出发绘制垂直线到x轴
segments(min_x_st, par("usr")[3], min_x_st, min_y_st, col = "limegreen", lty = "dotted" , lwd = 2)
# 从最小值点出发绘制水平线到y轴
segments(par("usr")[1], min_y_st, min_x_st, min_y_st, col = "limegreen", lty = "dotted" , lwd = 2)

# 添加图例
legend("bottomright",                          # 图例位置
       inset = c(0.05, 0.05),                  # 图例内缩，向图内移动
       legend = c("(373.3, 0.12)", "(655.5, -1.32)"), # 图例文本
       col = c("tan1", "limegreen"),               # 图例颜色
       pch = 20,                               # 图例点的类型
       pt.cex = 1.5,                           # 图例点的大小
       cex = 0.8,                              # 图例文本的大小
       bty = "n",                              # 不显示图例边框
       bg = 'white',                           # 图例背景颜色
       box.lwd = 1,                            # 图例边框宽度
       box.col = "grey50",                     # 图例边框颜色
       text.col = "black",                     # 文本颜色
       text.font = 2,                          # 文本字体
       horiz = FALSE)                          # 图例方向，FALSE为垂直，TRUE为水平

# 在第二个图形外部的左上角标记'B'
mtext("B", side = 3, line = 1, at = par("usr")[1], adj = 0, cex = 1.5, font = 2)
