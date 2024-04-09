#封装函数 代数模拟生成种群参考等位基因数目 根据二项分布

calculate_next_gen_num <- function(current_num, N, a1, b1, m, p3) {
  # 计算迁移和突变后的等位基因频率
  p0 <- current_num / (2 * N)
  p_after_migration <- p0 * (1 - m) + m * p3
  p_after_mutation <- p_after_migration * (1 - a1 - b1) + b1
  
  # 使用二项分布模拟下一代的等位基因数目
  rbinom(1, 2 * N, p_after_mutation)
}
#生成下一代的等位基因频率

#封装函数 存储每一代三个种群的参考等位基因数目

F_simulation_modified <- function(N, K, a1, b1, m_ij, p0, p3, num_generations) {
  # 初始化等位基因数目
  n_i <- matrix(nrow = K, ncol = num_generations)
  n_i[, 1] <- rep(2 * N * p0, K)
  n_i[3, ] <- rep(2 * N * p3, num_generations) # 第三个种群等位基因数目保持不变
  
  for (t in 1:(num_generations - 1)) {
    for (i in 1:2) {
      n_i[i, t + 1] <- calculate_next_gen_num(n_i[i, t], N, a1, b1, m_ij[i, 3], p3)
    }
  }
  
  # 计算每代的等位基因频率
  p_i <- n_i / (2 * N)
  return(p_i)
}


# 初始化参数 在该参数组合下从理论上来说F2是会出现拐点的
N <- 100  # 种群1和2的尺寸
a1 <- 0.0005  # A突变为a的概率
b1 <- 0.0002  # a突变为A的概率
m_ij <- matrix(0.015, nrow = 3, ncol = 3)
m_ij[1,2] <- m_ij[2,1] <- 0
diag(m_ij) <- 1 - 0.015 
m_ij[3,3] <- 1 - 0.015*2
p0 <- 0.2   # 前两个种群的初始等位基因频率
p3 <- 0.1121   # 第三个种群的固定等位基因频率
num_generations <- 2000  # 模拟的代数

# 修改后的模拟函数 更贴合理论部分的推导
repeat_simulation_2 <- function(N, num_simulations, K, a1, b1, m_ij, p0, p3, num_generations) {
  all_freqs <- matrix(0, nrow = num_generations, ncol = num_simulations)
  F2t <- matrix(0, nrow = num_generations, ncol = num_simulations) 
  for (sim in 1:num_simulations) {
    F_matrix <- F_simulation_modified(N, K, a1, b1, m_ij, p0, p3, num_generations)
    all_freqs[, sim] <- F_matrix[1, ] / (2 * N)  # 计算并存储种群1的等位基因频率
    F2t[, sim] <- (all_freqs[, sim] - p0) ^2
  }
  return(F2t)
}

set.seed(2)
num_simulations <- 1000  # 可以根据需要调整模拟次数
num_generations <- 500  # 模拟的代数
F2t_matrix <- repeat_simulation_2(N, num_simulations, 3, a1, b1, m_ij, p0, p3, num_generations)

# 计算F2值
F2_values <- rowMeans(F2t_matrix)

# 绘图
plot(1:num_generations, F2_values, type = 'l', xlab = "Generation", ylab = "Mean F2", lwd = 2,
     main = "Mean F2 over Generations")

# LOESS拟合“Mean F2 Over Generations”
generation_time <- 1:num_generations  # 创建一个时间（代数）序列
loess_fit_F2 <- loess(F2_values ~ generation_time, span = 0.28)  # 使用LOESS拟合

# 生成用于绘制拟合曲线的预测值
pred_time_vals <- seq(min(generation_time), max(generation_time), length.out = 1000)
pred_F2_vals <- predict(loess_fit_F2, newdata = data.frame(generation_time = pred_time_vals))

# 绘制Mean F2随代数变化的图，并添加LOESS拟合曲线
plot(generation_time, F2_values, type = 'l', xlab = "Generation" , ylab = expression(italic("Mean F")[bold(2)]), lwd =2,
     main = expression(italic("Mean F")[bold(2)] ~ "over Generations"), col = "slateblue")
lines(pred_time_vals, pred_F2_vals, col = "brown1", lwd = 2)  # 添加拟合曲线

# 最高点和最低点
time_interval <- pred_time_vals >= 200 & pred_time_vals <= 280
preds_interval <- pred_F2_vals[time_interval]
t_vals_interval <- pred_time_vals[time_interval]

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
       legend = c("(208, 0.03975635)", "(255, 0.03974821)"), # 图例文本
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
