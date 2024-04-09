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

repeat_simulation_2_modified <- function(N, num_simulations, K, a1, b1, m_ij, p0, p3, num_generations) {
  # 初始化存储模拟结果的矩阵
  all_diffs_up <- matrix(0, nrow = num_generations, ncol = num_simulations)
  all_diffs_down <- matrix(0, nrow = num_generations, ncol = num_simulations)
  for (sim in 1:num_simulations) {
    F_matrix <- F_simulation_modified(N, K, a1, b1, m_ij, p0, p3, num_generations)
    
    # 提取种群1和种群2的等位基因频率
    p1_freqs <- F_matrix[1, ] / (2 * N)
    p2_freqs <- F_matrix[2, ] / (2 * N)
    
    # 计算给定公式的值
    all_diffs_up[, sim] <- (p1_freqs - p2_freqs)^2
    all_diffs_down[, sim] <- (p1_freqs + p2_freqs - 2*p1_freqs*p2_freqs)
  }
  # 计算每个时间点上的均值
  mean_diffs <- rowMeans(all_diffs_up) / rowMeans(all_diffs_down)
  return(mean_diffs)
}

#封装函数 存储每一代三个种群的参考等位基因数目


# 初始化参数 在该参数组合下从理论上来说F2是会出现拐点的
N <- 1000  # 种群1和2的尺寸
a1 <- 0.0003  # A突变为a的概率
b1 <- 0.0002  # a突变为A的概率
m_ij <- matrix(0.005, nrow = 3, ncol = 3)
m_ij[1,2] <- m_ij[2,1] <- 0
diag(m_ij) <- 1 - 0.005 
m_ij[3,3] <- 1 - 0.005*2
p0 <- 0.4   # 前两个种群的初始等位基因频率
p3 <- 0.5   # 第三个种群的固定等位基因频率

num_generations <- 600
# 重新进行模拟，提取结果
num_simulations <- 2000  # 可以根据需要调整模拟次数

set.seed(2)
mean_diffs <- repeat_simulation_2_modified(N, num_simulations, 3, a1, b1, m_ij, p0, p3, num_generations)

generation_time2 <- 1:500
mean_diffs2 <- mean_diffs[1:500]
loess_fit2 <- loess(mean_diffs2 ~ generation_time2, span = 0.28)  # 使用LOESS拟合
# 生成用于绘制拟合曲线的预测值
pred_time_vals2 <- seq(min(generation_time2), max(generation_time2), length.out = 1000)
pred_vals2 <- predict(loess_fit2, newdata = data.frame(generation_time2 = pred_time_vals2))

# 绘制Mean Fst随代数变化的图，并添加LOESS拟合曲线
plot(generation_time2, mean_diffs2, type = 'l', xlab = "Generation", ylab = "Mean F2", lwd =2,
     main = "Mean Fst over Generations", col = "slateblue")
lines(pred_time_vals2, pred_vals2, col = "brown1", lwd = 2)  # 添加拟合曲线

# 最高点和最低点
time_interval2 <- pred_time_vals >= 200 & pred_time_vals <= 280
preds_interval2 <- pred_vals2[time_interval2]
t_vals_interval2 <- pred_time_vals2[time_interval2]

# 最高点和最低点
max_val2 <- max(preds_interval2)
min_val2 <- min(preds_interval2)
max_time2 <- t_vals_interval2[which.max(preds_interval2)]
min_time2 <- t_vals_interval2[which.min(preds_interval2)]
# 最大值点
max_x2 <- max_time2
max_y2 <- max_val2
points(max_x2, max_y2, col = "tan1", pch = 20, cex = 1.5)
max_x2
max_y2
# 从最大值点出发绘制垂直线到x轴
segments(max_x2, par("usr")[3], max_x2, max_y2, col = "tan1", lty = "dotted" , lwd = 2)
# 从最大值点出发绘制水平线到y轴
segments(par("usr")[1], max_y2, max_x2, max_y2, col = "tan1", lty = "dotted" , lwd = 2)

# 最高点和最低点
time_interval2 <- pred_time_vals2 >= 200 & pred_time_vals2 <= 400
preds_interval2 <- pred_vals2[time_interval2]
t_vals_interval2 <- pred_time_vals2[time_interval2]

min_val2 <- min(preds_interval2)
min_time2 <- t_vals_interval2[which.min(preds_interval2)]

# 最小值点
min_x2 <- min_time2
min_y2 <- min_val2
min_x2
min_y2

points(min_x2, min_y2, col = "limegreen", pch = 20, cex = 1.5)
# 从最小值点出发绘制垂直线到x轴
segments(min_x2, par("usr")[3], min_x2, min_y2, col = "limegreen", lty = "dotted" , lwd = 2)
# 从最小值点出发绘制水平线到y轴
segments(par("usr")[1], min_y2, min_x2, min_y2, col = "limegreen", lty = "dotted" , lwd = 2)

# 添加图例
legend("bottomright",                          # 图例位置
       legend = c("(233, 0.0000123)", "(310, 0.0000106)"), # 图例文本
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















