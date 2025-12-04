# ==============================================================================
# 主函数: 智能自动过滤 (Smart Auto Filter) - 更新版
# ==============================================================================
#' 自动识别分布形态并计算过滤阈值 (修复版)
#' @description
#' 这个函数用于根据 nFeature 或 nCount 的分布（双峰或单峰长尾），
#' 自动选择合适的算法（谷底法或拐点法）来计算过滤阈值。
#' 单峰模式下，采用“二阶导数最右侧局部峰值法”，寻找最紧贴主峰的拐角。
#'
#' @param values 数值向量。例如 Seurat 对象的 nFeature_RNA。
#' @param adjust 数值。核密度估计的平滑度，建议 2.0。
#' @param valley_strictness 数值 (0-1)。判定双峰的严格程度。
#' @param slope_threshold 数值。单峰模式下的斜率阈值。
#' @param show_plot 逻辑值。是否绘图。
#' @param title 字符串。图片标题。
#' @param adjust 平滑度。对于导数计算，平滑非常重要，建议 2.0 或更高。
#' @param min_peak_ratio 数值。副峰高度至少要是主峰高度的多少倍才被考虑？
#'   默认 1e-4 (0.01%)，用于处理极其悬殊的双峰分布。
#' @param valley_strictness 数值。谷底判定阈值，默认 0.99 (只要谷底比副峰低1%就算双峰)。
#' @return 一个数值 (cutoff)。  并可在Rplot中展示图像(show_plot=TRUE)
#'
#' @import ggplot2
#' @import Seurat
#' @importFrom stats density
#' 
#' @export
smart_qc_cutoff <- function(values, adjust = 2.0, valley_strictness = 0.99, min_peak_ratio = 1e-4, 
                            assumption = "gaussian", show_plot = TRUE, title = "") {
  
  # 1. 数据变换
  raw_values <- values
  if (assumption == "poisson") {
    values <- sqrt(values)
  }
  
  # 2. 计算密度
  d <- density(values, adjust = adjust)
  
  # 3. 极值点检测
  diff_sign <- diff(sign(diff(d$y)))
  min_indices <- which(diff_sign == 2) + 1 
  max_indices <- which(diff_sign == -2) + 1
  
  # 4. 双峰检测 (传入新的 min_peak_ratio 参数)
  valley_res <- .get_valley_cutoff(d, min_indices, max_indices, min_peak_ratio = min_peak_ratio)
  
  method_used <- "Unknown"
  calc_cutoff <- 0
  is_bimodal <- FALSE
  
  if (!is.null(valley_res)) {
    # 检查谷底深度
    if (valley_res$ratio < valley_strictness) {
      is_bimodal <- TRUE
      method_used <- "Bimodal (Deep Valley)"
      calc_cutoff <- valley_res$cutoff
    } else {
      # 找到了坑，但是不够深 (ratio > 0.99)，视为平缓的鞍部
      method_used <- "Unimodal (Last Knee) - Valley too shallow"
    }
  } else {
    method_used <- "Unimodal (Last Knee) - No left peak detected"
  }
  
  # 5. 单峰检测 (二阶导数法)
  if (!is_bimodal) {
    calc_cutoff <- .get_curvature_cutoff(d)
    
    # 保护措施：如果计算出的 cutoff 非常小 (比如 < 10)，可能是误判了左侧噪音
    # 这种情况下，尝试回退到谷底结果 (即使它很浅)
    # 这是一个启发式策略：如果二阶导数切在 0 附近，还不如切在那个浅浅的坑里
    if (calc_cutoff < 10 && !is.null(valley_res)) {
      calc_cutoff <- valley_res$cutoff
      method_used <- "Bimodal (Shallow Valley Fallback)"
      is_bimodal <- TRUE
    }
    
    if (method_used == "Unknown") method_used <- "Unimodal (Last Knee)"
  }
  
  # 6. 还原阈值
  final_cutoff <- calc_cutoff
  if (assumption == "poisson") {
    final_cutoff <- calc_cutoff^2
    if (final_cutoff < 0) final_cutoff <- 0
  }
  if (final_cutoff < 10) final_cutoff <- max(10, final_cutoff)
  
  # 7. 绘图
  if (show_plot) {
    d_plot <- density(raw_values, adjust = adjust) 
    df <- data.frame(x = d_plot$x, y = d_plot$y)
    
    p <- ggplot(df, aes(x = x)) +
      geom_area(aes(y = y), fill = ifelse(is_bimodal, "#a6cee3", "#fdbf6f"), alpha = 0.5) +
      geom_line(aes(y = y), size = 1, color = "black") +
      
      geom_vline(xintercept = final_cutoff, color = "red", linetype = "dashed", size = 1) +
      
      annotate("text", x = final_cutoff, y = max(d_plot$y)*0.6, 
               label = paste0("Cutoff: ", round(final_cutoff, 0)), 
               color = "red", angle = 90, vjust = -0.5) +
      
      theme_classic() +
      labs(title = paste0(title, "\nDetected: ", method_used),
           subtitle = paste0("Valley Ratio: ", ifelse(!is.null(valley_res), round(valley_res$ratio, 3), "NA")),
           x = "Value", y = "Density")
    
    print(p)
  }
  
  message(paste0("[", title, "] 模式: ", method_used, " | 阈值: ", round(final_cutoff, 2)))
  return(final_cutoff)
}

# ==============================================================================
# 子函数 1: 寻找谷底 (增强版：支持悬殊双峰)
# ==============================================================================
.get_valley_cutoff <- function(d, min_indices, max_indices, min_peak_ratio = 1e-4) {
  
  # 1. 确定全局主峰 (Main Peak)
  # 规则：找最靠右的、高度显著的峰
  global_max_height <- max(d$y)
  
  # 这里的阈值(0.05)是为了确定"谁是老大"，防止选中极右侧的微小拖尾
  main_candidates <- max_indices[d$y[max_indices] > (global_max_height * 0.05)]
  if (length(main_candidates) == 0) return(NULL)
  main_peak_idx <- main_candidates[which.max(d$x[main_candidates])]
  
  # 2. 寻找左侧所有候选副峰 (Candidate Left Peaks)
  # 【关键修改】大幅降低门槛，允许副峰非常矮 (例如只有主峰的 0.01%)
  # 只要它是局部极大值，且不是底噪，我们就考虑它
  potential_left_peaks <- max_indices[max_indices < main_peak_idx]
  
  # 过滤掉几乎贴地的噪音 (高度 < 主峰的 0.01%)
  left_peaks <- potential_left_peaks[d$y[potential_left_peaks] > (global_max_height * min_peak_ratio)]
  
  if (length(left_peaks) == 0) return(NULL)
  
  # 3. 遍历所有左侧峰，寻找“最佳分割点”
  # 我们不假设离主峰最近的才是副峰，可能是左边更远的那个。
  # 评判标准：谷底相对于副峰的深度 (Ratio)。Ratio 越小，说明坑越深，分割越可靠。
  
  best_result <- NULL
  best_ratio <- 1.0 # 初始值，越小越好
  
  for (sec_idx in left_peaks) {
    
    # 在该副峰(sec_idx) 和 主峰(main_peak_idx) 之间找谷底
    valid_valleys <- min_indices[min_indices > sec_idx & min_indices < main_peak_idx]
    
    if (length(valid_valleys) > 0) {
      # 找到这一对峰之间的最低点
      valley_idx <- valid_valleys[which.min(d$y[valid_valleys])]
      
      valley_h <- d$y[valley_idx]
      sec_peak_h <- d$y[sec_idx]
      
      # 计算相对深度：谷底 / 副峰高度
      # 注意：这里只除以副峰高度。因为在悬殊双峰中，副峰通常远小于主峰。
      # 只要谷底比副峰低，就是一个有效的切割位。
      current_ratio <- valley_h / sec_peak_h
      
      # 如果这个坑比之前的更深 (ratio更小)，记录下来
      # 或者如果 ratio 差不多，选靠右的那个 (保留更多细胞)
      if (current_ratio < best_ratio) {
        best_ratio <- current_ratio
        best_result <- list(cutoff = d$x[valley_idx], ratio = current_ratio)
      }
    }
  }
  
  return(best_result)
}

# ==============================================================================
# 核心算法: 寻找二阶导数最靠近主峰的那个峰值
# ==============================================================================
.get_curvature_cutoff <- function(d) {
  
  # 1. 寻找主峰 (Robust Main Peak Detection)
  # -----------------------------------------------------------
  global_max_height <- max(d$y)
  
  # 寻找所有局部极大值 (Peaks)
  peak_indices <- which(diff(sign(diff(d$y))) == -2) + 1
  
  # 过滤出显著的峰 (高度 > 10% max)
  significant_peaks <- peak_indices[d$y[peak_indices] > (global_max_height * 0.1)]
  
  if (length(significant_peaks) == 0) return(d$x[1]) 
  
  # 取最靠右的显著峰作为主峰
  main_peak_idx <- significant_peaks[which.max(d$x[significant_peaks])]
  
  # 2. 准备数据：只分析主峰左侧
  # -----------------------------------------------------------
  # 我们需要计算直到主峰为止的二阶导数
  focus_indices <- 1:main_peak_idx
  
  if (length(focus_indices) < 5) return(d$x[1])
  
  focus_y <- d$y[focus_indices]
  
  # 3. 计算二阶导数
  # -----------------------------------------------------------
  d1 <- diff(focus_y)
  d2 <- diff(d1) 
  
  # d2 的长度比 focus_y 少 2
  # 我们需要对齐索引。d2[i] 对应的原始索引大约是 i+2
  
  # 4. 寻找二阶导数的局部峰值 (Local Peaks of 2nd Derivative)
  # -----------------------------------------------------------
  # 寻找 d2 曲线上的所有“山峰”
  d2_peak_indices_local <- which(diff(sign(diff(d2))) == -2) + 1
  
  # 过滤：只保留正值的峰 (d2 > 0 代表向上弯曲，d2 < 0 代表向下弯曲)
  # 我们只关心向上弯曲的拐角
  d2_peak_indices_local <- d2_peak_indices_local[d2[d2_peak_indices_local] > 0]
  
  if (length(d2_peak_indices_local) == 0) {
    # 如果二阶导数没有峰值（比如一直平滑上升），则回退到找 d2 最大值
    # 这通常发生在曲线非常平滑且单调的时候
    best_local_idx <- which.max(d2)
  } else {
    # 【核心逻辑】选择最后一个局部峰值 (Last Local Peak)
    # 因为 d2_peak_indices_local 是按索引从小到大排的，
    # 最后一个就是最靠近主峰 (最右边) 的那个 d2 峰值
    best_local_idx <- tail(d2_peak_indices_local, 1)
    
    # 可选优化：如果这个“最后一个峰”极小（噪音），而前面有个巨大的 d2 峰，
    # 可能需要加一个高度阈值判断。
    # 但按照你的要求“确保在正态分布内”，紧贴主峰的微小拐角也是有意义的，所以直接取最后一个。
  }
  
  # 5. 映射回原始 X 轴坐标
  cutoff_idx <- focus_indices[best_local_idx + 1] # +1 或 +2 进行索引对齐
  
  return(d$x[cutoff_idx])
}