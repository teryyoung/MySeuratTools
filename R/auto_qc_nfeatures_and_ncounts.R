# ==============================================================================
# 主函数: 智能自动过滤 (Smart Auto Filter) - 更新版
# ==============================================================================
#' 自动识别分布形态并计算过滤阈值 (修复版)
#' @description
#' 这个函数用于根据 nFeature 或 nCount 或其他数值的分布（双峰或单峰长尾），
#' 自动选择合适的算法（谷底法或拐点法）来计算过滤阈值。
#' 双峰模式下，使用谷底作为阈值。
#' 单峰模式下，采用“二阶导数局部峰值法”，寻找分布函数的拐角。
#'
#' @param values 数值向量 (nFeature 或 nCount)。
#' @param filter_side 字符串。"low" (默认，切除左尾) 或 "high" (切除右尾)。
#' @param knee_selection 字符串。"closer" (默认，选离主峰最近的拐点) 或 "further" (选离主峰最远的拐点)。
#' @param adjust 平滑度 (默认 2.0)。
#' @param assumption "gaussian" 或 "poisson" (用于 counts 的 sqrt 变换)。
#' @param min_peak_ratio 双峰检测时，副峰的最小高度比例 (默认 1e-4)。
#' @param valley_strictness 双峰检测时，谷底深度的严格程度 (默认 0.99，即谷底低于某一峰值的0.99倍则判定为谷底)。
#' 
#' @export
smart_qc_cutoff <- function(values, 
                            filter_side = c("low", "high"), 
                            knee_selection = c("closer", "further"),
                            adjust = 2.0, 
                            assumption = "gaussian", 
                            min_peak_ratio = 1e-4, 
                            valley_strictness = 0.99,
                            show_plot = TRUE, 
                            title = "") {
  
  # 参数校验
  filter_side <- match.arg(filter_side)
  knee_selection <- match.arg(knee_selection)
  
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
  
  # 4. 寻找主峰 (根据 filter_side 自动判断左右)
  main_peak_idx <- .find_main_peak(d, max_indices, filter_side,)
  
  # 如果连主峰都找不到（极罕见），直接返回边缘
  if (is.null(main_peak_idx)) {
    warning("Cannot find main peak.")
    return(ifelse(filter_side == "low", min(raw_values), max(raw_values)))
  }
  
  # 5. 双峰检测 (Valley Method)
  valley_res <- .get_valley_cutoff(d, min_indices, max_indices, main_peak_idx, filter_side, min_peak_ratio)
  
  method_used <- "Unknown"
  calc_cutoff <- NA
  is_bimodal <- FALSE
  
  if (!is.null(valley_res)) {
    if (valley_res$ratio < valley_strictness) {
      is_bimodal <- TRUE
      method_used <- "Bimodal (Valley)"
      calc_cutoff <- valley_res$cutoff
    } else {
      method_used <- "Unimodal (Knee) - Valley too shallow"
    }
  } else {
    method_used <- "Unimodal (Knee) - No secondary peak"
  }
  
  # 6. 单峰检测 (Curvature/Knee Method)
  if (!is_bimodal) {
    calc_cutoff <- .get_curvature_cutoff(d, main_peak_idx, filter_side, knee_selection)
    
    # Fallback 逻辑: 如果单峰法切得太离谱，且存在浅谷底，回退到浅谷底
    # 离谱定义：
    # Low模式: Cutoff 极小 (<10)
    # High模式: Cutoff 极大 (接近 max)
    if (!is.null(valley_res)) {
       if ((filter_side == "low" && calc_cutoff < 10) || 
           (filter_side == "high" && calc_cutoff > max(d$x)*0.95)) {
         calc_cutoff <- valley_res$cutoff
         method_used <- "Bimodal (Shallow Valley Fallback)"
         is_bimodal <- TRUE
       }
    }
    
    if (method_used == "Unknown") {
      method_used <- paste0("Unimodal (Knee: ", knee_selection, ")")
    }
  }
  
  # 7. 还原阈值
  final_cutoff <- calc_cutoff
  if (assumption == "poisson") {
    final_cutoff <- calc_cutoff^2
    if (filter_side == "low" && final_cutoff < 0) final_cutoff <- 0
  }
  
  # 强制边界保护
  if (filter_side == "low") {
    final_cutoff <- max(0, final_cutoff)
  }
  
  # 8. 绘图
  if (show_plot) {
    d_plot <- density(raw_values, adjust = adjust) 
    df <- data.frame(x = d_plot$x, y = d_plot$y)
    
    # 确定主峰X坐标用于画图
    main_peak_x_plot <- d_plot$x[which.min(abs(d_plot$x - ifelse(assumption=="poisson", d$x[main_peak_idx]^2, d$x[main_peak_idx])))]
    
    p <- ggplot(df, aes(x = x)) +
      geom_area(aes(y = y), fill = ifelse(is_bimodal, "#a6cee3", "#fdbf6f"), alpha = 0.5) +
      geom_line(aes(y = y), size = 1, color = "black") +
      
      # 画出阈值线
      geom_vline(xintercept = final_cutoff, color = "red", linetype = "dashed", size = 1) +
      
      # 标出主峰位置 (蓝色点线)
      geom_vline(xintercept = main_peak_x_plot, color = "blue", linetype = "dotted", alpha=0.5) +
      
      annotate("text", x = final_cutoff, y = max(d_plot$y)*0.6, 
               label = paste0(filter_side, "-Cut: ", round(final_cutoff, 0)), 
               color = "red", angle = 90, vjust = -0.5) +
      
      theme_classic() +
      labs(title = paste0(title, " (", filter_side, " filter)"),
           subtitle = paste0("Method: ", method_used, "\nKnee Selection: ", knee_selection),
           x = "Value", y = "Density")
    
    print(p)
  }
  
  message(paste0("[", title, "] Side: ", filter_side, " | Mode: ", method_used, " | Threshold: ", round(final_cutoff, 2)))
  return(final_cutoff)
}


# ==============================================================================
# 内部工具: 寻找主峰 (根据过滤方向自动调整)
# ==============================================================================
.find_main_peak <- function(d, max_indices, filter_side) {
  global_max_height <- max(d$y)
  # 寻找显著峰 (高度 > 10% 全局最大值)
  significant_peaks <- max_indices[d$y[max_indices] > (global_max_height * 0.1)]
  
  if (length(significant_peaks) == 0) return(NULL)
  
  if (filter_side == "low") {
    # 低位过滤：假设正常细胞在右边，垃圾在左边。
    # 所以主峰是【最靠右】的那个显著峰。
    return(significant_peaks[which.max(d$x[significant_peaks])])
  } else {
    # 高位过滤：假设正常细胞在左边，双细胞在右边。
    # 所以主峰是【最靠左】的那个显著峰。
    return(significant_peaks[which.min(d$x[significant_peaks])])
  }
}

# ==============================================================================
# 内部工具: 寻找二阶导数拐点 (支持方向和远近选择)
# ==============================================================================
.get_curvature_cutoff <- function(d, main_peak_idx, filter_side, knee_selection) {
  
  # 1. 确定分析范围 (Focus Area)
  if (filter_side == "low") {
    # 分析：从最左端 -> 主峰
    focus_indices <- 1:main_peak_idx
  } else {
    # 分析：从主峰 -> 最右端
    focus_indices <- main_peak_idx:length(d$x)
  }
  
  if (length(focus_indices) < 5) return(d$x[main_peak_idx])
  
  focus_y <- d$y[focus_indices]
  
  # 2. 计算二阶导数
  d1 <- diff(focus_y)
  d2 <- diff(d1)
  
  # 3. 寻找 d2 的局部峰值 (Local Peaks of 2nd Derivative)
  # 注意：d2 的索引相对于 focus_indices 有偏移
  d2_peak_indices_local <- which(diff(sign(diff(d2))) == -2) + 1
  
  # 只保留正向弯曲 (Convexity) 的点
  d2_peak_indices_local <- d2_peak_indices_local[d2[d2_peak_indices_local] > 0]
  
  best_local_idx <- NULL
  
  if (length(d2_peak_indices_local) == 0) {
    # 如果没有明显的峰，取 d2 最大值
    best_local_idx <- which.max(d2)
  } else {
    # 【核心逻辑】根据用户选择挑选拐点
    
    # 定义：Focus 向量在 Low 模式下是 [尾巴 -> 峰]，在 High 模式下是 [峰 -> 尾巴]
    # 我们需要根据物理距离来判断
    
    if (filter_side == "low") {
      # 向量方向: [0, 10, 50, ... Peak]
      # index 越大，离 Peak 越近
      if (knee_selection == "closer") {
        best_local_idx <- tail(d2_peak_indices_local, 1) # 选 index 最大的 (近峰)
      } else {
        best_local_idx <- head(d2_peak_indices_local, 1) # 选 index 最小的 (远峰/近尾)
      }
    } else {
      # 向量方向: [Peak, ... 5000, 10000]
      # index 越小，离 Peak 越近
      if (knee_selection == "closer") {
        best_local_idx <- head(d2_peak_indices_local, 1) # 选 index 最小的 (近峰)
      } else {
        best_local_idx <- tail(d2_peak_indices_local, 1) # 选 index 最大的 (远峰/近尾)
      }
    }
  }
  
  # 4. 映射回原始坐标
  # d2 索引修正：diff 两次，位置大约偏移 +1
  cutoff_idx <- focus_indices[best_local_idx + 1]
  
  return(d$x[cutoff_idx])
}

# ==============================================================================
# 内部工具: 寻找谷底 (支持双向)
# ==============================================================================
.get_valley_cutoff <- function(d, min_indices, max_indices, main_peak_idx, filter_side, min_peak_ratio) {
  
  global_max_height <- max(d$y)
  
  # 1. 确定搜索区域的潜在副峰
  if (filter_side == "low") {
    # 找主峰左边的
    potential_sec_peaks <- max_indices[max_indices < main_peak_idx]
  } else {
    # 找主峰右边的
    potential_sec_peaks <- max_indices[max_indices > main_peak_idx]
  }
  
  # 过滤噪音峰
  sec_peaks <- potential_sec_peaks[d$y[potential_sec_peaks] > (global_max_height * min_peak_ratio)]
  
  if (length(sec_peaks) == 0) return(NULL)
  
  # 2. 遍历寻找最佳谷底
  best_result <- NULL
  best_ratio <- 1.0
  
  for (sec_idx in sec_peaks) {
    # 确定谷底搜索范围
    if (filter_side == "low") {
      valid_valleys <- min_indices[min_indices > sec_idx & min_indices < main_peak_idx]
    } else {
      valid_valleys <- min_indices[min_indices > main_peak_idx & min_indices < sec_idx]
    }
    
    if (length(valid_valleys) > 0) {
      # 找这一段里最低的点
      valley_idx <- valid_valleys[which.min(d$y[valid_valleys])]
      
      valley_h <- d$y[valley_idx]
      sec_peak_h <- d$y[sec_idx]
      current_ratio <- valley_h / sec_peak_h
      
      if (current_ratio < best_ratio) {
        best_ratio <- current_ratio
        best_result <- list(cutoff = d$x[valley_idx], ratio = current_ratio)
      }
    }
  }
  
  return(best_result)
}