#' 绘制分割的高亮聚类散点图
#' @description
#' 这个函数用于使用Seurat包的DimPlot函数绘制按照group_by分割的降维散点图。
#' 主要更新在于分割的散点图包含整体灰色背景
#'
#' @param object 一个 Seurat 对象。
#' @param group_by 字符型。元数据中用于分组的列名（例如 "seurat_clusters" 或 "orig.ident"）。
#' @param reduction 字符型。降维方式，默认为 "umap"（你代码中是 "corrected_umap"）。
#' @param colors 字符向量（可选）。自定义颜色向量。如果为空，则自动生成默认颜色。
#' @param pt_size 数值型。背景点的大小。
#' @param highlight_size 数值型。高亮点的大小。
#' @param ncol 整数型。拼图时的列数。
#' @param save_path 字符型（可选）。如果提供路径（以.png或.pdf结尾），则自动保存图片。
#' @param width 数值型。保存图片的宽度。
#' @param height 数值型。保存图片的高度。
#'
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @importFrom scales hue_pal
#'
#' @return 一个 ggplot/patchwork 对象。
#' @export
PlotSplitHighlight <- function(object,
                               group_by = "seurat_clusters",
                               reduction = "umap",
                               colors = NULL,
                               pt_size = 0.2,
                               highlight_size = 0.2,
                               ncol = 2,
                               save_path = NULL,
                               width = 8,
                               height = 12) {

  # 1. 检查 group_by 是否存在于 metadata 中
  if (!group_by %in% colnames(object@meta.data)) {
    stop(paste("Error: Column", group_by, "not found in meta.data"))
  }

  # 获取分组信息并排序
  # 确保转换为因子或字符，防止因NULL报错
  split_groups <- sort(unique(object@meta.data[[group_by]]))

  # 2. 定义颜色
  if (is.null(colors)) {
    # 如果用户没提供颜色，使用默认 hue_pal
    my_colors <- scales::hue_pal()(length(split_groups))
  } else {
    # 如果提供了颜色，确保长度足够
    if (length(colors) < length(split_groups)) {
      warning("提供的颜色数量少于分组数量，将循环使用颜色。")
    }
    my_colors <- colors
  }
  # 给颜色向量命名，确保一一对应
  names(my_colors) <- split_groups

  # 3. 使用 lapply 循环生成子图
  plot_list <- lapply(split_groups, function(current_group) {

    # 获取属于当前组的细胞ID (使用基础R筛选，比WhichCells expression更适合函数内部)
    cells_target <- rownames(object@meta.data)[object@meta.data[[group_by]] == current_group]

    # 绘图
    p <- Seurat::DimPlot(object,
                         reduction = reduction,
                         cells.highlight = cells_target,
                         cols.highlight = my_colors[[as.character(current_group)]], # 提取特定颜色
                         cols = "grey85",
                         pt.size = pt_size,
                         sizes.highlight = highlight_size
    ) +
      ggplot2::ggtitle(current_group) +
      Seurat::NoLegend() +
      ggplot2::theme(axis.title = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5))

    return(p)
  })

  # 4. 拼图
  # 获取整体标题
  main_title <- paste("Split by:", group_by)

  final_plot <- patchwork::wrap_plots(plot_list, ncol = ncol) +
    patchwork::plot_annotation(title = main_title,
                               theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5)))

  # 5. 保存 (如果提供了路径)
  if (!is.null(save_path)) {
    ggplot2::ggsave(plot = final_plot, filename = save_path, width = width, height = height)
    message(paste("Plot saved to:", save_path))
  }

  return(final_plot)
}

#' @description
#'
#' 绘制多基因共表达情况（柱状图 + UpSet图）
#'
#' @param obj Seurat 对象。
#' @param genes 字符向量。包含 2-6 个基因名称。
#' @param threshold 数值。判定表达阳性的阈值，默认为 0。
#' @param out_dir 字符。输出目录路径。
#' @param plot_prefix 字符。保存文件名的前缀（可选）。
#'
#' @import Seurat
#' @import ggplot2
#' @import ComplexHeatmap
#' @importFrom grid unit
#' @importFrom grDevices png dev.off
#'
#' @return 返回一个列表，包含柱状图对象(ggplot)和UpSet数据对象。
#' @export
PlotGeneCoexpression <- function(obj,
                                 genes,
                                 threshold = 0,
                                 out_dir = "./result/coexpression/",
                                 plot_prefix = NULL) {

  # --- 0. 输入检查 ---

  # 检查基因数量 (2-6个)
  if (length(genes) < 2 | length(genes) > 6) {
    stop("Error: 请提供 2 到 6 个基因名称！")
  }

  # 检查基因是否都在 Seurat 对象中
  # 注意：这里默认检查 RNA assay，如果用 integrated 或 SCT，请根据需要调整
  available_genes <- rownames(obj)
  missing_genes <- genes[!genes %in% available_genes]
  if (length(missing_genes) > 0) {
    stop(paste("Error: 以下基因在 Seurat 对象中未找到:", paste(missing_genes, collapse = ", ")))
  }

  # 检查/创建目录
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # 定义文件名前缀
  if (is.null(plot_prefix)) {
    plot_prefix <- paste(genes, collapse = "_")
  }

  # --- 1. 提取数据 ---

  # 提取表达矩阵
  df_expr <- Seurat::FetchData(obj, vars = genes)

  # --- 2. 生成组合标签 (Barplot 数据准备) ---

  # 定义一个辅助函数来生成标签
  get_status_label <- function(row_values, gene_names, thresh) {
    idx <- which(row_values > thresh)
    if (length(idx) == 0) {
      return("All Negative") # 或者 "None"
    } else {
      return(paste0(gene_names[idx], collapse = " + "))
    }
  }

  # 应用辅助函数 (使用 apply 逐行处理)
  df_expr$Status <- apply(df_expr[, genes], 1, function(x) get_status_label(x, genes, threshold))

  # 统计频数
  plot_data <- as.data.frame(table(df_expr$Status))
  colnames(plot_data) <- c("Combination", "Count")

  # 排序：数量多的排前面
  plot_data <- plot_data[order(plot_data$Count, decreasing = TRUE), ]
  plot_data$Combination <- factor(plot_data$Combination, levels = plot_data$Combination)

  # --- 3. 绘制柱状图 (ggplot2) ---

  p_bar <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Combination, y = Count)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.2, size = 3) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Cell Counts by Gene Combinations",
                  subtitle = paste("Genes:", paste(genes, collapse=", "), "\nThreshold >", threshold),
                  x = "Expression Combination",
                  y = "Number of Cells") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  # 保存柱状图 (ggsave 适用于 ggplot 对象)
  bar_filename <- file.path(out_dir, paste0(plot_prefix, "_barplot.png"))
  ggplot2::ggsave(filename = bar_filename, plot = p_bar, width = 10, height = 6)
  message(paste("Barplot saved to:", bar_filename))

  # --- 4. 绘制 UpSet 图 (ComplexHeatmap) ---

  # 准备二值化矩阵
  binary_mat <- as.matrix(ifelse(df_expr[, genes] > threshold, 1, 0))

  # 制作组合矩阵对象
  # suppressMessages 用于屏蔽 ComplexHeatmap 可能产生的大量提示信息
  comb_mat <- suppressMessages(ComplexHeatmap::make_comb_mat(binary_mat))

  # 生成 UpSet 图对象
  # 注意：UpSet() 返回的是一个 HeatmapList 对象，不是 ggplot 对象
  upset_plot <- ComplexHeatmap::UpSet(comb_mat,
                                      pt_size = grid::unit(3, "mm"),
                                      lwd = 2,
                                      comb_col = "black",
                                      # 顶部柱状图设置
                                      top_annotation = ComplexHeatmap::upset_top_annotation(
                                        comb_mat,
                                        add_numbers = TRUE,
                                        height = grid::unit(6, "cm") # 动态调整高度
                                      ),
                                      # 底部点阵图设置 (根据基因数量动态微调高度)
                                      height = grid::unit(length(genes) * 1, "cm"),
                                      column_title = paste("UpSet Plot (Threshold >", threshold, ")")
  )

  # 保存 UpSet 图 (关键修改：不能用 ggsave)
  upset_filename <- file.path(out_dir, paste0(plot_prefix, "_upset.png"))

  # 打开图形设备
  png(filename = upset_filename, width = 8, height = 6, units = "in", res = 300)
  # 绘制对象
  ComplexHeatmap::draw(upset_plot)
  # 关闭设备
  dev.off()

  message(paste("UpSet plot saved to:", upset_filename))

  # 返回结果以便后续查看
  return(list(barplot = p_bar, upset_data = comb_mat))
}

