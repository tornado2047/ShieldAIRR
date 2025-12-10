#' Repertoire summary helpers and plots
#'
#' 本文件包含针对单个 TCR/BCR AIRR 样本的统计量计算与可视化函数，
#' 包括 CDR3 长度分布、V/J 使用、V–J 配对、克隆丰度曲线以及综合概览图。
#'
#' 此文档块用于描述整个文件的功能，不对应单个函数。
#'
#' @name repertoire_summary
#' @keywords internal
#'
#' @import data.table
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom reshape2 dcast melt
#' @importFrom patchwork plot_annotation
#' @importFrom tibble column_to_rownames
#' @importFrom scales comma comma_format
NULL

# 内部小工具：NULL 合并（不导出）
`%||%` <- function(a, b) if (!is.null(a)) a else b

# 内部工具：大小写不敏感找列名（不导出）
find_col <- function(dt, candidates) {
  cols <- names(dt)
  idx  <- which(tolower(cols) %in% tolower(candidates))
  if (length(idx) == 0) {
    stop(paste("找不到列：", paste(candidates, collapse = " / ")))
  }
  cols[idx[1]]
}

# ---------------------------------------------------------------------
# 1. 基本统计量
# ---------------------------------------------------------------------

#' 计算 repertoire 基本统计量
#'
#' @param df data.frame 或 data.table，单个样本的 AIRR 数据。
#'   需包含 CDR3 列（例如 \code{junction_aa}）和克隆列
#'   （例如 \code{duplicate_count} 或 \code{clone_id}）。
#' @param locus 可选的位点标记，如 “TRA” / “TRB” 等。
#'
#' @return data.frame，包含总 reads 数、unique CDR3 数、克隆数、
#'   最高克隆比例等。
#' @export
getSummaryStats <- function(df, locus = NULL) {
  dt <- data.table::as.data.table(df)
  cdr3_col  <- find_col(dt, c("junction_aa", "cdr3_aa", "junctionaa"))
  clone_col <- find_col(dt, c("clone_id", "duplicate_count", "cloneid"))

  if (is.null(cdr3_col) || is.null(clone_col)) {
    stop("CDR3 column or clone column not found in the data.")
  }

  total      <- nrow(dt)
  # 按 CDR3 计数 unique 数量
  unique_seq <- dt[, .N, by = cdr3_col][, .N]
  clones     <- dt[, data.table::uniqueN(get(clone_col))]
  top_freq   <- max(table(dt[[clone_col]])) / total

  data.frame(
    total_sequences      = total,
    unique_sequences     = unique_seq,
    clone_count          = clones,
    top_clone_frequency  = round(top_freq, 5),
    locus                = locus %||% "Unknown",
    stringsAsFactors     = FALSE
  )
}

# ---------------------------------------------------------------------
# 2. CDR3 长度分布
# ---------------------------------------------------------------------

#' 获取 CDR3 长度分布
#'
#' @param df 单个样本 AIRR 数据。
#' @param type "aa" 表示氨基酸长度，"nt" 表示核苷酸长度。
#'
#' @return 整数向量，长度分布。
#' @export
getLengthDistribution <- function(df, type = "aa") {
  dt <- data.table::as.data.table(df)
  col <- find_col(
    dt,
    if (type == "aa") c("junction_aa", "cdr3_aa") else c("junction", "cdr3")
  )
  nchar(dt[[col]])
}

# ---------------------------------------------------------------------
# 3. V/J 基因使用频率
# ---------------------------------------------------------------------

#' 计算 V 基因使用频率
#'
#' @param df 单个样本 AIRR 数据。
#'
#' @return data.table，包含 \code{gene} 和 \code{frequency} 列。
#' @export
getVGeneDistributions <- function(df) {
  dt   <- data.table::as.data.table(df)
  vcol <- find_col(dt, c("v_call", "v.gene", "vcall", "v"))
  tab <- dt[
    ,
    list(count = .N),
    keyby = list(gene = get(vcol))
  ][order(-count)]
  tab[, frequency := count / sum(count) * 100]
  tab[, count := NULL]
}

#' 计算 J 基因使用频率
#'
#' @param df 单个样本 AIRR 数据。
#'
#' @return data.table，包含 \code{gene} 和 \code{frequency} 列。
#' @export
getJGeneDistributions <- function(df) {
  dt   <- data.table::as.data.table(df)
  jcol <- find_col(dt, c("j_call", "j.gene", "jcall", "j"))
  tab <- dt[
    ,
    list(count = .N),
    keyby = list(gene = get(jcol))
  ][order(-count)]
  tab[, frequency := count / sum(count) * 100]
  tab[, count := NULL]
}

# ---------------------------------------------------------------------
# 4. V–J 配对矩阵
# ---------------------------------------------------------------------

#' 计算 V–J 配对矩阵
#'
#' @param df 单个样本 AIRR 数据。
#'
#' @return 矩阵，行为 V 基因，列为 J 基因，元素为配对计数。
#' @export
getVJPairing <- function(df) {
  dt   <- data.table::as.data.table(df)
  vcol <- find_col(dt, c("v_call", "v"))
  jcol <- find_col(dt, c("j_call", "j"))
  tab <- dt[
    ,
    .N,
    by = list(V = get(vcol), J = get(jcol))
  ]
  mat <- reshape2::dcast(tab, V ~ J, value.var = "N", fill = 0)
  tibble::column_to_rownames(mat, "V") |> as.matrix()
}

# ---------------------------------------------------------------------
# 5. 克隆丰度与多样性
# ---------------------------------------------------------------------

#' 获取克隆丰度向量
#'
#' @param df 单个样本 AIRR 数据。
#'
#' @return 各克隆的 size 向量。
#' @export
getCloneAbundanceDistribution <- function(df) {
  dt        <- data.table::as.data.table(df)
  clone_col <- find_col(dt, c("clone_id", "duplicate_count", "clone"))
  as.vector(table(dt[[clone_col]]))
}

#' 计算多样性指数（Shannon / Simpson / Richness / Chao1）
#'
#' @param df 单个样本 AIRR 数据。
#'
#' @return data.frame，包含指数名称和取值。
#' @export
getDiversityIndex <- function(df) {
  sizes <- getCloneAbundanceDistribution(df)
  N <- sum(sizes)
  if (N == 0) {
    return(data.frame(Index = character(), Value = numeric()))
  }
  p <- sizes / N
  p <- p[p > 0]
  shannon  <- -sum(p * log(p))
  simpson  <- 1 - sum(p^2)
  richness <- length(sizes)
  f1 <- sum(sizes == 1)
  f2 <- sum(sizes == 2)
  chao1 <- richness +
    ifelse(
      f2 > 0,
      f1 * (f1 - 1) / (2 * (f2 + 1)),
      f1 * (f1 - 1) / 2
    )
  data.frame(
    Index = c("Shannon", "Simpson", "Richness", "Chao1"),
    Value = c(shannon, simpson, richness, chao1),
    stringsAsFactors = FALSE
  )
}

#' 计算 Gini 系数
#'
#' @param df 单个样本 AIRR 数据。
#'
#' @return 数值型，Gini 指数。
#' @export
getGiniIndex <- function(df) {
  sizes <- getCloneAbundanceDistribution(df)
  if (length(sizes) <= 1) return(0)
  sizes <- sort(sizes)
  n <- length(sizes)
  round(sum(sizes * (2 * (seq_len(n)) - n - 1)) / (n * sum(sizes)), 5)
}

# ---------------------------------------------------------------------
# 6. 可视化函数
# ---------------------------------------------------------------------

#' 绘制 CDR3 长度分布直方图
#'
#' @param df 单个样本 AIRR 数据。
#' @param type "aa" 或 "nt"。
#' @param title 图标题。
#'
#' @export
plotLengthDistribution <- function(df, type = "aa",
                                   title = "CDR3 Length Distribution") {
  lens <- getLengthDistribution(df, type)
  len_df <- data.frame(Length = lens)

  ggplot2::ggplot(len_df, ggplot2::aes(x = Length)) +
    ggplot2::geom_histogram(
      bins = 35,
      fill = "#1B9E77",
      color = "white"
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count)),
      color = "#0B5A47",
      linewidth = 1.2,
      alpha = 0.8
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste(
        "Mean =", round(mean(lens), 2),
        "| Median =", median(lens),
        "| Range =", paste(range(lens), collapse = "–")
      ),
      x = paste(
        "CDR3 Length (",
        ifelse(type == "aa", "amino acids", "nucleotides"),
        ")"
      ),
      y = "Count"
    ) +
    theme_shield(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' 绘制 V 或 J 基因使用条形图
#'
#' @param df 单个样本 AIRR 数据。
#' @param gene "V" 或 "J"。
#' @param top_n 显示前多少个基因。
#' @param title 可选标题。
#'
#' @export
plotGeneUsage <- function(df, gene = "V", top_n = 25, title = NULL) {

  dat <- if (toupper(gene) == "V") {
    getVGeneDistributions(df)
  } else {
    getJGeneDistributions(df)
  }

  dat <- head(dat, top_n)
  if (is.null(title)) title <- paste("Top", top_n, gene, "Gene Usage")

  ggplot2::ggplot(dat, ggplot2::aes(x = stats::reorder(gene, frequency),
                                    y = frequency)) +
    ggplot2::geom_col(
      fill  = ifelse(toupper(gene) == "V", "#1B9E77", "#D95F02"),
      alpha = 0.92,
      color = "white"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = title,
      subtitle = paste("Total unique", gene, "genes:", nrow(dat)),
      x        = paste(gene, "Gene"),
      y        = "Frequency (%)"
    ) +
    theme_shield(base_size = 13) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 10)
    )
}

#' 绘制 V–J 配对热图
#'
#' @param df 单个样本 AIRR 数据。
#' @param top_v 保留使用频率最高的 V 个数。
#' @param top_j 保留使用频率最高的 J 个数。
#' @param title 图标题。
#'
#' @export
plotVJPairing <- function(df, top_v = 20, top_j = 35,
                          title = "V–J Gene Pairing Landscape") {
  mat <- getVJPairing(df)
  v_keep <- names(sort(rowSums(mat), decreasing = TRUE)[seq_len(top_v)])
  j_keep <- names(sort(colSums(mat), decreasing = TRUE)[seq_len(top_j)])

  mat_sub <- mat[v_keep, j_keep, drop = FALSE]
  melt_mat <- reshape2::melt(
    mat_sub,
    varnames  = c("V", "J"),
    value.name = "Count"
  )

  ggplot2::ggplot(
    melt_mat,
    ggplot2::aes(
      x = J,
      y = factor(V, levels = rev(v_keep)),
      fill = Count
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(option = "plasma", name = "Usage Count") +
    ggplot2::labs(
      title    = title,
      subtitle = paste(
        top_v, "most used V ×", top_j, "most used J genes"
      ),
      x = "J Gene",
      y = "V Gene"
    ) +
    theme_shield(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size  = 9
      ),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "right"
    )
}

#' 绘制克隆丰度 rank–abundance 曲线
#'
#' @param df 单个样本 AIRR 数据。
#' @param title 图标题。
#'
#' @export
plotCloneAbundance <- function(df,
                               title = "Clonal Expansion Rank–Abundance Curve") {
  abund   <- sort(getCloneAbundanceDistribution(df), decreasing = TRUE)
  plot_df <- data.frame(Rank = seq_along(abund), Size = abund)

  stats <- getSummaryStats(df)
  top1_pct <- round(stats$top_clone_frequency * 100, 2)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Rank, y = Size)) +
    ggplot2::geom_line(color = "#6A3D9A", linewidth = 1.3) +
    ggplot2::geom_point(color = "#6A3D9A", size = 1.2, alpha = 0.8) +
    ggplot2::scale_x_log10(labels = scales::comma_format()) +
    ggplot2::scale_y_log10(labels = scales::comma_format()) +
    ggplot2::labs(
      title    = title,
      subtitle = paste(
        "Top clone =", scales::comma(abund[1]),
        paste0("(", top1_pct, "% of repertoire)"),
        "| Total clones =", scales::comma(length(abund))
      ),
      x = "Clone Rank (log10)",
      y = "Clone Size (log10)"
    ) +
    theme_shield(base_size = 14)
}

# ---------------------------------------------------------------------
# 7. 综合概览图
# ---------------------------------------------------------------------

#' 单样本 TCR/BCR repertoire 综合概览图
#'
#' @param df 单个样本 AIRR 数据。
#' @param sample_name 样本名称，用于标题。
#' @param output_pdf 是否自动保存 PDF/PNG 文件。
#'
#' @return 不可见返回 patchwork 对象 \code{final_plot}。
#' @export
summarizeRepertoirePlot <- function(df,
                                    sample_name = "Sample",
                                    output_pdf  = TRUE) {
  p1 <- plotLengthDistribution(df)
  p2 <- plotGeneUsage(df, "V", 25)
  p3 <- plotGeneUsage(df, "J", 35)
  p4 <- plotVJPairing(df, 20, 35)
  p5 <- plotCloneAbundance(df)

  stats <- getSummaryStats(df)

  final_plot <- (p1 | p2) /
    (p3 | p5) /
    p4 +
    patchwork::plot_annotation(
      title = "TCR/BCR Repertoire Summary",
      subtitle = paste(
        sample_name,
        "| Total reads:", scales::comma(stats$total_sequences),
        "| Unique CDR3:", scales::comma(stats$unique_sequences),
        "| Clones:",      scales::comma(stats$clone_count),
        "| Top clone:",   round(stats$top_clone_frequency * 100, 2), "%"
      ),
      theme = theme_shield(base_size = 20)
    )

  if (output_pdf) {
    pdf_name <- paste0(
      "Repertoire_Summary_",
      gsub("[ /]", "_", sample_name),
      ".pdf"
    )
    png_name <- sub("\\.pdf$", ".png", pdf_name)

    ggplot2::ggsave(pdf_name, final_plot,
                    width = 20, height = 16,
                    dpi = 300, device = grDevices::cairo_pdf)
    ggplot2::ggsave(png_name, final_plot,
                    width = 20, height = 16, dpi = 300)
  }

  invisible(final_plot)
}
