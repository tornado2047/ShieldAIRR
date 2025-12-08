# 内部小工具：NULL 合并
`%||%` <- function(a, b) if (!is.null(a)) a else b

# 内部工具：大小写不敏感找列名
find_col <- function(dt, candidates) {
  cols <- names(dt)
  idx  <- which(tolower(cols) %in% tolower(candidates))
  if (length(idx) == 0) {
    stop("找不到列：", paste(candidates, collapse = " / "))
  }
  cols[idx[1]]
}

# ======== 基本统计量 ========

getSummaryStats <- function(df, locus = NULL) {
  data.table::setDT(df)
  cdr3_col  <- find_col(df, c("junction_aa", "cdr3_aa", "junctionaa"))
  clone_col <- find_col(df, c("clone_id", "duplicate_count", "cloneid"))

  total      <- nrow(df)
  unique_seq <- df[, .N, by = get(cdr3_col)][, .N]
  clones     <- df[, data.table::uniqueN(get(clone_col))]
  top_freq   <- max(table(df[[clone_col]])) / total

  data.frame(
    total_sequences      = total,
    unique_sequences     = unique_seq,
    clone_count          = clones,
    top_clone_frequency  = round(top_freq, 5),
    locus                = locus %||% "Unknown",
    stringsAsFactors     = FALSE
  )
}

# ======== CDR3 长度分布 ========

getLengthDistribution <- function(df, type = "aa") {
  data.table::setDT(df)
  col <- find_col(df,
                  if (type == "aa")
                    c("junction_aa", "cdr3_aa")
                  else
                    c("junction", "cdr3"))
  nchar(df[[col]])
}

plotLengthDistribution <- function(df,
                                   type = "aa",
                                   title = "CDR3 Length Distribution") {
  lens <- getLengthDistribution(df, type)
  len_df <- data.frame(Length = lens)

  ggplot2::ggplot(len_df, ggplot2::aes(x = Length)) +
    ggplot2::geom_histogram(
      bins = 35, fill = "#1b9e77",
      color = "white", alpha = 0.9, size = 0.3
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count)),
      color = "#0b5a47", size = 1.2, alpha = 0.8
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste(
        "Mean =", round(mean(lens), 2),
        "| Median =", stats::median(lens),
        "| Range =", paste(range(lens), collapse = "–")
      ),
      x = paste("CDR3 Length (", ifelse(type == "aa", "amino acids", "nucleotides"), ")"),
      y = "Count"
    ) +
    theme_shield(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ======== V / J 基因使用 ========

getVGeneDistributions <- function(df) {
  data.table::setDT(df)
  vcol <- find_col(df, c("v_call", "v.gene", "vcall", "v"))
  tab  <- df[, .(count = .N), keyby = .(gene = get(vcol))][order(-count)]
  tab[, frequency := count / sum(count) * 100][, count := NULL][]
}

getJGeneDistributions <- function(df) {
  data.table::setDT(df)
  jcol <- find_col(df, c("j_call", "j.gene", "jcall", "j"))
  tab  <- df[, .(count = .N), keyby = .(gene = get(jcol))][order(-count)]
  tab[, frequency := count / sum(count) * 100][, count := NULL][]
}

plotGeneUsage <- function(df, gene = "V", top_n = 25, title = NULL) {
  dat <- if (toupper(gene) == "V") getVGeneDistributions(df) else getJGeneDistributions(df)
  dat <- head(dat, top_n)
  if (is.null(title)) title <- paste("Top", top_n, gene, "Gene Usage")

  ggplot2::ggplot(dat, ggplot2::aes(x = reorder(gene, frequency), y = frequency)) +
    ggplot2::geom_col(
      fill = ifelse(toupper(gene) == "V", "#1b9e77", "#d95f02"),
      alpha = 0.92, color = "white", size = 0.3
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = paste("Total unique", gene, "genes:", nrow(dat)),
      x = paste(gene, "Gene"), y = "Frequency (%)"
    ) +
    theme_shield(base_size = 13) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10)
    )
}

# ======== V-J 配对矩阵 & 热图 ========

getVJPairing <- function(df) {
  data.table::setDT(df)
  vcol <- find_col(df, c("v_call", "v"))
  jcol <- find_col(df, c("j_call", "j"))
  tab  <- df[, .N, by = .(V = get(vcol), J = get(jcol))]
  mat  <- reshape2::dcast(tab, V ~ J, value.var = "N", fill = 0)
  mat  <- tibble::column_to_rownames(mat, "V")
  as.matrix(mat)
}

plotVJPairing <- function(df,
                          top_v = 20,
                          top_j = 35,
                          title = "V–J Gene Pairing Landscape") {
  mat <- getVJPairing(df)
  v_keep <- names(sort(rowSums(mat), decreasing = TRUE)[1:top_v])
  j_keep <- names(sort(colSums(mat), decreasing = TRUE)[1:top_j])
  mat_sub <- mat[v_keep, j_keep, drop = FALSE]
  melt_mat <- reshape2::melt(mat_sub,
                             varnames = c("V", "J"),
                             value.name = "Count")

  ggplot2::ggplot(melt_mat,
                  ggplot2::aes(
                    x = J,
                    y = factor(V, levels = rev(v_keep)),
                    fill = Count
                  )) +
    ggplot2::geom_tile(color = "white", size = 0.4) +
    ggplot2::scale_fill_viridis_c(option = "plasma", name = "Usage Count") +
    ggplot2::labs(
      title = title,
      subtitle = paste(top_v, "most used V ×", top_j, "most used J genes"),
      x = "J Gene", y = "V Gene"
    ) +
    theme_shield(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "right"
    )
}

# ======== 克隆丰度、多样性、Gini ========

getCloneAbundanceDistribution <- function(df) {
  data.table::setDT(df)
  clone_col <- find_col(df, c("clone_id", "duplicate_count", "clone"))
  as.vector(table(df[[clone_col]]))
}

getDiversityIndex <- function(df) {
  sizes <- getCloneAbundanceDistribution(df)
  N <- sum(sizes)
  if (N == 0) {
    return(data.frame(Index = character(), Value = numeric()))
  }

  p <- sizes / N
  p <- p[p > 0]
  shannon <- -sum(p * log(p))
  simpson <- 1 - sum(p^2)
  richness <- length(sizes)
  f1 <- sum(sizes == 1)
  f2 <- sum(sizes == 2)
  chao1 <- richness + ifelse(f2 > 0, f1 * (f1 - 1) / (2 * (f2 + 1)), f1 * (f1 - 1) / 2)

  data.frame(
    Index = c("Shannon", "Simpson", "Richness", "Chao1"),
    Value = c(shannon, simpson, richness, chao1),
    stringsAsFactors = FALSE
  )
}

getGiniIndex <- function(df) {
  sizes <- getCloneAbundanceDistribution(df)
  if (length(sizes) <= 1) return(0)
  sizes <- sort(sizes)
  n <- length(sizes)
  round(sum(sizes * (2 * (1:n) - n - 1)) / (n * sum(sizes)), 5)
}

plotCloneAbundance <- function(df,
                               title = "Clonal Expansion Rank–Abundance Curve") {
  abund <- sort(getCloneAbundanceDistribution(df), decreasing = TRUE)
  plot_df <- data.frame(Rank = seq_along(abund), Size = abund)

  stats <- getSummaryStats(df)
  top1_pct <- round(stats$top_clone_frequency * 100, 2)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Rank, y = Size)) +
    ggplot2::geom_line(color = "#6a3d9a", size = 1.3) +
    ggplot2::geom_point(color = "#6a3d9a", size = 1.2, alpha = 0.8) +
    ggplot2::scale_x_log10(labels = scales::comma_format()) +
    ggplot2::scale_y_log10(labels = scales::comma_format()) +
    ggplot2::labs(
      title = title,
      subtitle = paste(
        "Top clone =", scales::comma(abund[1]),
        paste0("(", top1_pct, "% of repertoire)"),
        "| Total clones =", scales::comma(length(abund))
      ),
      x = "Clone Rank (log10)", y = "Clone Size (log10)"
    ) +
    theme_shield(base_size = 14) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 11.5, color = "gray40")
    )
}

# ======== 一键 summary：对外导出的主函数 ========

#' Repertoire summary plot for a single AIRR sample
#'
#' 一键生成 CDR3 长度分布、V/J 使用、克隆丰度和 V–J 配对的综合图，并可保存 pdf/png。
#'
#' @param df 单个样本的数据框，至少包含 CDR3、V/J、clone_id/duplicate_count 等列。
#' @param sample_name 样本名，用于标题和输出文件名。
#' @param output_pdf 是否保存 PDF/PNG 文件，默认 TRUE。
#'
#' @return patchwork 拼接的 ggplot 对象（invisible 返回）。
#' @export
summarizeRepertoirePlot <- function(df,
                                    sample_name = "Sample",
                                    output_pdf = TRUE) {
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
        "| Clones:", scales::comma(stats$clone_count),
        "| Top clone:", round(stats$top_clone_frequency * 100, 2), "%"
      ),
      theme = theme_shield(base_size = 20)
    )

  if (output_pdf) {
    pdf_name <- paste0("Repertoire_Summary_", gsub("[ /]", "_", sample_name), ".pdf")
    ggplot2::ggsave(
      pdf_name, final_plot,
      width = 20, height = 16, dpi = 300,
      device = grDevices::cairo_pdf
    )
    png_name <- sub("\\.pdf$", ".png", pdf_name)
    ggplot2::ggsave(
      png_name, final_plot,
      width = 20, height = 16, dpi = 300
    )
    message("完整报告已保存：\n  ", pdf_name, "\n  ", png_name)
  }

  invisible(final_plot)
}
