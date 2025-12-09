#' CDR3 physicochemical landscape analysis
#'
#' 给定 AIRR 风格的数据框（含 \code{junction_aa} / \code{duplicate_count}），
#' 计算 CDR3 的多种理化性质、Atchley/Kidera 因子，并拼成综合图。
#'
#' 依赖外部包 \code{sumrep}，该包未在 CRAN 上发布，需要用户自行安装或
#' 在当前会话中通过 \code{devtools::load_all()} 加载本地版本。
#'
#' @param df data.frame，至少包含列 \code{junction_aa} 和 \code{duplicate_count}。
#' @param sample_name 字符串，用于图标题和输出。
#' @param output_prefix 输出文件前缀（不含扩展名）。如果为 \code{NULL} 则不保存文件。
#'
#' @return patchwork 拼接的 ggplot 对象（invisible 返回）。
#' @export
shield_cdr3_landscape <- function(df,
                                  sample_name   = "Sample",
                                  output_prefix = "CDR3_Physicochemical_Landscape") {
  # -------- 0. 检查 sumrep 是否可用 --------
  if (!requireNamespace("sumrep", quietly = TRUE)) {
    stop(
      "本函数需要额外安装 'sumrep' 包。\n",
      "当前会话找不到 sumrep，请先在 R 里运行例如：\n\n",
      "  devtools::load_all('/Users/xfcheung/workspace/AIDeN/sumrep/')\n\n",
      "或者将 sumrep 安装成正式包后再调用 shield_cdr3_landscape()。",
      call. = FALSE
    )
  }
  # -------- 1. 提取并清洗 CDR3 序列 --------
  if (!("junction_aa" %in% names(df))) {
    stop("df 中找不到列 'junction_aa'，请确认数据格式。")
  }
  if (!("duplicate_count" %in% names(df))) {
    stop("df 中找不到列 'duplicate_count'，请确认数据格式。")
  }
  seqs <- df$junction_aa
  seqs <- seqs[!is.na(seqs) & seqs != "" & nchar(seqs) > 0]
  dup_count <- df$duplicate_count[match(seqs, df$junction_aa)]
  stopifnot(length(seqs) == length(dup_count))
  # -------- 2. 超安全理化性质计算函数 --------
  safe_calc <- function(func, seq_vector) {
    vapply(
      seq_vector,
      function(s) {
        s <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(s))
        if (nchar(s) < 2) return(NA_real_)
        res <- tryCatch(func(s), error = function(e) NA_real_)
        if (is.null(res) || length(res) == 0 || is.function(res)) return(NA_real_)
        as.numeric(res)
      },
      numeric(1)
    )
  }
  # -------- 3. 逐条计算理化性质 --------
  phys <- data.frame(
    Sequence      = seqs,
    Length        = nchar(seqs),
    Duplicate     = dup_count,
    GRAVY         = safe_calc(Peptides::gravy, seqs),
    Charge_pH7    = safe_calc(function(x) Peptides::charge(x, pH = 7), seqs),
    pI            = safe_calc(Peptides::pI, seqs),
    Aliphatic     = safe_calc(Peptides::aIndex, seqs),
    Hydrophob     = safe_calc(Peptides::hydrophobicity, seqs),
    Boman         = safe_calc(Peptides::boman, seqs),
    Aromatic_prop = vapply(
      seqs,
      function(x) {
        aa <- strsplit(x, "")[[1]]
        mean(aa %in% c("F", "Y", "W"))
      },
      numeric(1)
    ),
    Basic_prop = vapply(
      seqs,
      function(x) {
        aa <- strsplit(x, "")[[1]]
        mean(aa %in% c("K", "R", "H"))
      },
      numeric(1)
    ),
    Acidic_prop = vapply(
      seqs,
      function(x) {
        aa <- strsplit(x, "")[[1]]
        mean(aa %in% c("D", "E"))
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  phys$logDup <- log10(phys$Duplicate + 1)
  # -------- 4. Atchley 因子分布 --------
  AF <- sumrep::getAtchleyFactorDistributions(df)
  AF_df <- dplyr::bind_rows(
    lapply(seq_along(AF), function(i) {
      factor_name <- if (!is.null(names(AF)[i]) && names(AF)[i] != "") {
        names(AF)[i]
      } else {
        paste0("Factor", i)
      }
      data.frame(
        Factor = factor_name,
        Value  = as.numeric(AF[[i]]),
        stringsAsFactors = FALSE
      )
    })
  )
  p_atchley <- ggplot2::ggplot(AF_df, ggplot2::aes(x = Factor, y = Value)) +
    ggplot2::geom_violin(
      fill = "#1B9E77",
      color = "white",
      alpha = 0.85,
      scale = "width"
    ) +
    ggplot2::geom_boxplot(width = 0.2, outlier.alpha = 0) +
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size  = 3,
      fill  = "white"
    ) +
    ggplot2::labs(
      title = "Atchley Factors (5D)",
      x = "",
      y = "Score"
    ) +
    theme_shield(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  # -------- 5. Kidera 因子 + PCA --------
  KF <- sumrep::getKideraFactorDistributions(df)
  KF_mat <- do.call(cbind, KF)
  pca_kf <- stats::prcomp(KF_mat, scale. = TRUE)
  pca_df <- data.frame(
    PC1       = pca_kf$x[, 1],
    PC2       = pca_kf$x[, 2],
    Duplicate = dup_count,
    logDup    = log10(dup_count + 1)
  )
  var_exp <- summary(pca_kf)$importance[2, ]  # 第二行是 Proportion of Variance
  p_pca <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2, color = logDup, size = Duplicate)
  ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_color_viridis_c(name = "log10(Clone size + 1)") +
    ggplot2::scale_size_continuous(range = c(0.5, 6)) +
    ggplot2::labs(
      title    = "PCA of Kidera Factors (10D)",
      subtitle = paste0(
        "PC1: ", round(100 * var_exp[1], 1), "% | ",
        "PC2: ", round(100 * var_exp[2], 1), "%"
      ),
      x = "PC1",
      y = "PC2"
    ) +
    theme_shield(base_size = 14)
  # -------- 6. Aliphatic index vs 克隆扩增 --------
  p_ai <- ggplot2::ggplot(phys, ggplot2::aes(x = Aliphatic, y = Duplicate)) +
    ggplot2::geom_point(
      ggplot2::aes(color = logDup, size = Duplicate),
      alpha = 0.7
    ) +
    ggplot2::geom_smooth(
      method   = "gam",
      color    = "red",
      se       = FALSE,
      linetype = "dashed"
    ) +
    ggplot2::scale_y_log10(labels = scales::comma) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    ggplot2::labs(
      title = "Aliphatic Index vs Clone Expansion",
      x     = "Aliphatic Index",
      y     = "Clone size (log10)"
    ) +
    theme_shield(base_size = 14)
  # -------- 7. 关键理化性质的密度分布 --------
  phys_long <- phys |>
    dplyr::select(
      Length, GRAVY, Charge_pH7, pI,
      Aliphatic, Aromatic_prop, Basic_prop, Acidic_prop
    ) |>
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "Property",
      values_to = "Value"
    )
  p_dist <- ggplot2::ggplot(phys_long, ggplot2::aes(x = Value)) +
    ggplot2::geom_density(
      fill  = "#66C2A5",
      alpha = 0.7,
      color = "black"
    ) +
    ggplot2::facet_wrap(~Property, scales = "free", ncol = 4) +
    ggplot2::labs(
      title = "Key Physicochemical Properties Distribution",
      x     = NULL,
      y     = "Density"
    ) +
    ggplot2::theme_bw(base_size = 11)
  # -------- 8. Top 50 clone 理化性质热图（独立图） --------
  if (nrow(phys) >= 1) {
    top_n <- min(50L, nrow(phys))
    top50 <- phys |>
      dplyr::arrange(dplyr::desc(Duplicate)) |>
      dplyr::slice_head(n = top_n)
    mat <- top50 |>
      dplyr::select(
        GRAVY, Charge_pH7, pI, Aliphatic,
        Aromatic_prop, Basic_prop, Acidic_prop
      ) |>
      scale(center = TRUE, scale = TRUE) |>
      as.matrix()
    rownames(mat) <- paste0(
      "Clone", seq_len(top_n), " (n=",
      top50$Duplicate, ")"
    )
    pheatmap::pheatmap(
      mat,
      color          = grDevices::colorRampPalette(
        c("#2C7BB6", "white", "#D7191C")
      )(100),
      cluster_rows   = TRUE,
      cluster_cols   = TRUE,
      show_rownames  = TRUE,
      fontsize_row   = 7,
      main           = "Top Expanded Clones – Physicochemical Z-scores",
      device         = ragg::agg_png
    )
  }
  # -------- 9. 拼成主图并可选保存 --------
  main_title <- "Comprehensive Physicochemical Landscape of TCR/BCR CDR3 Repertoire"
  subtitle   <- paste(
    sample_name,
    "| Unique sequences:", nrow(phys),
    "| Total reads:", scales::comma(sum(phys$Duplicate))
  )
  final_plot <- (p_atchley | p_pca) /
    (p_ai | p_dist) +
    patchwork::plot_annotation(
      title  = main_title,
      subtitle = subtitle,
      theme  = theme_shield(base_size = 18)
    )
  if (!is.null(output_prefix)) {
    pdf_file <- paste0(output_prefix, ".pdf")
    png_file <- paste0(output_prefix, ".png")
    ggplot2::ggsave(
      filename = pdf_file,
      plot     = final_plot,
      width    = 18,
      height   = 13,
      dpi      = 300,
      device   = grDevices::cairo_pdf
    )
    ggplot2::ggsave(
      filename = png_file,
      plot     = final_plot,
      width    = 18,
      height   = 13,
      dpi      = 300
    )
  }
  invisible(final_plot)
}
