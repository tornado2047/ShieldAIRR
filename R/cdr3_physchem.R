#' CDR3 physicochemical landscape analysis
#'
#' 给定 AIRR 风格的数据框（含 junction_aa / duplicate_count），
#' 计算 CDR3 的理化性质、Atchley/Kidera 因子，并输出综合图。
#'
#' @param df data.frame，至少包含列 \code{junction_aa} 和 \code{duplicate_count}。
#' @param sample_name 字符串，用于图标题和输出。
#' @param output_prefix 输出文件前缀（不含扩展名）。如果为 NULL 则不保存文件。
#'
#' @return patchwork 拼接的 ggplot 对象（invisible 返回）。
#' @export
shield_cdr3_landscape <- function(df,
                                  sample_name = "Sample",
                                  output_prefix = "CDR3_Physicochemical_Landscape") {
  seqs <- df$junction_aa
  seqs <- seqs[!is.na(seqs) & seqs != "" & nchar(seqs) > 0]

  dup_count <- df$duplicate_count[match(seqs, df$junction_aa)]
  stopifnot(length(seqs) == length(dup_count))

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

  phys <- data.frame(
    Sequence     = seqs,
    Length       = nchar(seqs),
    Duplicate    = dup_count,
    GRAVY        = safe_calc(Peptides::gravy, seqs),
    Charge_pH7   = safe_calc(function(x) Peptides::charge(x, pH = 7), seqs),
    pI           = safe_calc(Peptides::pI, seqs),
    Aliphatic    = safe_calc(Peptides::aIndex, seqs),
    Hydrophob    = safe_calc(Peptides::hydrophobicity, seqs),
    Boman        = safe_calc(Peptides::boman, seqs),
    Aromatic_prop = vapply(seqs, function(x) {
      aa <- strsplit(x, "")[[1]]
      mean(aa %in% c("F", "Y", "W"))
    }, numeric(1)),
    Basic_prop   = vapply(seqs, function(x) {
      aa <- strsplit(x, "")[[1]]
      mean(aa %in% c("K", "R", "H"))
    }, numeric(1)),
    Acidic_prop  = vapply(seqs, function(x) {
      aa <- strsplit(x, "")[[1]]
      mean(aa %in% c("D", "E"))
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  phys$logDup <- log10(phys$Duplicate + 1)

  # Atchley factors
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
      fill = "#1b9e77", color = "white",
      alpha = 0.85, scale = "width"
    ) +
    ggplot2::geom_boxplot(width = 0.2, outlier.alpha = 0) +
    ggplot2::stat_summary(
      fun = mean, geom = "point",
      shape = 23, size = 3, fill = "white"
    ) +
    theme_shield(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "Atchley Factors (5D)",
      x = "", y = "Score"
    )

  # Kidera factors + PCA
  KF <- sumrep::getKideraFactorDistributions(df)
  KF_mat <- do.call(cbind, KF)
  pca_kf <- stats::prcomp(KF_mat, scale. = TRUE)

  pca_df <- data.frame(
    PC1       = pca_kf$x[, 1],
    PC2       = pca_kf$x[, 2],
    Duplicate = dup_count,
    logDup    = log10(dup_count + 1)
  )

  p_pca <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2, color = logDup, size = Duplicate)
  ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_color_viridis_c(name = "log10(Clone size + 1)") +
    ggplot2::scale_size_continuous(range = c(0.5, 6)) +
    theme_shield(base_size = 14) +
    ggplot2::labs(
      title = "PCA of Kidera Factors (10D)",
      subtitle = paste0(
        "PC1: ", round(100 * summary(pca_kf)$importance[2, 1], 1), "% | ",
        "PC2: ", round(100 * summary(pca_kf)$importance[2, 2], 1), "%"
      )
    )

  # Aliphatic vs abundance
  p_ai <- ggplot2::ggplot(phys, ggplot2::aes(x = Aliphatic, y = Duplicate)) +
    ggplot2::geom_point(
      ggplot2::aes(color = logDup, size = Duplicate),
      alpha = 0.7
    ) +
    ggplot2::geom_smooth(
      method = "gam", color = "red",
      se = FALSE, linetype = "dashed"
    ) +
    ggplot2::scale_y_log10(labels = scales::comma) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    theme_shield(base_size = 14) +
    ggplot2::labs(
      title = "Aliphatic Index vs Clone Expansion",
      x = "Aliphatic Index", y = "Clone size (log10)"
    )

  # Key properties distribution
  phys_long <- phys |>
    dplyr::select(
      Length, GRAVY, Charge_pH7, pI,
      Aliphatic, Aromatic_prop, Basic_prop, Acidic_prop
    ) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "Property", values_to = "Value"
    )

  p_dist <- ggplot2::ggplot(phys_long, ggplot2::aes(x = Value)) +
    ggplot2::geom_density(
      fill = "#66c2a5", alpha = 0.7,
      color = "black"
    ) +
    ggplot2::facet_wrap(~Property, scales = "free", ncol = 4) +
    ggplot2::labs(
      title = "Key Physicochemical Properties Distribution"
    ) +
    ggplot2::theme_bw(base_size = 11)

  # Top50 heatmap（使用 pheatmap）
  top50 <- phys |>
    dplyr::arrange(dplyr::desc(Duplicate)) |>
    dplyr::slice_head(n = 50)
  mat <- top50 |>
    dplyr::select(
      GRAVY, Charge_pH7, pI, Aliphatic,
      Aromatic_prop, Basic_prop, Acidic_prop
    ) |>
    scale(center = TRUE, scale = TRUE) |>
    as.matrix()
  rownames(mat) <- paste0(
    "Clone", 1:50, " (n=",
    top50$Duplicate, ")"
  )

  pheatmap::pheatmap(
    mat,
    color = grDevices::colorRampPalette(
      c("#2c7bb6", "white", "#d7191c")
    )(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = "Top 50 Expanded Clones – Physicochemical Z-scores",
    device = ragg::agg_png
  )

  final_plot <- (p_atchley | p_pca) /
    (p_ai | p_dist) +
    patchwork::plot_annotation(
      title = "Comprehensive Physicochemical Landscape of TCR/BCR CDR3 Repertoire",
      subtitle = paste(
        sample_name,
        "| Unique sequences:", nrow(phys),
        "| Total reads:", scales::comma(sum(phys$Duplicate))
      ),
      theme = theme_shield(base_size = 18)
    )

  if (!is.null(output_prefix)) {
    ggplot2::ggsave(
      paste0(output_prefix, ".pdf"),
      final_plot, width = 18, height = 13,
      dpi = 300, device = grDevices::cairo_pdf
    )
    ggplot2::ggsave(
      paste0(output_prefix, ".png"),
      final_plot, width = 18, height = 13,
      dpi = 300
    )
  }

  invisible(final_plot)
}
