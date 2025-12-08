# 内部工具：单样本按 V/J 汇总
summarise_vj_single_ <- function(df,
                                 gene_col = c("v_call", "j_call"),
                                 method = c("sum", "max", "mean", "ratio"),
                                 count_col = "duplicate_count") {
  gene_col <- rlang::arg_match(gene_col)
  method   <- rlang::arg_match(method)

  if (!gene_col %in% names(df)) {
    stop("数据中找不到列: ", gene_col)
  }
  if (!count_col %in% names(df)) {
    stop("数据中找不到列: ", count_col)
  }

  df <- dplyr::filter(df, !is.na(.data[[gene_col]]), !is.na(.data[[count_col]]))

  dplyr::group_by(df, .data[[gene_col]]) |>
    dplyr::summarise(
      stat = dplyr::case_when(
        method == "sum"  ~ sum(.data[[count_col]], na.rm = TRUE),
        method == "max"  ~ max(.data[[count_col]], na.rm = TRUE),
        method == "mean" ~ mean(.data[[count_col]], na.rm = TRUE),
        method == "ratio" ~ {
          vals <- sort(.data[[count_col]][!is.na(.data[[count_col]])])
          if (length(vals) == 0) {
            NA_real_
          } else {
            k <- max(1L, floor(length(vals) * 0.1))
            top_10    <- tail(vals, k)
            bottom_10 <- head(vals, k)
            if (sum(bottom_10) == 0) NA_real_ else sum(top_10) / sum(bottom_10)
          }
        },
        TRUE ~ NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::rename(gene = .data[[gene_col]])
}

#' Summarise V/J gene usage across multiple samples
#'
#' 给定一个样本列表（每个样本是 AIRR 风格的数据框），对 V/J 基因做
#' sum / max / mean / top:bottom ratio 等汇总，并返回：
#' \itemize{
#'   \item \code{summary_list}: 每个样本一个 data.frame，含 \code{gene} 和 \code{stat}
#'   \item \code{count_matrix}: 行为基因，列为样本的矩阵，可直接用于 DESeq2
#' }
#'
#' @param TRA_list 命名的 data.frame 列表，每个元素是一份样本的 TCR/BCR 表。
#' @param gene 选择 "v_call" 或 "j_call"。
#' @param method "sum", "max", "mean", "ratio" 之一。
#' @param count_col 表示丰度的列，默认 "duplicate_count"。
#'
#' @return 一个列表：\code{list(summary_list = ..., count_matrix = ..., gene = gene, method = method)}
#' @export
shield_vj_summarise <- function(TRA_list,
                                gene = c("v_call", "j_call"),
                                method = c("sum", "max", "mean", "ratio"),
                                count_col = "duplicate_count") {
  stopifnot(is.list(TRA_list), length(TRA_list) > 0)
  if (is.null(names(TRA_list))) {
    stop("TRA_list 必须是命名列表（每个样本一个名字）")
  }

  gene   <- rlang::arg_match(gene)
  method <- rlang::arg_match(method)

  # 每个样本单独汇总
  summary_list <- lapply(
    TRA_list,
    summarise_vj_single_,
    gene_col  = gene,
    method    = method,
    count_col = count_col
  )

  # 统一所有基因名，填 0
  all_genes <- sort(unique(unlist(lapply(summary_list, \(x) x$gene))))
  mat <- matrix(
    0,
    nrow = length(all_genes),
    ncol = length(summary_list),
    dimnames = list(all_genes, names(summary_list))
  )

  for (i in seq_along(summary_list)) {
    df_i <- summary_list[[i]]
    idx  <- match(df_i$gene, all_genes)
    mat[idx, i] <- df_i$stat
  }

  list(
    summary_list = summary_list,
    count_matrix = mat,
    gene         = gene,
    method       = method
  )
}

# 内部工具：从列名猜 group（Control/Patient）
infer_group_from_names_ <- function(sample_names) {
  grp <- rep(NA_character_, length(sample_names))

  grp[grepl("Control", sample_names, ignore.case = TRUE)] <- "Control"
  grp[grepl("Patient", sample_names, ignore.case = TRUE)] <- "Patient"

  if (any(is.na(grp))) {
    stop(
      "无法从样本名推断 group，请显式提供 group 向量。\n",
      "样本名: ", paste(sample_names, collapse = ", ")
    )
  }
  factor(grp, levels = c("Control", "Patient"))
}

#' Differential V/J gene usage using DESeq2
#'
#' 对 V/J 基因计数矩阵进行 DESeq2 差异分析，并输出：
#' DESeq2 结果、火山图、Cohen's d 图和归一化计数。
#'
#' @param count_matrix 矩阵，行为 V/J 基因，列为样本。
#' @param group 样本分组因子（长度 = 样本数，两个水平，例如 Control/Patient）。
#'   如果为 NULL，则尝试从列名中通过正则匹配 "Control"/"Patient" 推断。
#' @param alpha FDR 阈值，默认 0.05。
#' @param label_top_n 在火山图上标注的 top 基因数，按 |log2FC| 排序，默认 15。
#'
#' @return 列表：\code{list(dds, res, volcano, significant_genes, cohen_d, cohen_plot, normalized_counts)}
#' @export
shield_vj_deseq <- function(count_matrix,
                            group = NULL,
                            alpha = 0.05,
                            label_top_n = 15) {
  stopifnot(is.matrix(count_matrix))

  # group
  if (is.null(group)) {
    group <- infer_group_from_names_(colnames(count_matrix))
  } else {
    if (length(group) != ncol(count_matrix)) {
      stop("group 长度必须等于样本数（列数）")
    }
    group <- as.factor(group)
    if (length(levels(group)) != 2) {
      stop("目前仅支持两组比较（group 必须有且仅有两个水平）")
    }
  }

  # 构建 DESeq2 对象
  counts_round <- round(count_matrix)
  counts_round[is.na(counts_round)] <- 0

  col_data <- data.frame(condition = group)
  rownames(col_data) <- colnames(counts_round)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_round,
    colData   = col_data,
    design    = ~ condition
  )
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  significant <- res_df[!is.na(res_df$padj) & res_df$padj < alpha, , drop = FALSE]

  # 火山图数据
  volcano_df <- res_df |>
    dplyr::mutate(
      neg_log10_padj = -log10(padj),
      sig = dplyr::case_when(
        is.na(padj)           ~ "ns",
        padj >= alpha         ~ "ns",
        log2FoldChange > 0    ~ "up",
        log2FoldChange < 0    ~ "down",
        TRUE                  ~ "ns"
      )
    )

  to_label <- volcano_df |>
    dplyr::filter(sig != "ns") |>
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) |>
    dplyr::slice_head(n = label_top_n)

  volcano_plot <- ggplot2::ggplot(
    volcano_df,
    ggplot2::aes(x = log2FoldChange, y = neg_log10_padj)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = sig),
      alpha = 0.8,
      size = 2.5
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "ns"   = "grey80",
        "up"   = "#d73027",
        "down" = "#4575b4"
      ),
      breaks = c("up", "down", "ns"),
      labels = c("Up in group2", "Down in group2", "NS"),
      name   = NULL
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) +
    ggrepel::geom_text_repel(
      data = to_label,
      ggplot2::aes(label = gene),
      size = 3,
      max.overlaps = 50,
      box.padding = 0.3,
      point.padding = 0.2
    ) +
    ggplot2::labs(
      title = "Differential V/J gene usage (DESeq2)",
      x = "log2 fold change (group2 vs group1)",
      y = "-log10 adjusted p-value"
    ) +
    theme_shield(base_size = 14)

  # Cohen's d (归一化计数)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  lvl <- levels(group)
  g1 <- which(group == lvl[1])
  g2 <- which(group == lvl[2])

  cohen_d <- apply(norm_counts, 1, function(x) {
    x1 <- x[g1]
    x2 <- x[g2]
    if (all(x1 == 0) && all(x2 == 0)) return(NA_real_)
    m1 <- mean(x1)
    m2 <- mean(x2)
    s1 <- stats::sd(x1)
    s2 <- stats::sd(x2)
    n1 <- length(x1)
    n2 <- length(x2)
    pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    if (isTRUE(all.equal(pooled_sd, 0))) return(NA_real_)
    (m2 - m1) / pooled_sd
  })

  cohen_df <- data.frame(
    gene    = rownames(norm_counts),
    cohen_d = cohen_d,
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(abs(cohen_d)))

  cohen_plot <- ggplot2::ggplot(
    cohen_df,
    ggplot2::aes(x = stats::reorder(gene, cohen_d), y = cohen_d)
  ) +
    ggplot2::geom_col(
      fill = "#1b9e77",
      alpha = 0.9
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Effect size (Cohen's d) per V/J gene",
      x = "Gene",
      y = "Cohen's d (group2 - group1)"
    ) +
    theme_shield(base_size = 13)

  list(
    dds               = dds,
    res               = res_df,
    volcano           = volcano_plot,
    significant_genes = significant,
    cohen_d           = cohen_df,
    cohen_plot        = cohen_plot,
    normalized_counts = norm_counts
  )
}
