#' Summarise V/J usage for two groups of samples (list vs list)
#'
#' @param list_case    实验组样本列表，每个元素是一个 data.frame（一个样本）
#' @param list_control 对照组样本列表，每个元素是一个 data.frame（一个样本）
#' @param gene         "v_call" 或 "j_call"
#' @param method       汇总方式：\code{"sum"}, \code{"max"}, \code{"mean"}, \code{"ratio"}
#' @param case_label   实验组标签（默认 "Patient"）
#' @param control_label对照组标签（默认 "Control"）
#'
#' @return list，包含：
#'   \item{count_mat}{基因 × 样本 的矩阵（行为基因，列为样本）}
#'   \item{condition}{因子向量，长度等于列数，对应 Control / Patient}
#'
#' @export
shield_vj_summarise_lists <- function(list_case,
                                      list_control,
                                      gene   = c("v_call", "j_call"),
                                      method = c("sum", "max", "mean", "ratio"),
                                      case_label    = "Patient",
                                      control_label = "Control") {
  gene   <- match.arg(gene)
  method <- match.arg(method)

  if (is.null(names(list_case))) {
    names(list_case) <- paste0("case", seq_along(list_case))
  }
  if (is.null(names(list_control))) {
    names(list_control) <- paste0("ctrl", seq_along(list_control))
  }

  # 所有基因全集：两个列表一起
  all_genes <- unique(c(
    unlist(lapply(list_case,    function(df) df[[gene]])),
    unlist(lapply(list_control, function(df) df[[gene]]))
  ))
  all_genes <- all_genes[!is.na(all_genes)]

  # 内部函数：汇总一个列表 → 矩阵（基因 × 样本）
  summarise_one_list <- function(sample_list, prefix) {
    summarise_one_sample <- function(df) {
      if (!gene %in% names(df)) {
        stop("找不到列: ", gene, "\n当前列名为: ",
             paste(names(df), collapse = ", "))
      }
      if (!"duplicate_count" %in% names(df)) {
        stop("找不到列: duplicate_count")
      }

      dt <- data.table::as.data.table(df)

      if (method == "sum") {
        tab <- dt[
          ,
          list(value = sum(duplicate_count, na.rm = TRUE)),
          by = gene
        ]
      } else if (method == "max") {
        tab <- dt[
          ,
          list(value = max(duplicate_count, na.rm = TRUE)),
          by = gene
        ]
      } else if (method == "mean") {
        tab <- dt[
          ,
          list(value = mean(duplicate_count, na.rm = TRUE)),
          by = gene
        ]
      } else { # ratio
        tab <- dt[
          ,
          {
            vals <- sort(duplicate_count[!is.na(duplicate_count)])
            if (length(vals) == 0) {
              list(value = NA_real_)
            } else {
              n10 <- max(1, floor(length(vals) * 0.1))
              top_10    <- sum(tail(vals,   n10))
              bottom_10 <- sum(head(vals,   n10))
              list(value = ifelse(bottom_10 == 0, NA_real_, top_10 / bottom_10))
            }
          },
          by = gene
        ]
      }

      out <- rep(NA_real_, length(all_genes))
      names(out) <- all_genes
      idx <- match(tab[[gene]], all_genes)
      out[idx] <- tab[["value"]]
      out
    }

    mat <- sapply(sample_list, summarise_one_sample)
    colnames(mat) <- paste0(prefix, "_", names(sample_list))
    rownames(mat) <- all_genes
    mat
  }

  mat_case    <- summarise_one_list(list_case,    prefix = case_label)
  mat_control <- summarise_one_list(list_control, prefix = control_label)

  # 合并两个矩阵（行名已经统一是 all_genes）
  count_mat <- cbind(mat_control, mat_case)

  condition <- factor(
    c(
      rep(control_label, ncol(mat_control)),
      rep(case_label,    ncol(mat_case))
    ),
    levels = c(control_label, case_label)
  )

  list(
    count_mat = count_mat,
    condition = condition
  )
}

#' Differential V/J usage using DESeq2 (two-group lists)
#'
#' @param list_case    实验组样本列表，每个元素为 data.frame
#' @param list_control 对照组样本列表，每个元素为 data.frame
#' @param gene         "v_call" 或 "j_call"
#' @param method       "sum"/"max"/"mean"/"ratio"
#' @param case_label   实验组标签
#' @param control_label对照组标签
#'
#' @return list，包含：
#'   \item{dds}{DESeqDataSet 对象}
#'   \item{res}{DESeq2 结果表}
#'   \item{volcano_plot}{火山图}
#'   \item{cohend}{每个基因的 Cohen's d}
#'   \item{cohend_plot}{Cohen's d 条形图}
#'   \item{normalized_counts}{归一化计数矩阵}
#'
#' @export
shield_vj_deseq_lists <- function(list_case,
                                  list_control,
                                  gene   = c("v_call", "j_call"),
                                  method = c("sum", "max", "mean", "ratio"),
                                  case_label    = "Patient",
                                  control_label = "Control") {

  gene   <- match.arg(gene)
  method <- match.arg(method)

  # 先做汇总，得到矩阵 + condition
  summed <- shield_vj_summarise_lists(
    list_case    = list_case,
    list_control = list_control,
    gene         = gene,
    method       = method,
    case_label    = case_label,
    control_label = control_label
  )

  count_mat <- summed$count_mat
  condition <- summed$condition

  # 替换 NA 为 0，四舍五入为整数
  count_mat[is.na(count_mat)] <- 0
  count_mat <- round(count_mat)

  col_data <- data.frame(condition = condition)
  rownames(col_data) <- colnames(count_mat)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = col_data,
    design    = ~ condition
  )
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)

  volcano_data <- data.frame(
    log2_fc = res$log2FoldChange,
    p_value = res$padj,
    gene    = rownames(res)
  )

  volcano_plot <- ggplot2::ggplot(
    volcano_data,
    ggplot2::aes(x = log2_fc, y = -log10(p_value))
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = p_value < 0.05),
      size = 3, alpha = 0.8
    ) +
    ggplot2::scale_color_manual(values = c("gray80", "#D62728")) +
    ggplot2::labs(
      title = "V/J Differential Usage (DESeq2)",
      x     = "Log2 Fold Change",
      y     = "-log10(p-adj)"
    ) +
    theme_shield(base_size = 14) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = ifelse(p_value < 0.05, gene, "")),
      size = 3
    ) +
    ggplot2::theme(legend.position = "none")

  # ---- 计算 Cohen's d ----
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  control_idx <- which(condition == control_label)
  case_idx    <- which(condition == case_label)

  cohend_values <- apply(norm_counts, 1, function(x) {
    c_vals <- x[control_idx]
    p_vals <- x[case_idx]
    m1 <- mean(c_vals); m2 <- mean(p_vals)
    s1 <- stats::sd(c_vals); s2 <- stats::sd(p_vals)
    n1 <- length(c_vals); n2 <- length(p_vals)
    if (n1 + n2 <= 2 || (s1 == 0 && s2 == 0)) return(NA_real_)
    pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    (m2 - m1) / pooled_sd
  })

  cohend_df <- data.frame(
    gene   = rownames(norm_counts),
    cohend = cohend_values
  )

  cohend_plot <- ggplot2::ggplot(
    cohend_df,
    ggplot2::aes(x = gene, y = cohend)
  ) +
    ggplot2::geom_col(
      fill = ifelse(cohend_df$cohend > 0, "#D62728", "gray80")
    ) +
    ggplot2::labs(
      title = "Cohen's d for V/J usage",
      x     = "Gene",
      y     = "Cohen's d"
    ) +
    theme_shield(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5, size = 8
      )
    )

  list(
    dds               = dds,
    res               = res,
    volcano_plot      = volcano_plot,
    cohend            = cohend_values,
    cohend_plot       = cohend_plot,
    normalized_counts = norm_counts,
    count_mat         = count_mat,
    condition         = condition
  )
}
