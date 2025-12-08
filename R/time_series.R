#' Build long-form clonotype × time table
#'
#' @param dfs 命名的 data.frame 列表，每个元素为一个时间点的 AIRR 表格。
#'   列表名会被当作 time 变量（数值）使用。
#' @param clonotype_col 克隆定义列，默认 "junction_aa"。
#' @param abundance_col 丰度列，默认 "duplicate_count"。
#' @param min_count 过滤掉丰度小于该值的条目。
#'
#' @return data.frame，含列 junction_aa, duplicate_count, time, l_abundance。
#' @export
make_long <- function(
    dfs,
    clonotype_col = "junction_aa",
    abundance_col = "duplicate_count",
    min_count = 1) {

  out_list <- vector("list", length(dfs))
  i <- 0

  for (nm in names(dfs)) {
    i <- i + 1
    df <- dfs[[i]]

    if (!clonotype_col %in% names(df)) {
      stop(
        "在时间点 ", nm, " 的数据中找不到列: ", clonotype_col,
        "\n当前列名为: ", paste(names(df), collapse = ", ")
      )
    }
    if (!abundance_col %in% names(df)) {
      stop(
        "在时间点 ", nm, " 的数据中找不到列: ", abundance_col,
        "\n当前列名为: ", paste(names(df), collapse = ", ")
      )
    }

    clono <- df[[clonotype_col]]
    abund <- df[[abundance_col]]

    keep <- !is.na(abund) & abund >= min_count

    tmp <- data.frame(
      junction_aa     = clono[keep],
      duplicate_count = abund[keep],
      time            = as.numeric(nm),
      l_abundance     = log10(abund[keep] + 1),
      stringsAsFactors = FALSE
    )

    out_list[[i]] <- tmp
  }

  long <- do.call(rbind, out_list)
  rownames(long) <- NULL
  long
}

#' Summarise clonotype time-series features
#'
#' @param long 由 make_long() 生成的数据框。
#'
#' @return data.frame，每个 clonotype 一行，含 slope / persistence 等特征。
#' @export
summarise_clonotypes <- function(long) {
  n_time_levels <- dplyr::n_distinct(long$time)

  slope_fun <- function(x_time, x_lab) {
    if (length(unique(x_time)) < 2) return(NA_real_)
    stats::coef(stats::lm(x_lab ~ x_time))[2]
  }

  long |>
    dplyr::group_by(junction_aa) |>
    dplyr::summarise(
      n_time     = dplyr::n_distinct(time),
      first_time = min(time),
      last_time  = max(time),
      tot        = sum(duplicate_count),
      mean_ab    = mean(duplicate_count),
      l_ab_mean  = mean(l_abundance),
      slope      = slope_fun(time, l_abundance),
      peak       = max(duplicate_count),
      peak_time  = time[which.max(duplicate_count)],
      .groups    = "drop"
    ) |>
    dplyr::mutate(
      persistence = n_time / n_time_levels
    )
}

#' k-means clustering on clonotype features
#'
#' @param clono_features summarise_clonotypes() 的输出。
#' @param k 聚类数。
#'
#' @return 增加 cluster 列的 data.frame。
#' @export
cluster_clonotypes <- function(clono_features, k = 6) {
  feature_mat <- clono_features |>
    dplyr::select(l_ab_mean, slope, persistence, peak_time) |>
    scale() |>
    as.matrix()

  set.seed(123)
  km <- stats::kmeans(feature_mat, centers = k, nstart = 20)

  clono_features |>
    dplyr::mutate(cluster = paste0("C", km$cluster))
}

#' Plot clonotype abundance trajectories by clusters
#'
#' @param long make_long() 的输出。
#' @param clono_clustered cluster_clonotypes() 的输出。
#' @param k 聚类数，仅用于标题。
#'
#' @export
plot_cluster_traj <- function(long, clono_clustered, k) {
  long_plot <- long |>
    dplyr::inner_join(
      clono_clustered |>
        dplyr::select(junction_aa, cluster),
      by = "junction_aa"
    )

  mean_traj <- long_plot |>
    dplyr::group_by(cluster, time) |>
    dplyr::summarise(
      mean_l_ab = mean(l_abundance),
      .groups   = "drop"
    )

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = long_plot,
      ggplot2::aes(
        x = time, y = l_abundance,
        group = junction_aa, colour = cluster
      ),
      alpha = 0.1
    ) +
    ggplot2::geom_line(
      data = mean_traj,
      ggplot2::aes(
        x = time, y = mean_l_ab,
        colour = cluster
      ),
      size  = 1.2,
      alpha = 0.9
    ) +
    ggplot2::facet_wrap(~cluster, scales = "free_y") +
    ggplot2::labs(
      x     = "time",
      y     = "log10(duplicate_count + 1)",
      title = paste("Clone abundance trajectories (k =", k, ")")
    ) +
    theme_shield(base_size = 14) +
    ggplot2::theme(legend.position = "none")
}
