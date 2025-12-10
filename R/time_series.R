
#' Time-series: make long-format clonotype table
#'
#' @param dfs 命名列表，每个元素是一个 data.frame（一个时间点的 AIRR 数据）
#' @param clonotype_col 定义克隆的列名，默认 "junction_aa"
#' @param abundance_col 丰度列名，默认 "duplicate_count"
#' @param min_count 保留的最小丰度阈值
#'
#' @return data.frame，包含 junction_aa, duplicate_count, time, l_abundance
#' @export
make_long <- function(
    dfs,
    clonotype_col   = "junction_aa",
    abundance_col   = "duplicate_count",
    min_count       = 1
) {
  if (is.null(names(dfs)) || any(names(dfs) == "")) {
    names(dfs) <- seq_along(dfs) - 1L
  }

  out_list <- vector("list", length(dfs))
  i <- 0L

  for (nm in names(dfs)) {
    i <- i + 1L
    df <- dfs[[i]]

    if (!clonotype_col %in% names(df)) {
      stop("在时间点 ", nm, " 的数据中找不到列: ", clonotype_col,
           "\n当前列名为: ", paste(names(df), collapse = ", "))
    }
    if (!abundance_col %in% names(df)) {
      stop("在时间点 ", nm, " 的数据中找不到列: ", abundance_col,
           "\n当前列名为: ", paste(names(df), collapse = ", "))
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

#' Time-series: summarise clonotype features
#'
#' @param long make_long() 的输出
#'
#' @return data.frame，每个 clonotype 对应一行
#' @export
summarise_clonotypes <- function(long) {

  n_time_levels <- dplyr::n_distinct(long$time)

  slope_fun <- function(x_time, x_lab) {
    if (length(unique(x_time)) < 2) return(0)  # 少于 2 个时间点，斜率记为 0
    stats::coef(stats::lm(x_lab ~ x_time))[2]
  }

  long %>%
    dplyr::group_by(junction_aa) %>%
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
    ) %>%
    dplyr::mutate(
      persistence = n_time / n_time_levels
    )
}

#' Time-series: cluster clonotypes by trajectory features
#'
#' @param clono_features summarise_clonotypes() 的输出
#' @param k 聚类簇数
#' @param min_time 出现的最少时间点
#' @param min_tot 最小累计丰度
#'
#' @return data.frame，增加一列 cluster
#' @export
cluster_clonotypes <- function(clono_features,
                               k = 6,
                               min_time = 2,
                               min_tot = 0) {

  clono_for_cluster <- clono_features %>%
    dplyr::filter(
      n_time >= min_time,
      tot    >= min_tot
    )

  if (nrow(clono_for_cluster) < k) {
    stop("可用于聚类的克隆数量（",
         nrow(clono_for_cluster),
         "）少于 k = ", k,
         "，请降低 k 或放宽过滤条件。")
  }

  feature_mat <- clono_for_cluster %>%
    dplyr::select(l_ab_mean, slope, persistence, peak_time) %>%
    as.matrix()

  # 去掉 NA / Inf
  bad <- !is.finite(feature_mat)
  if (any(bad)) {
    keep_row <- apply(!bad, 1, all)
    feature_mat <- feature_mat[keep_row, , drop = FALSE]
    clono_for_cluster <- clono_for_cluster[keep_row, , drop = FALSE]
  }

  if (nrow(feature_mat) < k) {
    stop("去除 NA/Inf 后可用于聚类的克隆数量（",
         nrow(feature_mat),
         "）少于 k = ", k,
         "，请降低 k 或放宽过滤条件。")
  }

  feature_mat <- scale(feature_mat)

  set.seed(123)
  km <- stats::kmeans(feature_mat, centers = k, nstart = 20)

  clono_for_cluster %>%
    dplyr::mutate(cluster = paste0("C", km$cluster))
}

#' Time-series: plot cluster trajectories
#'
#' @param long make_long() 的输出
#' @param clono_clustered cluster_clonotypes() 的输出
#' @param k 聚类簇数
#'
#' @export
plot_cluster_traj <- function(long, clono_clustered, k) {

  long_plot <- long %>%
    dplyr::inner_join(
      clono_clustered %>% dplyr::select(junction_aa, cluster),
      by = "junction_aa"
    )

  mean_traj <- long_plot %>%
    dplyr::group_by(cluster, time) %>%
    dplyr::summarise(
      mean_l_ab = mean(l_abundance),
      .groups   = "drop"
    )

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = long_plot,
      ggplot2::aes(x = time, y = l_abundance,
                   group = junction_aa, colour = cluster),
      alpha = 0.1
    ) +
    ggplot2::geom_line(
      data = mean_traj,
      ggplot2::aes(x = time, y = mean_l_ab, colour = cluster),
      linewidth  = 1.2,
      alpha = 0.9
    ) +
    ggplot2::facet_wrap(~ cluster, scales = "free_y") +
    ggplot2::labs(
      x     = "time",
      y     = "log10(duplicate_count + 1)",
      title = paste("Clone abundance trajectories (k =", k, ")")
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none")
}
