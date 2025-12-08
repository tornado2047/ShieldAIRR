#' ShieldAIRR ggplot2 theme
#'
#' 一个统一的绘图主题，整体风格与 CDR3 理化性质图保持一致：
#' theme_minimal 基础，绿色标题、灰色副标题，去掉多余网格。
#'
#' @param base_size 基础字号，默认 14。
#'
#' @return 一个 ggplot2 theme 对象。
#' @export
theme_shield <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = base_size + 2,
        face = "bold",
        color = "#1b9e77"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 1,
        color = "gray30"
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
}
