# ShieldAIRR

**ShieldAIRR** 是一个面向 **TCR / BCR AIRR-seq 数据** 的 R 语言分析框架，  
用于系统性刻画 **免疫组库组成、基因使用差异以及克隆随时间的动态行为**。

该包强调 **统计可解释性、模块化设计与生态兼容性**，  
可与 **Immcantation / sumrep** 无缝衔接，  
适用于基础免疫学研究、纵向队列分析及临床免疫监测场景。

---

## ✨ Key Features

- 🧬 **单样本免疫组库概览**
  - CDR3 长度分布
  - V / J 基因使用
  - V–J 配对矩阵
  - Rank–abundance 克隆结构

- 📊 **组间差异分析**
  - 基于 DESeq2 的 V/J 基因差异检验
  - 同时报告效应量（Cohen's *d*）
  - 内置火山图与统计可视化

- ⏱ **时间序列克隆轨迹建模**
  - 克隆纵向特征提取
  - 克隆轨迹聚类
  - 动态扩增 / 收缩模式可视化

- 🔗 **生态兼容**
  - 原生支持 AIRR / Immcantation 数据结构
  - 与 `sumrep`、`alakazam`、`shazam` 生态兼容

---

## 📦 Installation

### 1. 从 GitHub 安装 ShieldAIRR

```r
install.packages("remotes")
remotes::install_github("tornado2047/ShieldAIRR")
library(ShieldAIRR)

2.（推荐）一键安装全部依赖
ShieldAIRR 依赖部分 非 CRAN 包（如 sumrep、CollessLike）。
推荐使用内置函数完成统一部署：
shield_install_deps(
    sumrep_path      = "/path/to/sumrep",
    colless_tar_path = "/path/to/CollessLike_2.0.tar.gz"
)
该函数将自动安装并配置：
CRAN 依赖
Bioconductor 依赖
本地 sumrep
本地 CollessLike
sumrep 将在首次调用相关功能时自动加载。

:inbox_tray: Input Data Format
ShieldAIRR 接受 AIRR 标准风格 的 data.frame，
至少包含以下字段：
ColumnDescriptionjunction_aaCDR3 氨基酸序列duplicate_count克隆丰度v_callV 基因注释j_callJ 基因注释
多样本数据组织方式
RA_Control <- list(sampleA_df, sampleB_df, ...)
RA_Patient <- list(sampleC_df, sampleD_df, ...)

:test_tube: Example 1：单样本免疫组库综合概览
df <- RA_Patient[[1]]

summarizeRepertoirePlot(
    df,
    sample_name = "Patient_1",
    output_pdf  = FALSE
)
该函数生成一个标准化的免疫组库摘要图，适合：
单样本 QC
Supplementary figure
不同实验批次快速对比

:chart_with_upwards_trend: Example 2：V / J 基因组间差异分析
res_v <- shield_vj_deseq_lists(
    list_control = RA_Control,
    list_case    = RA_Patient,
    gene         = "v_call"
)

res_j <- shield_vj_deseq_lists(
    list_control = RA_Control,
    list_case    = RA_Patient,
    gene         = "j_call"
)
可视化输出
res_v$volcano_plot
res_v$cohend_plot
差异结果表
head(res_v$res)

⏱ Example 3：时间序列克隆轨迹分析
# 构造 10 个时间点示例
set.seed(2025)
dfs <- RA_Patient[sample(seq_along(RA_Patient), 10)]
names(dfs) <- 0:9

# 转换为 long format
long <- make_long(dfs)

# 提取克隆时间特征
feat <- summarise_clonotypes(long)

# 轨迹聚类
clu <- cluster_clonotypes(
    clono_features = feat,
    k        = 6,
    min_time = 3,
    min_tot  = 100
)

# 可视化克隆轨迹
plot_cluster_traj(long, clu, k = 6)
该模块适用于：
疫苗接种纵向追踪
治疗前后免疫反应评估
克隆扩增动力学建模

:repeat: Recommended Workflow
# 单样本 QC
summarizeRepertoirePlot(RA_Patient[[1]])

# 差异分析
res_v <- shield_vj_deseq_lists(RA_Control, RA_Patient, gene = "v_call")

# 时间序列分析
long <- make_long(RA_Patient[1:10])
feat <- summarise_clonotypes(long)
clu  <- cluster_clonotypes(feat, k = 6)
plot_cluster_traj(long, clu)

:books: Function Overview
单样本分析
summarizeRepertoirePlot()
getSummaryStats()
CDR3 特征
shield_cdr3_landscape()
V/J 基因使用
getVGeneDistributions()
getJGeneDistributions()
差异分析
shield_vj_summarise_lists()
shield_vj_deseq_lists()
时间序列分析
make_long()
summarise_clonotypes()
cluster_clonotypes()
plot_cluster_traj()
依赖管理
shield_install_deps()

:handshake: Contributing
欢迎以下形式的贡献：
:ladybug: Bug reports
:sparkles: Feature requests
:twisted_rightwards_arrows: Pull requests
:brain: 新分析指标或方法建议

:page_facing_up: License
ShieldAIRR is released under the MIT License.

---

如果你接下来想做的是：

- 📘 **vignette / 教程文档**
- 🐳 **Docker / Conda 部署说明**
- 🧪 **论文 Methods 风格算法说明版 README**

直接告诉我即可，我可以在这个 README 基础上继续无缝扩展。