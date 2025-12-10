
## ----time-series-example, eval=FALSE-----------------------------------
library(ShieldAIRR)
library(dplyr)
# 假设前面已经构造好了 RA_Patient 列表
list1 <- RA_Patient
set.seed(2025)
dfs <- list1[sample(seq_along(list1), 10)]
names(dfs) <- 0:9
long <- make_long(dfs,
  clonotype_col = "junction_aa",
  abundance_col = "duplicate_count",
  min_count     = 1
)
feat <- summarise_clonotypes(long)
k   <- 6
clu <- cluster_clonotypes(
  clono_features = feat,
  k              = k,
  min_time       = 3,
  min_tot        = 100
)
p_traj <- plot_cluster_traj(long, clu, k = k)
p_traj
张雪飞
  下午 3:25
# ShieldAIRR
ShieldAIRR 是一个用于 **TCR/BCR AIRR-seq 数据分析** 的 R 包。
它提供从 **V/J 基因使用、组间差异、单样本概览、时间序列追踪** 的完整分析工作流。
---
## :扳手: 安装
### 1. 从 GitHub 安装包
```r
# 如有需要先安装 remotes：
# install.packages("remotes")
remotes::install_github("tornado2047/ShieldAIRR")
library(ShieldAIRR)
张雪飞
  下午 3:31
# ShieldAIRR
ShieldAIRR 是一个用于 TCR/BCR AIRR-seq 数据分析的 R 包，提供从单样本分析、组间差异分析，到时间序列克隆轨迹建模的全流程工具。
本包的设计理念是模块化、可扩展、与 sumrep、Immcantation 等生态兼容，方便进行高质量可视化与统计分析。
====================================================================
安装 Installation
====================================================================
1. 从 GitHub 安装 ShieldAIRR：
install.packages("remotes")
remotes::install_github("tornado2047/ShieldAIRR")
library(ShieldAIRR)
2. 推荐：一键安装全部依赖（包括 sumrep 和 CollessLike 的本地安装）
shield_install_deps(
    sumrep_path      = "/path/to/sumrep",
    colless_tar_path = "/path/to/CollessLike_2.0.tar.gz"
)
该函数会自动安装：
- CRAN 依赖
- Bioconductor 依赖
- 本地 sumrep
- 本地 CollessLike
====================================================================
输入数据格式 Input Format
====================================================================
ShieldAIRR 支持所有 AIRR 风格数据框（data.frame），至少需包含以下列：
junction_aa        CDR3 氨基酸序列
duplicate_count     克隆丰度
v_call              V 基因注释
j_call              J 基因注释
多个样本应组织成如下结构：
RA_Control <- list(sampleA_df, sampleB_df, ...)
RA_Patient <- list(sampleC_df, sampleD_df, ...)
====================================================================
示例 1：单样本 Repertoire 综合分析图
====================================================================
df_demo <- RA_Patient[[1]]
summarizeRepertoirePlot(
    df_demo,
    sample_name = "Demo sample",
    output_pdf  = FALSE
)
输出包含：
- CDR3 长度分布
- V / J 基因使用
- V-J 配对矩阵
- Rank-Abundance 克隆丰度曲线
- 统一主题 theme_shield()
====================================================================
示例 2：对照组 vs 实验组 的 V/J 基因差异分析
====================================================================
res_v <- shield_vj_deseq_lists(
    list_control = RA_Control,
    list_case    = RA_Patient,
    gene         = "v_call",
    method       = "mean"
)
res_j <- shield_vj_deseq_lists(
    list_control = RA_Control,
    list_case    = RA_Patient,
    gene         = "j_call",
    method       = "mean"
)
查看火山图和 Cohen's d 图：
res_v$volcano_plot
res_j$volcano_plot
res_v$cohend_plot
res_j$cohend_plot
查看差异结果：
head(res_v$res)
====================================================================
示例 3：时间序列克隆轨迹建模（10 个时间点模拟）
====================================================================
# 从 RA_Patient 随机选择 10 个样本模拟 10 个时间点
set.seed(2025)
dfs <- RA_Patient[sample(seq_along(RA_Patient), 10)]
names(dfs) <- 0:9
# Step 1: 转换为 long-format 数据
long <- make_long(dfs)
# Step 2: 提取克隆的时间特征
feat <- summarise_clonotypes(long)
# Step 3: 聚类轨迹（例如 k = 6）
clu <- cluster_clonotypes(
    clono_features = feat,
    k        = 6,
    min_time = 3,
    min_tot  = 100
)
# Step 4: 绘制轨迹聚类图
plot_cluster_traj(long, clu, k = 6)
====================================================================
典型完整工作流（推荐使用）
====================================================================
# 1. 单样本总结：
summarizeRepertoirePlot(RA_Patient[[1]], "Patient_1")
# 2. 差异分析：
res_v <- shield_vj_deseq_lists(RA_Control, RA_Patient, gene="v_call")
res_j <- shield_vj_deseq_lists(RA_Control, RA_Patient, gene="j_call")
# 3. 时间序列分析：
dfs  <- RA_Patient[sample(seq_along(RA_Patient), 10)]
names(dfs) <- 0:9
long <- make_long(dfs)
feat <- summarise_clonotypes(long)
clu  <- cluster_clonotypes(feat, k = 6)
plot_cluster_traj(long, clu, k = 6)
====================================================================
函数索引 Function Index
====================================================================
CDR3 理化性质：
  shield_cdr3_landscape()
单样本分析：
  summarizeRepertoirePlot()
  getSummaryStats()
V/J 基因使用分析：
  getVGeneDistributions()
  getJGeneDistributions()
差异分析（DESeq2）：
  shield_vj_summarise_lists()
  shield_vj_deseq_lists()
时间序列分析：
  make_long()
  summarise_clonotypes()
  cluster_clonotypes()
  plot_cluster_traj()
可视化主题：
  theme_shield()
依赖安装：
  shield_install_deps()
====================================================================
贡献指南 Contributing
====================================================================
欢迎提交：
- Bug report
- Feature request
- Pull request
- 新功能建议
====================================================================
版权 License
====================================================================
ShieldAIRR is released under the MIT License.










