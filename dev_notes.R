
getwd()


usethis::use_git()      # 初始化 git
usethis::use_github()   # 直接在 GitHub 建仓并 push


usethis::use_mit_license("XF Zhang")
usethis::use_r("theme_shield.R")
usethis::use_r("repertoire_summary")
usethis::use_r("cdr3_physchem")
usethis::use_r("time_series")
usethis::use_r("bulk_vj_deseq")
usethis::use_r("dev_notes")

## 未执行
devtools::document()
devtools::check()
devtools::install()


usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)



install.packages(c("alakazam", "ape", "CollessLike", "data.table", "dplyr", "entropy", "jsonlite", "magrittr", "Peptides", "RecordLinkage", "shazam", "seqinr", "stringdist", "stringr", "testthat", "textmineR", "yaml"))
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library(devtools)

devtools::load_all("/Users/xfcheung/workspace/AIDeN/sumrep/")
library(sumrep)
