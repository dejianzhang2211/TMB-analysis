## ----------------------------------------------------------
## 获取luad和lusc的表达谱及其对应得clinical信息
## ----------------------------------------------------------

# 清空当前工作环境
rm(list=ls())
gc()
options(stringsAsFactors = F)



## -----------------------------------------
## 用R包TCGAbiolinks下载表达谱数据：FPKM
## -----------------------------------------
if(F){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  if(!require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks",ask = F,update = F)
  library(TCGAbiolinks)
  library(dplyr)
  library(DT)
  
  query <- GDCquery(project = "TCGA-LUAD", 
                    legacy = F, 
                    experimental.strategy = "RNA-Seq", 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - FPKM")   # HTSeq - Counts(原始count数)  HTSeq - FPKM(FPKM值/表达量值)
  GDCdownload(query)
  expr = GDCprepare(query)
  expr2 = assay(expr)  # 获取表达矩阵
  save(expr2, file='Rdata/TCGA_LUAD_RNAseq_expr_FPKM.Rdata')
}


## -----------------------------------
## 获取RNA-seq数据中所有患者的id
## -----------------------------------
if(F){
  ## load expr
  load('Rdata/TCGA_LUAD_RNAseq_expr_FPKM.Rdata')
  expr = expr_RNA_FPKM
  
  group_list = ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
  table(group_list)
  expr = as.data.frame(expr)
  expr_tumor  = expr[,group_list=='tumor']
  dim(expr_tumor)
  LUAD_pid_RNA = substr(colnames(expr_tumor),1,12)
  
  save(LUAD_pid_RNA, file='Rdata/TCGA_LUAD_pid_RNA.Rdata')
  
}


