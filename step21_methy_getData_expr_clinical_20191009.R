
## -----------------------------------------------------
## 获取TCGA的 LUAD 和 LUSC的 癌组织 methylation数据
## -----------------------------------------------------
## 推荐两种方式：TCGAbiolinks包，firehose网站
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
if(!require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks",ask = F,update = F)
if(!require("DT")) install.packages("DT",update = F,ask = F)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)


#下载甲基化数据
query_met <- GDCquery(project ='TCGA-LUAD',
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query_met)
expdat <- GDCprepare(query = query_met)
expr_methy_luad = assay(expdat)  # 将数据转换成matrix，以便后续使用
# save(expr_methy_luad, file='DNA_methylation_450k_LUAD.Rdata')
save(expr_methy_luad, file='TCGA_LUAD_DNA_methylation_450k_expr.Rdata')


query_met <- GDCquery(project ='TCGA-LUSC',
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query_met)
expdat <- GDCprepare(query = query_met)
expr_methy_lusc = assay(expdat)
# save(expr_methy_lusc, file='DNA_methylation_450k_LUSC.Rdata')
save(expr_methy_lusc, file='TCGA_LUSC_DNA_methylation_450k_expr.Rdata')



## ------------------------------------------------
## 表达谱数据预处理
## ------------------------------------------------
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

# 载入下载好的原始表达谱，beta value
load('Rdata/TCGA_LUAD_DNA_methylation_450k_expr.Rdata')
dim(expr_luad)
expr_luad = na.omit(expr_luad)
## 去掉detection pvalue<0.01的探针
## 用TCGAbiolinks下载的methylation expr数据是不是已经进行了beta value>0.01的过滤

## 去掉性染色体上的探针
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
dim(ann450k)

keep <- !(rownames(expr_luad) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
expr_luad <- expr_luad[keep,]



## ------------------------------------------------
## 获取methy数据中所有患者的id
## ------------------------------------------------
if(F){
  load('Rdata/TCGA_LUAD_DNA_methylation_450k_expr_clean.Rdata')
  expr = expr_luad
  
  group_list = ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
  table(group_list)
  expr = as.data.frame(expr)
  expr_tumor  = expr[,group_list=='tumor']
  dim(expr_tumor)
  LUAD_pid_cpg = substr(colnames(expr_tumor),1,12)
  
  save(LUAD_pid_cpg, file='Rdata/TCGA_LUAD_pid_cpg.Rdata')
}



## ---------------------------------------------
## 数据可视化
## ---------------------------------------------
if(F){
  library(ggplot2)
  
  ### 绘制单个样品的甲基化程度
  expr_luad = as.data.frame(expr_luad)
  # 样品：TCGA-55-8621-01A-11D-2398-05
  p = ggplot(expr_luad, aes(x=expr_luad[,1]))+   # expr_luad需要是data.frame格式
    geom_histogram(binwidth=0.02)+
    labs(x="% methylation per base",title = "Histogram of % CpG methylation\nSampleID:TCGA-55-8621-01A-11D-2398-05")+
    theme(text=element_text(size=20), axis.text.x=element_text(angle=30, hjust=1))
  ggsave(p, file='plots/methy_histogram_of_CpG_methylation4.png',width=15, height=10, dpi=300)
  
  
  ### 绘制多个样品甲基化程度小提琴图
  tmp <- expr_luad[,c(1,2,3,4,17,20,48,57)]   # 前4个tumor，后4个normal
  tmp <- data.frame(tmp)
  tmp2 <- stack(tmp)
  tmp2 <- data.frame(tmp2)
  colnames(tmp2) <- c("CpG methylation","SampleID")
  p = ggplot(tmp2,aes(x=SampleID,y=`CpG methylation`,fill=SampleID))+
    geom_violin()+
    theme(text=element_text(size=20), axis.text.x=element_text(angle=30, hjust=1))
  ggsave(p, file='plots/methy_violin_multi_sample2.png',width=15, height=10, dpi=300)
}

