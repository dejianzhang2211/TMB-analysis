## --------------------------------------------------------------
## prepare miRNA expr -- RPM
## --------------------------------------------------------------
rm(list=ls())
gc()


## ---------------------------------------------
## 用R包TCGAbiolinks下载LUAD患者的miRNA表达谱数据
## ---------------------------------------------
if(F){
  # 加载R包
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  if(!require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks",ask = F,update = F)
  library(TCGAbiolinks)
  library(dplyr)
  library(DT)
  
  
  ## step1:下载miRNA成熟体数据
  query <- GDCquery(project = c('TCGA-LUAD'), 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Isoform Expression Quantification")  # 包含miRNA成熟体，但是需要自己整理
  GDCdownload(query)
  expr_iso = GDCprepare(query)
  
  
  ## step2:整理下载的原始表达谱: 将长数据转化为宽数据
  library(dplyr)
  library(reshape2)
  
  expr_iso2 = as.data.frame(expr_iso)
  expr_iso3 = expr_iso2 %>% group_by(barcode, miRNA_region) %>%  # 将每个样品中的相同的成熟miRNA的counts数合并
    summarize(read_count = sum(read_count))
  expr_iso4 = as.data.frame(expr_iso3)
  expr_iso5 = dcast(expr_iso4, miRNA_region ~ barcode, value.var='read_count')  # 将长数据转化为宽数据，转换成表达谱
  expr_iso5 = expr_iso5[-c(2214,2215,2216),]  # 去掉不是miRNA-id的行
  expr_iso5$miRNA_region = substr(expr_iso5$miRNA_region,8,19)  # 整理miRNA_region列的字符
  
  
  ## step3:将accession number转换为miRNA_ID
  # 获取miR_accession和miR_id之间的对应
  if(!require('miRBaseVersions.db')) BiocManager::install('miRBaseVersions.db')
  miR_region = expr_iso5$miRNA_region
  items = select(miRBaseVersions.db, 
                 keys=miR_region,
                 keytype = 'MIMAT',
                 columns = '*')
  acce2name = items[items$VERSION == 21.0, c('ACCESSION','NAME')] # 只取miRBase v21的结果
  # 将表达谱中的accession转换为id
  expr_iso5 = expr_iso5[order(expr_iso5$miRNA_region),]
  acce2name = acce2name[order(acce2name$ACCESSION),]
  rownames(expr_iso5) = acce2name[match(expr_iso5$miRNA_region, acce2name$ACCESSION),'NAME']
  
  expr_miR_luad = expr_iso5[,-1]
  save(expr_miR_luad, file='Rdata/TCGA_LUAG_miRNA_expr_TCGAbiolinks.Rdata')
  
}
