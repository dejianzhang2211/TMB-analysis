## -------------------------------------------
## cpg位点mapping到基因
## -------------------------------------------
rm(list=ls())
gc()

library(GenomicAlignments)
library(tidyr)


### 获取mapping到基因组promoter区域的cpg位点
if(F){
  hm450_hg19 = readRDS('Rdata/hm450.hg19.manifest.rds')   # 载入处理过的450k芯片注释文件
  hm450_hg19_2 = granges(keepSeqlevels(hm450_hg19, paste0("chr", 1:22), pruning.mode="coarse"))   # 只保留1-22号染色上的cpg位点
  region_promoter = readRDS('Rdata/region_hs_promoter_TSS_2000_1000.rds')   # 载入promoter区域的granges文件
  cpg_promoter = hm450_hg19_2[queryHits(findOverlaps(hm450_hg19_2,region_promoter))]
  cpg_promoter$names = names(cpg_promoter)
  cpg_promoter = unique(cpg_promoter)
  save(cpg_promoter, file='Rdata/cpg_in_promoter_TSS_2000_1000.Rdata')
}
load('Rdata/cpg_in_promoter_TSS_2000_1000.Rdata')
cpg_promoter_df = as.data.frame(cpg_promoter)  # cpg_promoter预处理
colnames(cpg_promoter_df)[6] = 'Name'

# load DEG
load('Rdata/TCGA_LUAD_training_cohort_cpg_DEG_anno_TMBH_TMBL_20200603.Rdata')
DEG_cpg_TMBH_TMBL_ann = as.data.frame(DEG_cpg_TMBH_TMBL_ann)
DEG_cpg_TMBH_TMBL_ann = subset(DEG_cpg_TMBH_TMBL_ann, adj.P.Val < 0.05 & abs(logFC) > 0.15)

## 从差异甲基化cpg位点中提取出mapping到promoter区域的位点
# DEG_cpg_TMBH_TMBL_ann_promoter = merge(DEG_cpg_TMBH_TMBL_ann, cpg_promoter_df, by='Name', all=F)
# dim(DEG_cpg_TMBH_TMBL_ann_promoter)


## 获取cpg 在gene上的位置
if(F){
  cpg2geneSite = DEG_cpg_TMBH_TMBL_ann[, c(1,16)]
  # 去掉没有mapping到基因上的cpg位点
  cpg2geneSite[cpg2geneSite==''] = NA
  cpg2geneSite = na.omit(cpg2geneSite)
  # 将gene_symbol列的内容按照；分隔开，拆分成新的列
  cpg2geneSite2 = cpg2geneSite %>% separate_rows(UCSC_RefGene_Group, sep=';')
  # 去掉cpg和gene name完全重复的行
  cpg2geneSite3 = cpg2geneSite2[!duplicated(cpg2geneSite2),]
  ## 重命名并保存
  cpg2geneSite_TMBH_TMBL = cpg2geneSite3
  
}


## 获取cpg site位于哪些基因
if(F){
  cpg2gene = DEG_cpg_TMBH_TMBL_ann[, c(1,14)]
  # 去掉没有mapping到基因上的cpg位点
  cpg2gene[cpg2gene==''] = NA
  cpg2gene = na.omit(cpg2gene)
  # 将gene_symbol列的内容按照；分隔开，拆分成新的列
  cpg2gene2 = cpg2gene %>% separate_rows(UCSC_RefGene_Name, sep=';')
  # 去掉cpg和gene name完全重复的行
  cpg2gene3 = cpg2gene2[!duplicated(cpg2gene2),]
  ## 重命名并保存
  cpg2gene_TMBH_TMBL = cpg2gene3
  
}

## 获取cpg site位于哪些基因,什么位置
if(F){
  cpg2gene = DEG_cpg_TMBH_TMBL_ann[, c(1,14,16)]
  # 去掉没有mapping到基因上的cpg位点
  cpg2gene[cpg2gene==''] = NA
  cpg2gene = na.omit(cpg2gene)
  # 将gene_symbol列的内容按照；分隔开，拆分成新的列
  cpg2gene2 = cpg2gene %>% separate_rows(UCSC_RefGene_Name, sep=';') %>% separate_rows(UCSC_RefGene_Group, sep=';')
  # 去掉cpg和gene name完全重复的行
  cpg2gene3 = cpg2gene2[!duplicated(cpg2gene2),]
  ## 重命名并保存
  cpg2gene_TMBH_TMBL = cpg2gene3
  
  save(cpg2gene_TMBH_TMBL, file='Rdata/TCGA_LUAD_cohort1_methy_DEG_cpg_to_gene.Rdata')
  
  # barplot: cpg在gene上的区域分布
  table(cpg2gene_TMBH_TMBL$UCSC_RefGene_Group)
  
  df = data.frame(group=c("1stExon", "3'UTR", "5'UTR", "Body", "TSS1500", "TSS200"), number=c(3, 3, 4, 27, 4, 1))
  library(ggplot2)
  ggplot(df, aes(x=group,y=number,fill=group))+
    geom_bar(stat="identity",position="identity")+
    theme_bw() + theme(panel.grid =element_blank())+
    theme(text=element_text(size=15))+
    geom_text(aes(label = number, vjust = -0.8, hjust = 0.5))   ## 显示柱条上的数字
  
}



## -----------------------------------------
## 获取基因TSS区域和promoter区域  get TSS region
## -----------------------------------------
if(F){
  if(!require("Homo.sapiens")) BiocManager::install("Homo.sapiens",ask = F,update = F)
  
  # 获取所有基因
  all_genes = genes(Homo.sapiens)  
  # 获取基因转录起始位点TSS
  all_gene_TSS = resize(all_genes, 1)
  head(start(all_gene_TSS),n=10)  # 查看前10个基因的TSS
  # 获得TSS上下游 各100 bp的位置
  TSS_100 = promoters(all_gene_TSS,100,100)
  # 获得TSS上游2000bp和下游100bp的位置 -- 基因的启动子部分
  promoter_hs_2000_1000 = promoters(all_gene_TSS,2000,1000)
  # 获取TSS上游2000bp和下游100bp的的序列
  library(BSgenome.Hsapiens.UCSC.hg19)
  seq = BSgenome.Hsapiens.UCSC.hg19
  promoter_seq = getSeq(seq, promoter_hs_2000_1000)
  
  save(promoter_hs_2000_1000, file='Rdata/region_hs_promoter_TSS_2000_1000.Rdata')
  saveRDS(promoter_hs_2000_1000, file='Rdata/region_hs_promoter_TSS_2000_1000.rds')
}


### 获取mapping到基因组promoter区域的cpg位点
if(T){
  hm450_hg19 = readRDS('Rdata/hm450.hg19.manifest.rds')   # 载入处理过的450k芯片注释文件
  hm450_hg19_2 = granges(keepSeqlevels(hm450_hg19, paste0("chr", 1:22), pruning.mode="coarse"))   # 只保留1-22号染色上的cpg位点
  region_promoter = readRDS('Rdata/region_hs_promoter_TSS_2000_1000.rds')   # 载入promoter区域的granges文件
  cpg_promoter = hm450_hg19_2[queryHits(findOverlaps(hm450_hg19_2,region_promoter))]
  cpg_promoter$names = names(cpg_promoter)
  cpg_promoter = unique(cpg_promoter)
  save(cpg_promoter, file='Rdata/cpg_in_promoter_TSS_2000_1000.Rdata')
}

## 从差异甲基化cpg位点中提取出mapping到promoter区域的位点
if(T){
  cpg_promoter_df = as.data.frame(cpg_promoter)  # cpg_promoter预处理
  colnames(cpg_promoter_df)[6] = 'Name'
  diff_CpG_sign_ann_promoter = merge(diff_CpG_sign_ann, cpg_promoter_df, by='Name', all=F)
  dim(diff_CpG_sign_ann_promoter)
}


