## -------------------------------------------
## 甲基化位点注释
## -------------------------------------------
rm(list=ls())
gc()


library(dplyr)
library(tidyr)

# 准备450k芯片注释文件
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k_self = ann450k[c(1:4,9,18,19,24,25,26,28,29)]   # 挑选需要的内容


# 为差异甲基化分析结果添加注释
load('Rdata/TCGA_LUAD_training_cohort_cpg_DEG_TMBH_TMBL_20200603.Rdata')
DEG_cpg_TMBH_TMBL$Name = rownames(DEG_cpg_TMBH_TMBL)
DEG_cpg_TMBH_TMBL_ann = merge(DEG_cpg_TMBH_TMBL, ann450k_self, by='Name', all.x=T)
DEG_cpg_TMBH_TMBL_ann = DEG_cpg_TMBH_TMBL_ann[order(DEG_cpg_TMBH_TMBL_ann$adj.P.Val),]
save(DEG_cpg_TMBH_TMBL_ann, file='Rdata/TCGA_LUAD_training_cohort_cpg_DEG_anno_TMBH_TMBL_20200603.Rdata')


## 获取显著差异表达的cpg sites对应的gene
load('Rdata/TCGA_LUAD_training_cohort_cpg_DEG_anno_TMBH_TMBL_20200603.Rdata')
DEG_cpg_TMBH_TMBL_ann = as.data.frame(DEG_cpg_TMBH_TMBL_ann)
DEG_cpg_TMBH_TMBL_ann = subset(DEG_cpg_TMBH_TMBL_ann, adj.P.Val < 0.05 & abs(logFC) > 0.15)
## 获取初始的cpg to gene
cpg2gene = DEG_cpg_TMBH_TMBL_ann[,c(1,14)]


# 去掉没有mapping到基因上的cpg位点
cpg2gene[cpg2gene==''] = NA
cpg2gene = na.omit(cpg2gene)
# 将gene_symbol列的内容按照；分隔开，拆分成新的列
cpg2gene2 = cpg2gene %>% separate_rows(UCSC_RefGene_Name, sep=';')
# 去掉cpg和gene name完全重复的行
cpg2gene3 = cpg2gene2[!duplicated(cpg2gene2),]

## 重命名并保存
cpg2gene_TMBH_TMBL = cpg2gene3


save(cpg2gene_TMBH_TMBL, file='Rdata/TCGA_LUAD_training_cohort_cpg_DEG_cpg2gene_TMBH_TMBL_20200603.Rdata')





## -----------------------------------------------------------
## 可视化
## -----------------------------------------------------------
if(F){
  ## 火山图：padj<0.05, delta> \0.2\
  logFC_cutoff = 0.13
  diff_CpG_ann$change = as.factor(ifelse(diff_CpG_ann$adj.P.Val < 0.05 & abs(diff_CpG_ann$logFC) > logFC_cutoff,
                                         ifelse(diff_CpG_ann$logFC > logFC_cutoff ,'up','down'),'not'))
  this_tile <- paste0('Cutoff for logFC is ',logFC_cutoff,
                      '\nThe number of up gene is ',nrow(diff_CpG_ann[diff_CpG_ann$change =='up',]) ,
                      '\nThe number of down gene is ',nrow(diff_CpG_ann[diff_CpG_ann$change =='down',]))
  ## 绘图
  library(ggplot2)
  png('plot/methy_volcano_TMB_H_L.png',width=950,height=750,res=100)
  ggplot(data=as.data.frame(diff_CpG_ann), 
         aes(x=logFC, y=-log10(adj.P.Val), color=change)) +
    # geom_point(alpha=0.4, size=1.75) +
    geom_point() +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle(this_tile) + theme(plot.title = element_text(size=15,hjust = 0.5), panel.grid = element_blank())+
    scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
  dev.off()
  
  
  ## 箱线图：基因组不同组分上，差异甲基化cpg位点在normal组和tumor组样品中分布情况
  library(reshape2)
  
  load('Rdata/TCGA_LUAD_methy_expr_group_mean_TMB_H_L.Rdata')
  mean_expr_TMB_H_sign = mean_expr_TMB_H[match(diff_CpG_sign_ann$Name,mean_expr_TMB_H$cpg_site),][,c(2,3)]  # 获取差异表达cpg位点的表达谱
  
  exprSet = mean_expr_TMB_H_sign
  
  exprSet_L = melt(exprSet)
  colnames(exprSet_L) = c('probe','sample','value')
  exprSet_L$group = rep(group_list, each=nrow(exprSet))
  head(exprSet_L)
  ### ggplot2
  library(ggplot2)
  png('methy_boxplot_total_TMB_H.png',width=950,height=750,res=100)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()  
  dev.off()
  
  
  
  ## 圆形图：每条染色体上甲基化位点分布情况
  
  
  ## 金字塔图
  
  
  
}

