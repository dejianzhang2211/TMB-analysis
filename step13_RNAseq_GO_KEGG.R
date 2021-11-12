## --------------------------------------------------------------
## GO / KEGG / PPI for DEG of RNAseq between TMBH and TMBL
## --------------------------------------------------------------
rm(list=ls())
gc()

## ------------------------------------------
## GO / KEGG
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# 载入数据
load('Rdata/TCGA_LUAD_training_cohort_RNA_DEG_TMBH_TMBL_20200603_morethan0FPKM.Rdata')

# DEG  total
if(F){
  # 将 ensemble id 转换成 ENTREZID
  DEG_RNA_TMBH_TMBL$ID = rownames(DEG_RNA_TMBH_TMBL)
  gene = as.character(subset(DEG_RNA_TMBH_TMBL, adj.P.Val < 0.01 & abs(logFC) >= 3)[,'ID'])
  gene_df = bitr(gene, OrgDb = org.Hs.eg.db,
                 fromType = 'ENSEMBL',
                 toType = 'ENTREZID')
  gene = gene_df$ENTREZID
  
  # GO
  result_go <- enrichGO(gene          = gene,  # 选择所有差异表达基因进行富集分析
                        OrgDb         = "org.Hs.eg.db",
                        ont           = "all",  #同时进行'MF','BP','CC'三种分析
                        pvalueCutoff  = 0.5,   #实际上是 padj，选0.05即可
                        readable      = T)
  head(result_go)[,1:6]
  # 保存结果
  save(result_go, file='Rdata/TCGA_LUAD_cohort1_RNA_DEG_pvalue01_FC1_result_go_total.Rdata')
  write.csv(as.data.frame(result_go), file="csv/TCGA_LUAD_cohort1_RNA_DEG_pvalue01_FC1_result_go_total.csv")  #结果中'MF','BP','CC'三种都有
  # 绘图
  # png('plot/result_GO_analysis_RNAseq_DEG_total_TMBH_TMBL.png', width = 1300, height = 1100, res=140)
  dotplot(result_go, split="ONTOLOGY", font.size=18) + facet_grid(ONTOLOGY~.,scale="free")  # 三个一起画图
  dev.off()
  
  # KEGG
  result_kegg <- enrichKEGG(gene          = gene,
                            organism      = 'hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.5)   #官方建议，推荐！
  head(result_kegg)[,1:6]
  # 保存结果
  save(result_kegg, file='Rdata/TCGA_LUAD_cohort1_RNA_DEG_pvalue01_FC1_result_kegg_total.Rdata')
  write.csv(as.data.frame(result_kegg), file="Rdata/TCGA_LUAD_cohort1_RNA_DEG_pvalue01_FC1_result_kegg_total.csv")
  # KEGG结果可视化
  # png('plot/result_KEGG_analysis_RNAseq_DEG_total_TMBH_TMBL.png', width = 850, height = 750, res=100)
  dotplot(result_kegg, font.size=18, showCategory=10)
  dev.off()
}









