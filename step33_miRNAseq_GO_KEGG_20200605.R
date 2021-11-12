## --------------------------------------------------------------
## GO/KEGG for target gene from miR
## --------------------------------------------------------------
rm(list=ls())
gc()


## GO / KEGG
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

##### 准备数据
df1 = read.csv('csv/miRWalk_miRNA_Targets_01_10.csv', header = T, fill = T)
df2 = read.csv('csv/miRWalk_miRNA_Targets_11_20.csv', header = T, fill = T)
df3 = read.csv('csv/miRWalk_miRNA_Targets_21_30.csv', header = T, fill = T)
df4 = read.csv('csv/miRWalk_miRNA_Targets_31_40.csv', header = T, fill = T)
df5 = read.csv('csv/miRWalk_miRNA_Targets_41_50.csv', header = T, fill = T)
df = rbind(df1, df2)
df = rbind(df,  df3)
df = rbind(df,  df4)
df = rbind(df,  df5)
# 筛选
# data = subset(df, miRDB == 1 & TargetScan == 1 & validated != '')
data = subset(df, miRDB == 1 & TargetScan == 1)  # 两个数据库共有的靶基因
gene = unique(as.character(data$genesymbol))


# id转换，目标 ENTREZID
gene_df = bitr(gene, OrgDb = org.Hs.eg.db,
               fromType = 'SYMBOL',
               toType = 'ENTREZID')
gene = gene_df$ENTREZID

# GO
result_go <- enrichGO(gene          = gene,  # 选择所有差异表达基因进行富集分析
                      OrgDb         = "org.Hs.eg.db",
                      ont           = "all",  #同时进行'MF','BP','CC'三种分析
                      pvalueCutoff  = 0.1,   #实际上是 padj，选0.05即可
                      readable      = T)
head(result_go)[,1:6]
# 绘图
# png('plot/plot_GO_analysis_miRNAseq_total_result_go_of_targetGene.png', width = 1300, height = 1100, res=140)
dotplot(result_go, split="ONTOLOGY", font.size=15) + facet_grid(ONTOLOGY~.,scale="free")  # 三个一起画图
# 保存结果
save(result_go, file='Rdata/TCGA_LUAD_training_cohort_miRNA_DEG50_adjp0.05_FC0.35_result_go_of_targetGene_20200605.Rdata')
# write.csv(as.data.frame(result_go), file="csv/TCGA_LUAD_cohort1_miRNA_DEG47_result_go_of_targetGene_20200224.csv")  #结果中'MF','BP','CC'三种都有


# KEGG
result_kegg <- enrichKEGG(gene          = gene,
                          organism      = 'hsa',
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.1)   #官方建议，推荐！
head(result_kegg)[,1:6]
# KEGG结果可视化
# png('plot/plot_KEGG_analysis_miRNAseq_total_result_kegg_of_targetGene.png', width = 850, height = 750, res=100)
dotplot(result_kegg, font.size=15, showCategory=10)
# 保存结果
save(result_kegg, file='Rdata/TCGA_LUAD_training_cohort_miRNA_DEG50_adjp0.05_FC0.35_result_kegg_of_targetGene_20200605.Rdata')
# write.csv(as.data.frame(result_kegg), file="csv/TCGA_LUAD_cohort1_miRNA_DEG47_result_KEGG_of_targetGene_20200224.csv")


