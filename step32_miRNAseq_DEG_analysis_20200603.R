## --------------------------------------------------
## DEG分析
## --------------------------------------------------
rm(list=ls())
gc()
options(stringsAsFactors = F)


## --------------------------------------------------
## TMBH  TMBL
## --------------------------------------------------
library(limma)
load('Rdata/TCGA_LUAD_miR_expr_RPM_TMBH_TMBL_training.Rdata')
expr_miR_TMBH_TMBL_training[is.na(expr_miR_TMBH_TMBL_training)] = 0


# 处理表达谱：miRNA要在10%以上的患者中有表达
expr1 = expr_miR_TMBH_TMBL_training
expr2 = expr1[1,][-1,]  ## 准备输出expr
cutoff = ceiling(length(group_list_miR_training) * 0.9)
i= 1
for (i in 1:nrow(expr1)){
  df = as.data.frame(table(as.numeric(expr1[i,])))   # 统计第1行所有数值出现的频次
  df[,1] = as.numeric(as.character(df[,1]))   # 统计最小的数值出现的频次
  if(df[1,1]==0 & df[1,2] <= cutoff){  # 如果最小的数值是0，而且0的数量超过了患者总数的90%(277 = 307 x 0.9)，即只有10%的患者有表达
    expr2 = rbind(expr2, expr1[i,]) 
  } 
  # cat(i,'\n')
}

## 数据预处理
exprset = expr2
group_list = group_list_miR_training

# 表达矩阵: exprset
dim(exprset)
# 分组矩阵: design
cat('limma analysis begin...','\n')
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprset)
# 差异比较矩阵: contrast.matrix
contrast.matrix <- makeContrasts("TMBH - TMBL",levels = design)   # 声明把tumor组跟healthy组进行比较
## 开始差异分析
##step1
fit <- lmFit(exprset,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
DEG = na.omit(tempOutput)
DEG <- DEG[order(DEG$adj.P.Val),]
# 筛选显著差异的cpg位点

DEG_miR_TMBH_TMBL = DEG
nrow(subset(DEG_miR_TMBH_TMBL, adj.P.Val < 0.05 & abs(logFC)>=0.35))  ## log(1.5) == 0.4054

save(DEG_miR_TMBH_TMBL, file='Rdata/TCGA_LUAD_training_cohort_miR_DEG_TMBH_TMBL_20200603.Rdata')



## -------------------------------------
## 导出差异表达的miRNA
## -------------------------------------
if(F){
  load('Rdata/TCGA_LUAD_training_cohort_miR_DEG_TMBH_TMBL_20200603.Rdata')
  head(DEG_miR_TMBH_TMBL)
  DEG_miR_TMBH_TMBL$ID = rownames(DEG_miR_TMBH_TMBL)
  miR_DEG = subset(DEG_miR_TMBH_TMBL, adj.P.Val < 0.05 & abs(logFC)>=0.35)[, 'ID']
  
  write.csv(miR_DEG, file='csv/TCGA_LUAD_training_cohort_miR_DEG_id.csv')
}



## -----------------------------------------
## 差异表达基因可视化: heatmap, volcanoplot
## -----------------------------------------
if(F){
  rm(list=ls())
  gc()
  
  load('Rdata/TCGA_LUAD_training_cohort_miR_DEG_TMBH_TMBL_20200603.Rdata')
  load('Rdata/TCGA_LUAD_miR_expr_RPM_TMBH_TMBL_training.Rdata')
  
  ## 数据预处理
  DEG = DEG_miR_TMBH_TMBL
  
  expr_miR_TMBH_TMBL_training[is.na(expr_miR_TMBH_TMBL_training)] = 0
  expr = expr_miR_TMBH_TMBL_training
  
  group_list = group_list_miR_training
  
  nrow(subset(DEG, adj.P.Val < 0.05 & abs(logFC)>=0.35))
  p_cutoff = 0.05
  logFC_cutoff = 0.35
  
  
  ## heatmap
  if(F){
    library(pheatmap)
    ## 用前50个分子的表达谱
    choose_gene=head(rownames(DEG),50) ## 选前50个分子作图
    choose_matrix=expr[choose_gene,]
    choose_matrix[1:4,1:4]
    ## 标准化
    choose_matrix=t(scale(t(log2(choose_matrix+1))))  # 标准化
    range(choose_matrix)
    ## 设置数值的上下限，让图表好看些
    choose_matrix[choose_matrix > 6] = 6  ## 这个数值要自定义
    ## 准备图表中的参数
    annotation_col = data.frame( group_list=group_list  )
    rownames(annotation_col)=colnames(expr)
    # 绘图
    # png(paste0('plot/plot_heatmap_DEG_top50_',n,'.png'),width=750,height=750,res=100)
    pheatmap(choose_matrix,
             show_colnames = F,  # 是否显示列名
             show_rownames = T,  # 是否显示行名
             cluster_rows = T,  # 是否对行进行聚类
             cluster_cols = F,  # 是否对列进行聚类
             annotation_col = annotation_col)
  }
  
  ## volcano plot
  if(F){
    library(ggplot2)
    DEG = subset(DEG, abs(logFC) < 20)  ## 为了作图好看，去掉logFC值太大的分子
    DEG$change = as.factor(ifelse(DEG$adj.P.Val < p_cutoff & abs(DEG$logFC) > logFC_cutoff,
                                  ifelse(DEG$logFC > logFC_cutoff ,'up','down'),'not'))
    this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                        '\nThe number of up molecular is ',nrow(DEG[DEG$change =='up',]) ,
                        '\nThe number of down molecular is ',nrow(DEG[DEG$change =='down',]))
    
    g = ggplot(data=DEG, 
               aes(x=logFC, y=-log10(adj.P.Val), color=change)) +
      geom_point() +
      theme_bw(base_size=20) + theme(panel.grid =element_blank())+
      xlab("log2 FoldChange") + ylab("-log10 adj.p-value") +
      theme(text=element_text(size=20)) +
      geom_hline(yintercept=-log10(p_cutoff),linetype=5,col="black")+
      geom_vline(xintercept = c(logFC_cutoff,- logFC_cutoff), lty = 5, col='black')+
      ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5), panel.grid = element_blank())+
      scale_colour_manual(values = c('darkgreen','grey','red')) ## corresponding to the levels(res$change)
    print(g)
    
      
    }
    
  }
  
}



source('function_draw_heatmap_volcano_pca_20191006.R')

## TMBH  TMBL
# 载入数据
load('Rdata/TCGA_LUAD_group1_miR_expr_TMBH_TMBL_3omics.Rdata')
load('Rdata/TCGA_LUAD_group1_miR_DEG_TMBH_TMBL_3omics.Rdata')
draw_h_v_p(expr_miR_TMBH_TMBL_3omics1, DEG_TMBH_TMBL, 'miRNA_TMBH_TMBL', group_list_miR, 1)



## ----------------------------------------
## 获取上调和下调的miR
## -----------------------------------------
load('Rdata/TCGA_LUAD_miR_DEG_TMBH_TMBL_3omics.Rdata.Rdata')
dim(DEG_TMBH_TMBL)
head(DEG_TMBH_TMBL)
miR_up   = as.character(subset(DEG_TMBH_TMBL, padj<0.05 & log2FoldChange>  1)$ID)
miR_down = as.character(subset(DEG_TMBH_TMBL, padj<0.05 & log2FoldChange< -1)$ID)
save(miR_up,miR_down, file='Rdata/TCGA_LUAD_miRNA_DEG_up_down_TMBH_TMBL_3omics.Rdata')
write.csv(miR_up, file='csv/TCGA_LUAD_miRNA_up_TMBH_TMBL_3omics.csv')
write.csv(miR_down, file='csv/TCGA_LUAD_miRNA_down_TMBH_TMBL_3omics.csv')

if(F){
  load('Rdata/TCGA_LUAD_miRNA_DEG_TMBH_norm.Rdata')
  dim(DEG_TMBH_norm)
  head(DEG_TMBH_norm)
  miR_up   = as.character(subset(DEG_TMBH_norm, padj<0.05 & log2FoldChange>  1)$ID)
  miR_down = as.character(subset(DEG_TMBH_norm, padj<0.05 & log2FoldChange< -1)$ID)
  save(miR_up,miR_down, file='Rdata/TCGA_LUAD_miRNA_DEG_up_down_TMBH_norm.Rdata')
  write.csv(miR_up, file='csv/TCGA_LUAD_miRNA_up_TMBH_norm.csv')
  write.csv(miR_down, file='csv/TCGA_LUAD_miRNA_down_TMBH_norm.csv')
  
  load('Rdata/TCGA_LUAD_miRNA_DEG_TMBL_norm.Rdata')
  dim(DEG_TMBL_norm)
  head(DEG_TMBL_norm)
  miR_up   = as.character(subset(DEG_TMBL_norm, padj<0.05 & log2FoldChange>  1)$ID)
  miR_down = as.character(subset(DEG_TMBL_norm, padj<0.05 & log2FoldChange< -1)$ID)
  save(miR_up,miR_down, file='Rdata/TCGA_LUAD_miRNA_DEG_up_down_TMBL_norm.Rdata')
  write.csv(miR_up, file='csv/TCGA_LUAD_miRNA_up_TMBL_norm.csv')
  write.csv(miR_down, file='csv/TCGA_LUAD_miRNA_down_TMBL_norm.csv')
}






