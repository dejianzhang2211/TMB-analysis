## ---------------------------------------------------
## DEG分析
## ---------------------------------------------------
rm(list=ls())
gc()
options(stringsAsFactors = F)



## ----------------------------------------------
## TMBH - TMBL
## ----------------------------------------------
library(limma)
## 数据预处理
load('Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_training.Rdata')

# 处理表达谱：RNA要在10%以上的患者中有表达
expr_RNA_TMBH_TMBL_training[is.na(expr_RNA_TMBH_TMBL_training)] = 0
expr1 = expr_RNA_TMBH_TMBL_training
expr2 = expr1[1,][-1,]  ## 准备输出expr
cutoff = ceiling(length(group_list_RNA_training) * 0.9)
i= 1
for (i in 1:nrow(expr1)){
  df = as.data.frame(table(as.numeric(expr1[i,])))   # 统计第1行所有数值出现的频次
  df[,1] = as.numeric(as.character(df[,1]))   # 统计最小的数值出现的频次
  if(df[1,1]==0 & df[1,2] <= cutoff){  # 如果最小的数值是0，而且0的数量超过了患者总数的90%(277 = 307 x 0.9)，即只有10%的患者有表达
    expr2 = rbind(expr2, expr1[i,]) 
  }
  temp = i %% 2000  # 查看脚本进度
  if (temp == 0){
    cat(i, '\n')
  }
  # cat(i,'\n')
}

## 数据预处理
exprset = expr2
group_list = group_list_RNA_training


## 准备3个输入矩阵
# 表达矩阵：exprset
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

DEG_RNA_TMBH_TMBL = DEG
nrow(subset(DEG_RNA_TMBH_TMBL, adj.P.Val < 0.01 & abs(logFC)>=1))

save(DEG_RNA_TMBH_TMBL, file='Rdata/TCGA_LUAD_training_cohort_RNA_DEG_TMBH_TMBL_20200603.Rdata')



## ----------------------------------------------
## 导出ensembl ID
## ----------------------------------------------
if(F){
  ## data 1 morethan 0 FPKM
  load('Rdata/TCGA_LUAD_training_cohort_RNA_DEG_TMBH_TMBL_20200603_morethan0FPKM.Rdata')
  head(DEG_RNA_TMBH_TMBL)
  DEG_RNA_TMBH_TMBL$ID = rownames(DEG_RNA_TMBH_TMBL)
  RNA = as.character(subset(DEG_RNA_TMBH_TMBL, adj.P.Val < 0.01 & abs(logFC) > 3)[,'ID'])
  write.csv(RNA, file='csv/TCGA_LUAD_training_cohort_ensembl_from_RNA_DEG_morethan0FPKM_20200604.csv')
  
  ## data 2 morethan 10 percent
  load('Rdata/TCGA_LUAD_training_cohort_RNA_DEG_TMBH_TMBL_20200604_morethan10percent.Rdata')
  head(DEG_RNA_TMBH_TMBL)
  DEG_RNA_TMBH_TMBL$ID = rownames(DEG_RNA_TMBH_TMBL)
  RNA = as.character(subset(DEG_RNA_TMBH_TMBL, adj.P.Val < 0.05 & abs(logFC) > 1)[,'ID'])
  write.csv(RNA, file='csv/TCGA_LUAD_training_cohort_ensembl_from_RNA_DEG_morethan10percent_20200604.csv')
  
}



## -----------------------------------------
## 差异表达基因可视化: heatmap, volcanoplot
## -----------------------------------------
if(F){
  rm(list=ls())
  gc()
  
  load('Rdata/TCGA_LUAD_training_cohort_RNA_DEG_TMBH_TMBL_20200603_morethan0FPKM.Rdata')
  load('Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_training.Rdata')
  
  ## 数据预处理
  DEG = DEG_RNA_TMBH_TMBL
  
  expr_RNA_TMBH_TMBL_training[is.na(expr_RNA_TMBH_TMBL_training)] = 0
  expr = expr_RNA_TMBH_TMBL_training
  
  group_list = group_list_RNA_training
  
  nrow(subset(DEG, adj.P.Val < 0.01 & abs(logFC)>=3))
  # nrow(subset(DEG, adj.P.Val < 0.05 & abs(logFC)>=1))
  p_cutoff = 0.01
  logFC_cutoff = 3
  
  
  ## heatmap
  if(F){
    library(pheatmap)
    ## 用前50个分子的表达谱
    choose_gene=head(rownames(DEG),50) ## 选前50个分子作图
    # choose_gene = substr(choose_gene,1,15)
    choose_matrix=expr[choose_gene,]
    choose_matrix[1:4,1:4]
    ## 标准化
    choose_matrix=t(scale(t(log2(choose_matrix+1))))  # 标准化
    range(choose_matrix)
    ## 设置数值的上下限，让图表好看些
    choose_matrix[choose_matrix < -2] = -2  ## 这个数值要自定义
    choose_matrix[choose_matrix > 4]  = 4
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
    DEG = subset(DEG, abs(logFC) < 50)  ## 为了作图好看，去掉logFC值太大的分子
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
    
    ## another method
    if(F){
      p=ggplot(DEG, aes(x = logFC, y=-log10(adj.P.Val),color = factor(change)))+
        scale_color_manual(values = c("#2f5688","#BBBBBB","#CC0000"))+    #定义颜色
        geom_point(alpha=1,size = 0.5,pch=16)+        #定义图片类型
        scale_x_continuous(expand = c(0, 0),limits=c(-3,3),breaks=c(0,-1,-2,1,2))+         #定义x轴
        scale_y_continuous(limits=c(0,10),breaks = seq(0,10,by = 2))+                        #定义y轴
        geom_vline(xintercept = c(-1,1), size = 0.25,lty = 2)+             #在x轴的虚线
        geom_hline(yintercept = -log10(0.05), size = 0.25,lty = 2)+        #在Y轴的虚线
        theme_bw()+ #清楚背景
        theme(
          legend.position = 'none',
          #legend.text = element_text(size = 10, face = "bold"),
          axis.line = element_blank(),
          axis.ticks.x=element_line(colour = 'black',size = 0.25), ###显示x轴刻度线
          axis.ticks.y=element_line(colour = 'black',size = 0.25), ###显示y轴刻度线
          axis.ticks.length=unit(0.075,"lines"),##设置X轴上的刻度上的标尺
          axis.text.x = element_text(size = 3,color = 'black' ,vjust = 0.5, hjust = 0.5),
          axis.text.y=element_text(size = 3, color = 'black',vjust = 0.5, hjust = 0.5),
          axis.title.x = element_text(size = 3.5, color = 'black',vjust = 0, hjust = 0.5),
          axis.title.y = element_text(size = 3.5, color = 'black',vjust = 2, hjust = 0.5),
          legend.title = element_text(face ="bold"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_rect(colour = 'black',size = 0.375),
          panel.border = element_blank())+
        labs(x = "Log2FoldChange",
             y = "-log10(Adjust P-value)")
      
    }
    
  }
  
}


