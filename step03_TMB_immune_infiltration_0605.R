## -------------------------------------------------
## tumor immune infiltration  肿瘤免疫浸润
## -------------------------------------------------
## 备注：用网页工具TIMER获得的每个患者的8种免疫细胞浸润情况进行分析
rm(list=ls())
gc()

## 
## step0: 使用上一次分析使用的immune score数据
##
if(F){
  ## expr prepare
  load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
  immune_expr = read.csv('csv/immune_estimation_440patients_overlap_in_3omics_20200606.csv', header = T, row.names = 1)
  immune_expr_TMBH = immune_expr[rownames(immune_expr) %in% pid_TMBH,]
  immune_expr_TMBL = immune_expr[rownames(immune_expr) %in% pid_TMBL,]

  
  ## t.test
  p_B_cell     = t.test(immune_expr_TMBH$B_cell,     immune_expr_TMBL$B_cell)$p.value
  p_T_cell.CD4 = t.test(immune_expr_TMBH$T_cell.CD4, immune_expr_TMBL$T_cell.CD4)$p.value
  p_T_cell.CD8 = t.test(immune_expr_TMBH$T_cell.CD8, immune_expr_TMBL$T_cell.CD8)$p.value
  p_Neutrophil = t.test(immune_expr_TMBH$Neutrophil, immune_expr_TMBL$Neutrophil)$p.value
  p_Macrophage = t.test(immune_expr_TMBH$Macrophage, immune_expr_TMBL$Macrophage)$p.value
  p_DC         = t.test(immune_expr_TMBH$DC,         immune_expr_TMBL$DC)$p.value
  
  
  # visualization1: barplot   用graphpad绘制
  if(F){
    write.csv(immune_expr_TMBH, file='csv/immune_estimation_148TMBH_overlap_3omics_6.csv')
    write.csv(immune_expr_TMBL, file='csv/immune_estimation_292TMBL_overlap_3omics_6.csv')
  }
  
  
  # visualization2：heatmap
  if(F){
    # prepare expr
    data = rbind(immune_expr_TMBH, immune_expr_TMBL)
    data = as.data.frame(t(data))  ## 转置
    data = data[c(2,4),]
    # data = scale(data)
    # data = t(scale(t(log2(data + 1))))
    data[data < 0.1] = 0.1
    data[data > 0.18] = 0.18
    # range(data)
    
    # prepare group_list
    load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
    group_list = ifelse(colnames(data) %in% pid_TMBH, 'TMB-H', 'TMB-L')
    annotation_col = data.frame(group_list=group_list)
    rownames(annotation_col)=colnames(data)
    # 绘图
    # png(paste0('plot/plot_heatmap_DEG_top50_',n,'.png'),width=750,height=750,res=100)
    pheatmap(data,
             show_colnames = F,  # 是否显示列名
             show_rownames = T,  # 是否显示行名
             cluster_rows = T,  # 是否对行进行聚类
             cluster_cols = F,  # 是否对列进行聚类
             annotation_col = annotation_col)
  }
    
  
  # visualization3：correlation plot
  if(F){
    # 准备数据
    data = rbind(immune_expr_TMBH, immune_expr_TMBL)
    # data = data[,-5]
    # 
    if(!require("corrplot")) BiocManager::install("corrplot",ask = F,update = F)
    library(corrplot)
    # 处理数据
    res = cor(data)
    # 绘图 -- 推荐
    corrplot(res, 
             method = 'number', # 用数值展示
             col = 'black')     # 数值的颜色是黑色
    corrplot(res, 
             add = T,  # 此图画在第一张图的上面
             type = 'upper',  # 只画右上部分
             tl.pos = 'n')    # 不画坐标轴
    
    
    # corrplot(res)   # 简单画图
    # corrplot(res, order='AOE', addCoef.col = 'black')  # 图 + 数值
    # corrplot(res, method='number', order='hclust', col = 'black')  # 只有数值
    # corrplot.mixed(res, lower.col = 'black', number.cex = 0.7)
    
  }
}






###
### step1: 导出TMBH和TMBL患者RNA表达谱的csv文件
###
## 获取expr的ensembl id
if(F){
  load('Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_overlap.Rdata')
  expr_RNA_TMBH_TMBL_overlap[1:5,1:5]
  id_ensembl = rownames(expr_RNA_TMBH_TMBL_overlap)
  write.csv(id_ensembl, file='csv/TCGA_LUAD_training_cohort_RNA_id_ensembl.csv')
  ## 将ENSEMBL id转换为 gene symbol
  if(!require('clusterProfiler')) BiocManger::install('clusterProfiler')
  if(!require('org.Hs.eg.db')) BiocManger::install('org.Hs.eg.db')
  id_ensembl = rownames(expr_RNA_TMBH_TMBL_overlap)
  id = bitr(id_ensembl, OrgDb = org.Hs.eg.db,
            fromType = 'ENSEMBL',
            toType = 'SYMBOL')
  head(id)
  length(id$SYMBOL)
  length(unique(id$SYMBOL))
  id = id[!duplicated(id[,c(2)]),]
  ## 将表达谱的ensembl id替换为gene symbol
  expr = expr_RNA_TMBH_TMBL_overlap
  expr$ENSEMBL = rownames(expr)
  expr = merge(id, expr, by='ENSEMBL', all = F)
  rownames(expr) = expr$SYMBOL
  expr = expr[,c(-1,-2)]  ## expr TMBH AND TMBL
  write.csv(expr, file='csv/TCGA_LUAD_RNA_expr_FPKM_geneSymbol_for_immune.csv')
  ## 将TMBH和TMBL患者的表达谱单独保存成csv
  load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
  expr_TMBH = expr[, pid_TMBH]
  expr_TMBL = expr[, pid_TMBL]
  write.csv(expr_TMBH, file='csv/TCGA_LUAD_RNA_expr_FPKM_TMBH_genesymbol_for_immune.csv')
  write.csv(expr_TMBL, file='csv/TCGA_LUAD_RNA_expr_FPKM_TMBL_genesymbol_for_immune.csv')
}


##
## 在网页工具TIMER中载入csv文件，计算免疫得分
##
## 此环节在网页工具TIMMER完成


##
## 将TMBH和TMBL患者的肿瘤免疫浸润得分进行可视化分析
##
rm(list=ls())
gc()

## 数据预处理
if(F){
  immune_expr_TMBH1 = read.csv('csv/estimation_matrix_immune_TMBH_LUAD_RNA_overlap.csv', header = T, row.names = 1)
  immune_expr_TMBH = as.data.frame(t(immune_expr_TMBH1))[,101:106]
  colnames(immune_expr_TMBH) = c('B_cell','C_fibroblast','CD4_T_cell','CD8_T_cell','Endothelial','Macrophage')
  
  immune_expr_TMBL1 = read.csv('csv/estimation_matrix_immune_TMBL1_LUAD_RNA_overlap.csv', header = T, row.names = 1)
  immune_expr_TMBL2 = read.csv('csv/estimation_matrix_immune_TMBL2_LUAD_RNA_overlap.csv', header = T, row.names = 1)
  immune_expr_TMBL11 = as.data.frame(t(immune_expr_TMBL1))[,101:106]
  immune_expr_TMBL21 = as.data.frame(t(immune_expr_TMBL2))[,101:106]
  colnames(immune_expr_TMBL11) = c('B_cell','C_fibroblast','CD4_T_cell','CD8_T_cell','Endothelial','Macrophage')
  colnames(immune_expr_TMBL21) = c('B_cell','C_fibroblast','CD4_T_cell','CD8_T_cell','Endothelial','Macrophage')
  immune_expr_TMBL = rbind(immune_expr_TMBL11, immune_expr_TMBL21)
  
  # difference of 2 dataset
  p_B_cell        = t.test(immune_expr_TMBH$B_cell,     immune_expr_TMBL$B_cell)$p.value
  p_C_fibroblast  = t.test(immune_expr_TMBH$C_fibroblast,     immune_expr_TMBL$C_fibroblast)$p.value
  p_CD4_T_cell    = t.test(immune_expr_TMBH$CD4_T_cell,     immune_expr_TMBL$CD4_T_cell)$p.value
  p_CD8_T_cell    = t.test(immune_expr_TMBH$CD8_T_cell,     immune_expr_TMBL$CD8_T_cell)$p.value
  p_Endothelial   = t.test(immune_expr_TMBH$Endothelial,     immune_expr_TMBL$Endothelial)$p.value
  p_Macrophage    = t.test(immune_expr_TMBH$Macrophage,     immune_expr_TMBL$Macrophage)$p.value
}


# visualization1：boxplot
if(F){
  # B_cell
  if(F){
    df1 = data.frame(score_of_abundances = immune_expr_TMBH$B_cell, group='TMBH')
    df2 = data.frame(score_of_abundances = immune_expr_TMBL$B_cell, group='TMBL')
    df = rbind(df1, df2)
    p = ggplot(df, aes(x=group,y=score_of_abundances, fill=group)) + 
      geom_boxplot()+  #箱线图
      geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('B_cell', 'pvalue', p_B_cell))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/tumor_immune_infiltration_B_cell.png', width = 8, height = 6)
  }
  
  # T_cell.CD4
  if(F){
    df1 = data.frame(score_of_abundances = immune_expr_TMBH$T_cell.CD4, group='TMBH')
    df2 = data.frame(score_of_abundances = immune_expr_TMBL$T_cell.CD4, group='TMBL')
    df = rbind(df1, df2)
    p = ggplot(df, aes(x=group,y=score_of_abundances, fill=group)) + 
      geom_boxplot()+  #箱线图
      geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('T_cell.CD4', 'pvalue', p_T_cell.CD4))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/tumor_immune_infiltration_T_cell.CD4.png', width = 8, height = 6)
  }
  
  # T_cell.CD8
  if(F){
    df1 = data.frame(score_of_abundances = immune_expr_TMBH$T_cell.CD8, group='TMBH')
    df2 = data.frame(score_of_abundances = immune_expr_TMBL$T_cell.CD8, group='TMBL')
    df = rbind(df1, df2)
    p = ggplot(df, aes(x=group,y=score_of_abundances, fill=group)) + 
      geom_boxplot()+  #箱线图
      geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('T_cell.CD8', 'pvalue', p_T_cell.CD8))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/tumor_immune_infiltration_T_cell.CD8.png', width = 8, height = 6)
  }
  
  # Neutrophil
  if(F){
    df1 = data.frame(score_of_abundances = immune_expr_TMBH$Neutrophil, group='TMBH')
    df2 = data.frame(score_of_abundances = immune_expr_TMBL$Neutrophil, group='TMBL')
    df = rbind(df1, df2)
    p = ggplot(df, aes(x=group,y=score_of_abundances, fill=group)) + 
      geom_boxplot()+  #箱线图
      geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('Neutrophil', 'pvalue', p_Neutrophil))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/tumor_immune_infiltration_Neutrophil.png', width = 8, height = 6)
  }
  
  # Macrophage
  if(F){
    df1 = data.frame(score_of_abundances = immune_expr_TMBH$Macrophage, group='TMBH')
    df2 = data.frame(score_of_abundances = immune_expr_TMBL$Macrophage, group='TMBL')
    df = rbind(df1, df2)
    p = ggplot(df, aes(x=group,y=score_of_abundances, fill=group)) + 
      geom_boxplot()+  #箱线图
      geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('Macrophage', 'pvalue', p_Macrophage))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/tumor_immune_infiltration_Macrophage.png', width = 8, height = 6)
  }
  
  # DC
  if(F){
    df1 = data.frame(score_of_abundances = immune_expr_TMBH$DC, group='TMBH')
    df2 = data.frame(score_of_abundances = immune_expr_TMBL$DC, group='TMBL')
    df = rbind(df1, df2)
    p = ggplot(df, aes(x=group,y=score_of_abundances, fill=group)) + 
      geom_boxplot()+  #箱线图
      geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('DC', 'pvalue', p_DC))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/tumor_immune_infiltration_DC.png', width = 8, height = 6)
  }
  
}   ## 绘图1：R语言
if(F){
  write.csv(immune_expr_TMBH, file='csv/immune_estimation_TMBH_overlap_in_3omics.csv')
  write.csv(immune_expr_TMBL, file='csv/immune_estimation_TMBL_overlap_in_3omics.csv')
}


# visualization2：heatmap   不好看，还是别画了
if(F){
  library(pheatmap)
  if(F){
    # step1：准备数据
    data = rbind(immune_expr_TMBH, immune_expr_TMBL)
    data = t(data)
    # step2：将数据中心化和标准化
    data = scale(data, center = T, scale = T)
    # step3：查看数据的范围，最大是多少，最小是多少
    ran = range(data)
    # step4：根据数据的范围，设置颜色分布范围，使得 白色 位于 0 的位置（这样画图最好看）
    bk <- c(seq(round(ran[1]),-0.01,by=0.01), seq(0,round(ran[2]),by=0.01))
    # step5：开始画图
    pheatmap(data,
             scale = "none",
             color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),  # 颜色分布
             breaks = bk,
             show_colnames = F,  # 是否显示列名
             show_rownames = T,  # 是否显示行名
             cluster_rows = T,  # 是否对行进行聚类
             cluster_cols = F,   # 是否对列进行聚类
             fontsize = 13,
             main = " ")
  }

  
  if(F){
    # prepare expr
    data = rbind(immune_expr_TMBH, immune_expr_TMBL)
    data = as.data.frame(t(data))  ## 转置
    data = data[-5,]
    colnames(data) = gsub('[.]', '-', colnames(data))
    data = t(scale(t(log2(data + 1))))
    data[data > 4] = 4
    range(data)
    
    # prepare group_list
    load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
    group_list = ifelse(colnames(data) %in% pid_TMBH, 'TMB-H', 'TMB-L')
    annotation_col = data.frame( group_list=group_list  )
    rownames(annotation_col)=colnames(data)
    # 绘图
    # png(paste0('plot/plot_heatmap_DEG_top50_',n,'.png'),width=750,height=750,res=100)
    pheatmap(data,
             show_colnames = F,  # 是否显示列名
             show_rownames = T,  # 是否显示行名
             cluster_rows = T,  # 是否对行进行聚类
             cluster_cols = F,  # 是否对列进行聚类
             annotation_col = annotation_col)
  }

  
}


# visualization3：correlation plot
if(F){
  # 准备数据
  data = rbind(immune_expr_TMBH, immune_expr_TMBL)
  data = data[,-5]
  # 
  if(!require("corrplot")) BiocManager::install("corrplot",ask = F,update = F)
  library(corrplot)
  # 处理数据
  res = cor(data)
  # 绘图 -- 推荐
  corrplot(res, 
           method = 'number', # 用数值展示
           col = 'black')     # 数值的颜色是黑色
  corrplot(res, 
           add = T,  # 此图画在第一张图的上面
           type = 'upper',  # 只画右上部分
           tl.pos = 'n')    # 不画坐标轴
  
  
  # corrplot(res)   # 简单画图
  # corrplot(res, order='AOE', addCoef.col = 'black')  # 图 + 数值
  # corrplot(res, method='number', order='hclust', col = 'black')  # 只有数值
  # corrplot.mixed(res, lower.col = 'black', number.cex = 0.7)
  
}





