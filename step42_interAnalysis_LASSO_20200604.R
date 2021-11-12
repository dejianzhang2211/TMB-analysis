## -----------------------------------------
## LASSO 机器学习
## -----------------------------------------
## top 45, top 45, top 45, 15+15+15

rm(list=ls())
gc()

library(pROC)
library(ggplot2)
library(glmnet)



## ----------------------------------------
## RNA
## ----------------------------------------
if(F){
  rm(list=ls())
  gc()
  
  # 准备marker
  load('Rdata/TCGA_LUAD_training_cohort_RNA_DEG_TMBH_TMBL_20200603_morethan0FPKM.Rdata')
  # RNA = as.character(subset(DEG_TMBH_TMBL, adj.P.Val<0.05 & abs(logFC)>1)[,1])
  DEG_RNA_TMBH_TMBL$ID = rownames(DEG_RNA_TMBH_TMBL)
  RNA = as.character(DEG_RNA_TMBH_TMBL[1:45, 'ID'])
  write.csv(RNA, file='csv/TCGA_LUAD_training_cohort_RNA_DEG_top45_morethan0FPKM.csv')  ## 用这个csv去网页g:Profiler转换id
  # 准备expr
  load('Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_training.Rdata')
  data = expr_RNA_TMBH_TMBL_training
  data = data[rownames(data) %in% RNA,]
  data2 = as.data.frame(t(data))  # 转置
  # change name
  gene_from_RNAseq = read.csv('csv/TCGA_LUAD_training_cohort_genesymbol_from_RNA_DEG_top45_morethan0FPKM.csv')  ## 这是用前面的csv转换得来的
  colnames(data2) = gene_from_RNAseq[match(colnames(data2), gene_from_RNAseq$initial_alias), 3]  # 用gene name替换ensembl
  expr_RNA_TMBH_TMBL_3omics = data2
  dim(expr_RNA_TMBH_TMBL_3omics)
  # 添加正确结果类
  load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
  expr_RNA_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_RNA_TMBH_TMBL_3omics) %in% TMBH_group1, 'TMBH', 'TMBL')
  expr_RNA_TMBH_TMBL_3omics$outcome = ifelse(expr_RNA_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
  
  save(expr_RNA_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_RNA_expr_RPKM_top45DEG_20200604.Rdata')
  
  
  ## LASSO回归分析
  if(F){
    ## data prepare
    load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_RNA_expr_RPKM_top45DEG_20200604.Rdata')
    x = as.matrix(subset(expr_RNA_TMBH_TMBL_3omics, select = - c(sample, outcome)))  # 准备自变量
    y = expr_RNA_TMBH_TMBL_3omics[, 'outcome']  # 准备因变量
    ## logistic回归
    set.seed(123)
    cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc", nfolds = 10)
    # 可视化
    plot(cvfit)
    # 获取最优模型分别包括哪些分子，及它们的系数 coef
    coef(cvfit, s='lambda.min')
    cvfit
    save(cvfit, file='Rdata/TCGA_LUAD_training_cohort_result_of_LASSO_RNA_top45.Rdata')
  }
  

}


## ----------------------------------------
## miR 
## ----------------------------------------
if(F){
  rm(list=ls())
  gc()
  
  # 准备marker
  load('Rdata/TCGA_LUAD_training_cohort_miR_DEG_TMBH_TMBL_20200603.Rdata')
  DEG_miR_TMBH_TMBL$ID = rownames(DEG_miR_TMBH_TMBL)
  miR = as.character(DEG_miR_TMBH_TMBL[1:45, 'ID'])
  # miR = as.character(DEG_TMBH_TMBL[1:50, 'ID'])
  
  # 准备expr
  load('Rdata/TCGA_LUAD_miR_expr_RPM_TMBH_TMBL_training.Rdata')
  data = expr_miR_TMBH_TMBL_training
  data = data[rownames(data) %in% miR,]
  data2 = as.data.frame(t(data))  # 转置
  data2[is.na(data2)] = 0
  expr_miR_TMBH_TMBL_3omics = data2
  dim(expr_miR_TMBH_TMBL_3omics)
  
  # 添加正确结果类
  load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
  expr_miR_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_miR_TMBH_TMBL_3omics) %in% TMBH_group1, 'TMBH', 'TMBL')
  expr_miR_TMBH_TMBL_3omics$outcome = ifelse(expr_miR_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
  
  save(expr_miR_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_miR_expr_RPM_top45DEG_20200604.Rdata')
  
  ## LASSO
  if(F){
    ## data prepare
    load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_miR_expr_RPM_top45DEG_20200604.Rdata')
    x = as.matrix(subset(expr_miR_TMBH_TMBL_3omics, select = - c(sample, outcome)))  # 准备自变量
    y = expr_miR_TMBH_TMBL_3omics[, 'outcome']  # 准备因变量
    ## logistic回归
    set.seed(123)
    cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc", nfolds = 10)
    # 可视化
    plot(cvfit)
    # 获取最优模型分别包括哪些分子，及它们的系数 coef
    coef(cvfit, s='lambda.min')
    cvfit
    save(cvfit, file='Rdata/TCGA_LUAD_training_cohort_result_of_LASSO_miR_top45.Rdata')
  }
  
  
  ## -----------------------------------------------------------------
  ### 计算AUC
  if(F){
    
    ## 单项观测指标
    df_auc = data.frame(marker=NA, auc=NA)[-1,]
    for (i in 1:54){
      cat(i,'\n')
      auc = auc(roc(expr_miR_TMBH_TMBL_3omics$outcome, expr_miR_TMBH_TMBL_3omics[,i]))
      df_auc[i,1] = colnames(expr_miR_TMBH_TMBL_3omics)[i]
      df_auc[i,2] = as.numeric(auc)
    }
    df_auc2 = df_auc[order(-df_auc$auc),]
    
    save(expr_miR_TMBH_TMBL_3omics, df_auc2, file='Rdata/TCGA_LUAD_group1_data_for_miR_ROC.Rdata')
    
    
    ### 计算AUC - 多项观测指标
    load('Rdata/TCGA_LUAD_group1_data_for_miR_ROC.Rdata')
    data = expr_miR_TMBH_TMBL_3omics
    
    # 拟合模型
    if(F){
      model_full = glm(data$outcome ~ `hsa-miR-552-3p` + `hsa-miR-194-5p` + `hsa-miR-1-3p` + `hsa-miR-508-3p` + `hsa-miR-552-5p` + `hsa-miR-506-3p` 
                       + `hsa-miR-301b-5p` + `hsa-miR-514a-3p` + `hsa-miR-7702` + `hsa-miR-509-3-5p` + `hsa-miR-676-5p` + `hsa-miR-147b` + `hsa-miR-144-5p` 
                       + `hsa-miR-192-5p` + `hsa-miR-192-3p` + `hsa-miR-508-5p` + `hsa-miR-133a-3p` + `hsa-miR-486-5p` + `hsa-miR-676-3p` + `hsa-miR-184` 
                       + `hsa-miR-194-3p` + `hsa-miR-509-3p` + `hsa-miR-3065-3p` + `hsa-miR-202-5p` + `hsa-miR-216a-5p` + `hsa-miR-451a` + `hsa-miR-675-3p` 
                       + `hsa-miR-215-5p` + `hsa-miR-556-3p` + `hsa-miR-891a-5p` + `hsa-miR-149-5p` + `hsa-miR-495-3p` + `hsa-miR-519a-5p` + `hsa-miR-323a-3p`
                       + `hsa-miR-431-5p` + `hsa-miR-23c` + `hsa-miR-372-3p` + `hsa-miR-1269b` + `hsa-miR-496` + `hsa-miR-485-3p` 
                       + `hsa-miR-369-3p` + `hsa-miR-487b-3p` + `hsa-miR-379-3p` + `hsa-miR-935` + `hsa-miR-129-2-3p` + `hsa-miR-1304-3p` + `hsa-miR-377-3p` + `hsa-miR-539-3p`
                       + `hsa-miR-541-3p` + `hsa-miR-376a-3p` + `hsa-miR-514b-5p` + `hsa-miR-206` + `hsa-miR-514a-5p` + `hsa-miR-513c-5p`, 
                       data=data, family='binomial')
    }  # total miR
    if(F){
      model_full = glm(data$outcome ~ `hsa-miR-506-3p` + `hsa-miR-508-5p` + `hsa-miR-486-5p` + `hsa-miR-184` 
                       + `hsa-miR-301b-5p` + `hsa-miR-7702` + `hsa-miR-676-5p` + `hsa-miR-147b` 
                       + `hsa-miR-509-3p` + `hsa-miR-508-3p` , 
                       data=data, family='binomial')
    }  # group1  0.8044
    if(F){
      model_full = glm(data$outcome ~ `hsa-miR-194-5p` 
                       + `hsa-miR-215-5p` + `hsa-miR-556-3p` + `hsa-miR-149-5p`
                       + `hsa-miR-372-3p` + `hsa-miR-935`
                       + `hsa-miR-7702` + `hsa-miR-676-5p` 
                       + `hsa-miR-192-3p` + `hsa-miR-486-5p`, 
                       data=data, family='binomial')
    }  # group2  0.8330
    
    model_full = glm(data$outcome ~ `hsa-miR-194-5p` 
                     + `hsa-miR-215-5p` + `hsa-miR-556-3p` + `hsa-miR-149-5p`
                     + `hsa-miR-372-3p` + `hsa-miR-935`
                     + `hsa-miR-7702` + `hsa-miR-676-5p` 
                     + `hsa-miR-192-3p` + `hsa-miR-486-5p`, 
                     data=data, family='binomial')
    summary(model_full)
    
    model_step = step(model_full)
    summary(model_step)
    
    # 将模型代入表达谱
    data$pre_full = predict(model_full, type='response')  # 建立预测变量
    # 计算auc
    auc(roc(data$outcome, data$pre_full))  # roc函数建立ROC曲线
    # 画图
    roc_multi = roc(data$outcome, data$pre_full)
    plot(roc_multi, print.auc=T, auc.polygon=T, main='data1-10miR-miR-194-5p/miR-215-5p/miR-556-3p/miR-149-5p/
         miR-372-3p/miR-935/miR-7702/miR-676-5p/miR-192-3p/miR-486-5p',
         grid=c(0.1, 0.2), grid.col=c("green", "red"), 
         max.auc.polygon=T, auc.polygon.col="skyblue", 
         print.thres=T)
    
    
    ## -------------------------------------------
    ## 绘制top5的ROC曲线
    ## -------------------------------------------
    if(F){
      # 载入数据
      load('Rdata/TCGA_LUAD_group1_data_for_miR_ROC.Rdata')
      head(df_auc2)
      data = expr_miR_TMBH_TMBL_3omics
      # 开始绘图
      roc_1 = plot.roc(data$outcome, data$`hsa-miR-552-3p`, percent=T, col='1')
      roc_2 = lines.roc(data$outcome, data$`hsa-miR-194-5p`, percent=T, col='2')
      roc_3 = lines.roc(data$outcome, data$`hsa-miR-1-3p`, percent=T, col='3')
      roc_4 = lines.roc(data$outcome, data$`hsa-miR-508-3p`, percent=T, col='4')
      roc_5 = lines.roc(data$outcome, data$`hsa-miR-552-5p`, percent=T, col='5')
      # 添加图例
      legend('bottomright', legend=c('hsa-miR-552-3p','hsa-miR-194-5p','hsa-miR-1-3p','hsa-miR-508-3p','hsa-miR-552-5p'), 
             col=c('1','2','3','4','5'), lwd=2)
  }
  

  

  

    
    
  }
  
}


## ----------------------------------------
## cpg
## ----------------------------------------
if(F){
  rm(list=ls())
  gc()
  
  load('Rdata/TCGA_LUAD_training_cohort_cpg_DEG_TMBH_TMBL_20200603.Rdata')
  DEG_cpg_TMBH_TMBL$ID = rownames(DEG_cpg_TMBH_TMBL)
  cpg = DEG_cpg_TMBH_TMBL[1:45, 'ID']
  # 准备expr
  load('Rdata/TCGA_LUAD_cpg_expr_betaValue_TMBH_TMBL_training.Rdata')
  data = expr_cpg_TMBH_TMBL_training
  data = data[rownames(data) %in% cpg,]
  data2 = as.data.frame(t(data))  # 转置
  expr_cpg_TMBH_TMBL_3omics = data2
  dim(expr_cpg_TMBH_TMBL_3omics)
  
  # 添加正确结果类
  load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
  expr_cpg_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_cpg_TMBH_TMBL_3omics) %in% TMBH_group1, 'TMBH', 'TMBL')
  expr_cpg_TMBH_TMBL_3omics$outcome = ifelse(expr_cpg_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
  save(expr_cpg_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_cpg_expr_betaValue_top45DEG_20200604.Rdata')
  
  
  ## LASSO
  if(F){
    ## data prepare
    load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_cpg_expr_betaValue_top45DEG_20200604.Rdata')
    x = as.matrix(subset(expr_cpg_TMBH_TMBL_3omics, select = - c(sample, outcome)))  # 准备自变量
    y = expr_cpg_TMBH_TMBL_3omics[, 'outcome']  # 准备因变量
    ## logistic回归
    set.seed(123)
    cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc", nfolds = 10)
    # 可视化
    plot(cvfit)
    # 获取最优模型分别包括哪些分子，及它们的系数 coef
    coef(cvfit, s='lambda.min')
    cvfit
    save(cvfit, file='Rdata/TCGA_LUAD_training_cohort_result_of_LASSO_cpg_top45.Rdata')
  }
  
}


## ----------------------------------------
## 组合：RNA + cpg + miR
## ----------------------------------------
if(F){
  rm(list=ls())
  gc()
  
  ## 准备3个组学的top15分子的表达谱
  if(F){
    ## RNA
    if(F){
      rm(list=ls())
      gc()
      
      # 准备marker
      load('Rdata/TCGA_LUAD_training_cohort_RNA_DEG_TMBH_TMBL_20200603_morethan0FPKM.Rdata')
      # RNA = as.character(subset(DEG_TMBH_TMBL, adj.P.Val<0.05 & abs(logFC)>1)[,1])
      DEG_RNA_TMBH_TMBL$ID = rownames(DEG_RNA_TMBH_TMBL)
      RNA = as.character(DEG_RNA_TMBH_TMBL[1:15, 'ID'])
      write.csv(RNA, file='csv/TCGA_LUAD_training_cohort_RNA_DEG_top15_morethan0FPKM.csv')  ## 用这个csv去网页g:Profiler转换id
      # 准备expr
      load('Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_training.Rdata')
      data = expr_RNA_TMBH_TMBL_training
      data = data[rownames(data) %in% RNA,]
      data2 = as.data.frame(t(data))  # 转置
      # change name
      gene_from_RNAseq = read.csv('csv/TCGA_LUAD_training_cohort_genesymbol_from_RNA_DEG_top15_morethan0FPKM.csv')  ## 这是用前面的csv转换得来的
      colnames(data2) = gene_from_RNAseq[match(colnames(data2), gene_from_RNAseq$initial_alias), 3]  # 用gene name替换ensembl
      expr_RNA_TMBH_TMBL_3omics = data2
      dim(expr_RNA_TMBH_TMBL_3omics)
      # 添加正确结果类
      load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
      expr_RNA_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_RNA_TMBH_TMBL_3omics) %in% TMBH_group1, 'TMBH', 'TMBL')
      expr_RNA_TMBH_TMBL_3omics$outcome = ifelse(expr_RNA_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
      
      save(expr_RNA_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_RNA_expr_RPKM_top15DEG_20200604.Rdata')
    }
    
    ## miR
    if(F){
      rm(list=ls())
      gc()
      
      # 准备marker
      load('Rdata/TCGA_LUAD_training_cohort_miR_DEG_TMBH_TMBL_20200603.Rdata')
      DEG_miR_TMBH_TMBL$ID = rownames(DEG_miR_TMBH_TMBL)
      miR = as.character(DEG_miR_TMBH_TMBL[1:15, 'ID'])
      # miR = as.character(DEG_TMBH_TMBL[1:50, 'ID'])
      
      # 准备expr
      load('Rdata/TCGA_LUAD_miR_expr_RPM_TMBH_TMBL_training.Rdata')
      data = expr_miR_TMBH_TMBL_training
      data = data[rownames(data) %in% miR,]
      data2 = as.data.frame(t(data))  # 转置
      data2[is.na(data2)] = 0
      expr_miR_TMBH_TMBL_3omics = data2
      dim(expr_miR_TMBH_TMBL_3omics)
      
      # 添加正确结果类
      load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
      expr_miR_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_miR_TMBH_TMBL_3omics) %in% TMBH_group1, 'TMBH', 'TMBL')
      expr_miR_TMBH_TMBL_3omics$outcome = ifelse(expr_miR_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
      
      save(expr_miR_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_miR_expr_RPM_top15DEG_20200604.Rdata')
    }
    
    ## cpg
    if(F){
      rm(list=ls())
      gc()
      
      load('Rdata/TCGA_LUAD_training_cohort_cpg_DEG_TMBH_TMBL_20200603.Rdata')
      DEG_cpg_TMBH_TMBL$ID = rownames(DEG_cpg_TMBH_TMBL)
      cpg = DEG_cpg_TMBH_TMBL[1:15, 'ID']
      # 准备expr
      load('Rdata/TCGA_LUAD_cpg_expr_betaValue_TMBH_TMBL_training.Rdata')
      data = expr_cpg_TMBH_TMBL_training
      data = data[rownames(data) %in% cpg,]
      data2 = as.data.frame(t(data))  # 转置
      expr_cpg_TMBH_TMBL_3omics = data2
      dim(expr_cpg_TMBH_TMBL_3omics)
      
      # 添加正确结果类
      load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
      expr_cpg_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_cpg_TMBH_TMBL_3omics) %in% TMBH_group1, 'TMBH', 'TMBL')
      expr_cpg_TMBH_TMBL_3omics$outcome = ifelse(expr_cpg_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
      save(expr_cpg_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_cpg_expr_betaValue_top15DEG_20200604.Rdata')
      
    }
    
  }
  
  
  # 加载
  load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_RNA_expr_RPKM_top15DEG_20200604.Rdata')
  expr_RNA_TMBH_TMBL_3omics = expr_RNA_TMBH_TMBL_3omics[,c(-16, -17)]
  load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_miR_expr_RPM_top15DEG_20200604.Rdata')
  expr_miR_TMBH_TMBL_3omics = expr_miR_TMBH_TMBL_3omics[,c(-16, -17)]
  load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_cpg_expr_betaValue_top15DEG_20200604.Rdata')
  expr_cpg_TMBH_TMBL_3omics = expr_cpg_TMBH_TMBL_3omics[,c(-16, -17)]
  
  # 3个组学数据合并
  expr_cpg_TMBH_TMBL_3omics$sample = rownames(expr_cpg_TMBH_TMBL_3omics)
  expr_RNA_TMBH_TMBL_3omics$sample = rownames(expr_RNA_TMBH_TMBL_3omics)
  expr_miR_TMBH_TMBL_3omics$sample = rownames(expr_miR_TMBH_TMBL_3omics)
  expr_TMBH_TMBL_3omics = merge(expr_cpg_TMBH_TMBL_3omics, expr_RNA_TMBH_TMBL_3omics, by='sample')
  expr_TMBH_TMBL_3omics = merge(expr_TMBH_TMBL_3omics,       expr_miR_TMBH_TMBL_3omics, by='sample')
  rownames(expr_TMBH_TMBL_3omics) = expr_TMBH_TMBL_3omics$sample
  expr_TMBH_TMBL_3omics = expr_TMBH_TMBL_3omics[,-1]
  
  # 添加正确结果类
  load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
  expr_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_TMBH_TMBL_3omics) %in% TMBH_group1, 'TMBH', 'TMBL')
  expr_TMBH_TMBL_3omics$outcome = ifelse(expr_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
  expr_TMBH_TMBL_3omics_training = expr_TMBH_TMBL_3omics
  save(expr_TMBH_TMBL_3omics_training, file='Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_3omics_expr_top15_15_15.Rdata')
  
  ## LASSO 回归
  if(F){
    ## data prepare
    load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_3omics_expr_top15_15_15.Rdata')
    x = as.matrix(subset(expr_TMBH_TMBL_3omics_training, select = - c(sample, outcome)))  # 准备自变量
    y = expr_TMBH_TMBL_3omics_training[, 'outcome']  # 准备因变量
    ## logistic回归
    set.seed(123)
    cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc", nfolds = 10)
    # 可视化
    plot(cvfit)
    # 获取最优模型分别包括哪些分子，及它们的系数 coef
    coef(cvfit, s='lambda.min')
    cvfit
    save(cvfit, file='Rdata/TCGA_LUAD_training_cohort_result_of_LASSO_15RNA+15miR+15CpG.Rdata')
  }
  
  }


## ----------------------------------------------------
## 准备表达谱 test cohort
## ----------------------------------------------------
if(F){
  ## 准备3个组学的top15分子的表达谱
  if(F){
    ## RNA
    if(F){
      rm(list=ls())
      gc()
      
      # 准备marker
      load('Rdata/TCGA_LUAD_training_cohort_RNA_DEG_TMBH_TMBL_20200603_morethan0FPKM.Rdata')
      # RNA = as.character(subset(DEG_TMBH_TMBL, adj.P.Val<0.05 & abs(logFC)>1)[,1])
      DEG_RNA_TMBH_TMBL$ID = rownames(DEG_RNA_TMBH_TMBL)
      RNA = as.character(DEG_RNA_TMBH_TMBL[1:15, 'ID'])
      write.csv(RNA, file='csv/TCGA_LUAD_training_cohort_RNA_DEG_top15_morethan0FPKM.csv')  ## 用这个csv去网页g:Profiler转换id
      # 准备expr
      load('Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_test.Rdata')
      data = expr_RNA_TMBH_TMBL_test
      data = data[rownames(data) %in% RNA,]
      data2 = as.data.frame(t(data))  # 转置
      # change name
      gene_from_RNAseq = read.csv('csv/TCGA_LUAD_training_cohort_genesymbol_from_RNA_DEG_top15_morethan0FPKM.csv')  ## 这是用前面的csv转换得来的
      colnames(data2) = gene_from_RNAseq[match(colnames(data2), gene_from_RNAseq$initial_alias), 3]  # 用gene name替换ensembl
      expr_RNA_TMBH_TMBL_3omics = data2
      dim(expr_RNA_TMBH_TMBL_3omics)
      # 添加正确结果类
      load('Rdata/TCGA_LUAD_pid_cohort_test.Rdata')
      expr_RNA_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_RNA_TMBH_TMBL_3omics) %in% TMBH_group2, 'TMBH', 'TMBL')
      expr_RNA_TMBH_TMBL_3omics$outcome = ifelse(expr_RNA_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
      
      save(expr_RNA_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_RNA_expr_RPKM_top15DEG_20200604.Rdata')
    }
    
    ## miR
    if(F){
      rm(list=ls())
      gc()
      
      # 准备marker
      load('Rdata/TCGA_LUAD_training_cohort_miR_DEG_TMBH_TMBL_20200603.Rdata')
      DEG_miR_TMBH_TMBL$ID = rownames(DEG_miR_TMBH_TMBL)
      miR = as.character(DEG_miR_TMBH_TMBL[1:15, 'ID'])
      # miR = as.character(DEG_TMBH_TMBL[1:50, 'ID'])
      
      # 准备expr
      load('Rdata/TCGA_LUAD_miR_expr_RPM_TMBH_TMBL_test.Rdata')
      data = expr_miR_TMBH_TMBL_test
      data = data[rownames(data) %in% miR,]
      data2 = as.data.frame(t(data))  # 转置
      data2[is.na(data2)] = 0
      expr_miR_TMBH_TMBL_3omics = data2
      dim(expr_miR_TMBH_TMBL_3omics)
      
      # 添加正确结果类
      load('Rdata/TCGA_LUAD_pid_cohort_test.Rdata')
      expr_miR_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_miR_TMBH_TMBL_3omics) %in% TMBH_group2, 'TMBH', 'TMBL')
      expr_miR_TMBH_TMBL_3omics$outcome = ifelse(expr_miR_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
      
      save(expr_miR_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_miR_expr_RPM_top15DEG_20200604.Rdata')
    }
    
    ## cpg
    if(F){
      rm(list=ls())
      gc()
      
      load('Rdata/TCGA_LUAD_training_cohort_cpg_DEG_TMBH_TMBL_20200603.Rdata')
      DEG_cpg_TMBH_TMBL$ID = rownames(DEG_cpg_TMBH_TMBL)
      cpg = DEG_cpg_TMBH_TMBL[1:15, 'ID']
      # 准备expr
      load('Rdata/TCGA_LUAD_cpg_expr_betaValue_TMBH_TMBL_test.Rdata')
      data = expr_cpg_TMBH_TMBL_test
      data = data[rownames(data) %in% cpg,]
      data2 = as.data.frame(t(data))  # 转置
      expr_cpg_TMBH_TMBL_3omics = data2
      dim(expr_cpg_TMBH_TMBL_3omics)
      
      # 添加正确结果类
      load('Rdata/TCGA_LUAD_pid_cohort_test.Rdata')
      expr_cpg_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_cpg_TMBH_TMBL_3omics) %in% TMBH_group2, 'TMBH', 'TMBL')
      expr_cpg_TMBH_TMBL_3omics$outcome = ifelse(expr_cpg_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
      save(expr_cpg_TMBH_TMBL_3omics, file='Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_cpg_expr_betaValue_top15DEG_20200604.Rdata')
      
    }
    
  }
  
  
  # 加载
  load('Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_RNA_expr_RPKM_top15DEG_20200604.Rdata')
  expr_RNA_TMBH_TMBL_3omics = expr_RNA_TMBH_TMBL_3omics[,c(-16, -17)]
  load('Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_miR_expr_RPM_top15DEG_20200604.Rdata')
  expr_miR_TMBH_TMBL_3omics = expr_miR_TMBH_TMBL_3omics[,c(-16, -17)]
  load('Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_cpg_expr_betaValue_top15DEG_20200604.Rdata')
  expr_cpg_TMBH_TMBL_3omics = expr_cpg_TMBH_TMBL_3omics[,c(-16, -17)]
  
  # 3个组学数据合并
  expr_cpg_TMBH_TMBL_3omics$sample = rownames(expr_cpg_TMBH_TMBL_3omics)
  expr_RNA_TMBH_TMBL_3omics$sample = rownames(expr_RNA_TMBH_TMBL_3omics)
  expr_miR_TMBH_TMBL_3omics$sample = rownames(expr_miR_TMBH_TMBL_3omics)
  expr_TMBH_TMBL_3omics = merge(expr_cpg_TMBH_TMBL_3omics, expr_RNA_TMBH_TMBL_3omics, by='sample')
  expr_TMBH_TMBL_3omics = merge(expr_TMBH_TMBL_3omics,       expr_miR_TMBH_TMBL_3omics, by='sample')
  rownames(expr_TMBH_TMBL_3omics) = expr_TMBH_TMBL_3omics$sample
  expr_TMBH_TMBL_3omics = expr_TMBH_TMBL_3omics[,-1]
  
  # 添加正确结果类
  load('Rdata/TCGA_LUAD_pid_cohort_test.Rdata')
  expr_TMBH_TMBL_3omics$sample = ifelse(rownames(expr_TMBH_TMBL_3omics) %in% TMBH_group2, 'TMBH', 'TMBL')
  expr_TMBH_TMBL_3omics$outcome = ifelse(expr_TMBH_TMBL_3omics$sample == 'TMBH', 1, 0)
  expr_TMBH_TMBL_3omics_test = expr_TMBH_TMBL_3omics
  save(expr_TMBH_TMBL_3omics_test, file='Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_3omics_expr_top15_15_15.Rdata')
}
