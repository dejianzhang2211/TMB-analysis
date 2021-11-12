## ----------------------------------------
## 用AUC值最大的分子组合构建模型
## ----------------------------------------
rm(list=ls())
gc()

## load expr
load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_3omics_expr_top15_15_15.Rdata')

## the coef of molecule group with max AUC
load('Rdata/TCGA_LUAD_training_cohort_result_of_LASSO_15RNA+15miR+15CpG.Rdata')
library(glmnet)
coef(cvfit, s='lambda.min')
cvfit
write.csv(as.matrix(coef(cvfit, s='lambda.min')), file='csv/TCGA_LUAD_training_cohort_result_of_lasso_3omics.csv')

## 构建模型
data = expr_TMBH_TMBL_3omics_training
data$score = -1.555696454 * data$cg02031308 -0.939485314 * data$cg03286742 -0.532855695 * data$cg04046889 -1.603385472 * data$cg12095807 -1.171295176 * data$cg16794961 -1.341848062 * data$cg24553235 + 0.203290638 * data$YBX2 + 0.000323171 * data$HLTF + 0.355814358 * data$KLC3 + 0.017454209 * data$WRNIP1 + 0.010739241 * data$CKS1B + 0.013056543 * data$RNF26 + 0.039397451 * data$ZYG11A + 0.582628142 * data$`hsa-miR-571` + 3.954182602 * data$`hsa-miR-586` + 0.068239671 * data$`hsa-miR-151b` + 0.000724033 * data$`hsa-miR-378i` + 0.25824073 * data$`hsa-miR-6727-5p` -0.731679875 * data$`hsa-miR-502-3p` -0.007119299 * data$`hsa-miR-6798-3p`
  
expr_TMBH_TMBL_3omics_training = data
save(expr_TMBH_TMBL_3omics_training, file='Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_3omics_expr_top15_15_15_with_modelScore1.Rdata')


## model score 与 TMB 数值的关系
if(F){
  ## 1. 将TMB和model score 合并在一起
  load('Rdata/TCGA_LUAD_WES_TMB.Rdata')
  # 处理TMB数据
  head(TMB_LUAD)
  colnames(TMB_LUAD)[1] = 'sample_id'
  TMB_LUAD$sample_id = substr(TMB_LUAD$sample_id, 1, 12)
  # 处理expr
  data$sample_id = rownames(data)
  # 合并
  data2 = merge(data, TMB_LUAD, by='sample_id', all=F)
  rownames(data2) = data2$sample_id
  
  
  ## 2. 绘制散点图
  library(ggplot2)
  ggplot(data2, aes(x=score, y=`LUAD-mutect`))+geom_point() + stat_smooth(method=lm) +
    theme_bw() + theme(panel.grid =element_blank())+
    theme(text=element_text(size=20))
  #该回归直线的置信区间默认的置信度是95%，我们也可以对其进行修改。
  # ggplot(data2, aes(x=score, y=`LUAD-mutect`))+ geom_point() + stat_smooth(method=lm, level=0.99)#99%的置信度
  #也可以不显示置信区间
  # ggplot(data2, aes(x=score, y=`LUAD-mutect`))+ geom_point(colour="grey60") + stat_smooth(method=lm, se=FALSE)
  
  
  ## 3. 相关性分析
  cor = cor.test(data2$score, data2$`LUAD-mutect`, alternative = 'two.sided', method = 'pearson', conf.level = 0.95)
  pvalue = cor[3]
  pvalue
  # pvalue = 3.399279e-48
  
  ## 4. 一元线性回归
  lm.sol <- lm(data2$`LUAD-mutect` ~ 1 + data2$score)
  summary(lm.sol)
  # Multiple R-squared:  0.5028,	Adjusted R-squared:  0.5011
  
}


## ROC analysis : model score 区分TMB高低的效果
if(F){
  library(pROC)
  
  auc(roc(data$outcome, data$score))  # roc函数建立ROC曲线
  
  roc_multi = roc(data$outcome, data$score)
  # 画图美化
  plot(roc_multi, print.auc=T, auc.polygon=T,
       grid=c(0.1, 0.2), grid.col=c("green", "red"), 
       max.auc.polygon=T, auc.polygon.col="skyblue", 
       print.thres=T)
  
  ## 获取 sensitivity, specificity, PPV, NPV, accuracy
  roc1 <- roc(data$outcome, data$score, percent=TRUE)
  coords(roc1, "local maximas", ret=c("threshold", "sens", "spec", "ppv", "npv", "accuracy"))
}


## -------------------------------------
## analysis for cohot2
## -------------------------------------
if(F){
  rm(list=ls())
  gc()
  
  
  ## data prepare
  if(F){
    load('Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_3omics_expr_top15_15_15.Rdata')
    
    ## 构建模型
    data = expr_TMBH_TMBL_3omics_test
    data$score = -1.555696454 * data$cg02031308 -0.939485314 * data$cg03286742 -0.532855695 * data$cg04046889 -1.603385472 * data$cg12095807 -1.171295176 * data$cg16794961 -1.341848062 * data$cg24553235 + 0.203290638 * data$YBX2 + 0.000323171 * data$HLTF + 0.355814358 * data$KLC3 + 0.017454209 * data$WRNIP1 + 0.010739241 * data$CKS1B + 0.013056543 * data$RNF26 + 0.039397451 * data$ZYG11A + 0.582628142 * data$`hsa-miR-571` + 3.954182602 * data$`hsa-miR-586` + 0.068239671 * data$`hsa-miR-151b` + 0.000724033 * data$`hsa-miR-378i` + 0.25824073 * data$`hsa-miR-6727-5p` -0.731679875 * data$`hsa-miR-502-3p` -0.007119299 * data$`hsa-miR-6798-3p`
    
    expr_TMBH_TMBL_3omics_test = data
    save(expr_TMBH_TMBL_3omics_test, file='Rdata/TCGA_LUAD_test_cohort_data_for_LASSO_3omics_expr_top15_15_15_with_modelScore.Rdata')
    
  }

  ## ROC analysis
  if(F){
    library(pROC)
    
    auc(roc(data$outcome, data$score))  # roc函数建立ROC曲线
    
    roc_multi = roc(data$outcome, data$score)
    # 画图美化
    plot(roc_multi, print.auc=T, auc.polygon=T,
         grid=c(0.1, 0.2), grid.col=c("green", "red"), 
         max.auc.polygon=T, auc.polygon.col="skyblue", 
         print.thres=T)
    
    ## 获取 sensitivity, specificity, PPV, NPV, accuracy
    roc1 <- roc(data$outcome, data$score, percent=TRUE)
    coords(roc1, "local maximas", ret=c("threshold", "sens", "spec", "ppv", "npv", "accuracy"))
    
    # threshold       -3.3655
    # sensitivity     91.1111
    # specificity     75.0000
    # ppv             65.0794
    # npv             94.2857
    # accuracy        80.4511
    
  }
  
  ## point plot of model score and TMB
  if(F){
    ## 1. 将TMB和model score 合并在一起
    load('Rdata/TCGA_LUAD_WES_TMB.Rdata')
    # 处理TMB数据
    head(TMB_LUAD)
    colnames(TMB_LUAD)[1] = 'sample_id'
    TMB_LUAD$sample_id = substr(TMB_LUAD$sample_id, 1, 12)
    # 处理expr
    data$sample_id = rownames(data)
    # 合并
    data2 = merge(data, TMB_LUAD, by='sample_id', all=F)
    rownames(data2) = data2$sample_id
    
    ## 2. 绘制散点图
    library(ggplot2)
    ggplot(data2, aes(x=score, y=`LUAD-mutect`))+geom_point() + stat_smooth(method=lm) +
      theme_bw() + theme(panel.grid =element_blank())+
      theme(text=element_text(size=20))
    #该回归直线的置信区间默认的置信度是95%，我们也可以对其进行修改。
    # ggplot(data2, aes(x=score, y=`LUAD-mutect`))+ geom_point() + stat_smooth(method=lm, level=0.99)#99%的置信度
    #也可以不显示置信区间
    # ggplot(data2, aes(x=score, y=`LUAD-mutect`))+ geom_point(colour="grey60") + stat_smooth(method=lm, se=FALSE)
    
  }
  
  
  ## 3. 相关性分析
  cor = cor.test(data2$score, data2$`LUAD-mutect`, alternative = 'two.sided', method = 'pearson', conf.level = 0.95)
  pvalue = cor[3]
  pvalue
  # pvalue = 1.194545e-14
  
  ## 4. 一元线性回归
  lm.sol <- lm(data2$`LUAD-mutect` ~ 1 + data2$score)
  summary(lm.sol)
  # Multiple R-squared:  0.3663,	Adjusted R-squared:  0.3615

  
}
