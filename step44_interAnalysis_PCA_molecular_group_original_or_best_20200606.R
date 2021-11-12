## ------------------------------------------------------
## PCA分析：原始组合 和 最佳组合
## ------------------------------------------------------
rm(list=ls())
gc()

## 1. 准备表达谱
load('Rdata/TCGA_LUAD_training_cohort_data_for_LASSO_3omics_expr_top15_15_15.Rdata')
expr_orign = expr_TMBH_TMBL_3omics_training[, 1:45]
expr_best  = expr_TMBH_TMBL_3omics_training[, c('cg02031308','cg03286742','cg04046889','cg12095807','cg16794961','cg24553235',
                                                'YBX2','HLTF','KLC3','WRNIP1','CKS1B','RNF26','ZYG11A',
                                                'hsa-miR-571','hsa-miR-586','hsa-miR-151b','hsa-miR-378i','hsa-miR-6727-5p','hsa-miR-502-3p','hsa-miR-6798-3p')]
# 添加group
load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
expr_orign$group = ifelse(rownames(expr_orign) %in% pid_TMBH, 'TMBH', 'TMBL')
expr_best$group  = ifelse(rownames(expr_best)  %in% pid_TMBH, 'TMBH', 'TMBL')


## 2. PCA分析：原始组合 147个分子
library(ggfortify)
df = expr_orign
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()


## 3. PCA分析：最佳组合   个Fenix
df = expr_best
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()

