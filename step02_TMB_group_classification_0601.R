## ------------------------------------------------------------------------------------
## 选择来自TCGAbiolinks的临床信息，根据TMB将患者分成TMB-H, TMB-L两组
## ------------------------------------------------------------------------------------
rm(list=ls())
gc()
options(stringsAsFactors = F)

library(survival)
library(survminer)
library(ggplot2)


# 加载整理好的，TMB生存分析数据
# load('Rdata/TCGA_LUAD_clinical_from_TCGAbiolinks_for_survival.Rdata')
# dim(clinical_luad)
load('Rdata/TCGA_LUAD_WES_TMB.Rdata')
dim(TMB_LUAD)


# 查看tumor和normal的分布情况
group_list=ifelse(as.numeric(substr(TMB_LUAD$Var1,14,15)) < 10,'tumor','normal')
table(group_list)  ## 全是tumor， 567


## 用 LUAD-mutect 将患者分组
cutoff = 10   # 设置阈值
TMB_group = ifelse(as.numeric(TMB_LUAD[,2]) >= cutoff,'TMB_H','TMB_L')
table(TMB_group)

TMB_LUAD$TMB_group = TMB_group
head(TMB_LUAD)

sample_list_TMB_H = substr(as.character(subset(TMB_LUAD, TMB_group=='TMB_H')[,1]),1,12)
sample_list_TMB_L = substr(as.character(subset(TMB_LUAD, TMB_group=='TMB_L')[,1]),1,12)


save(sample_list_TMB_H, sample_list_TMB_L, file='Rdata/TCGA_LUAD_wes_sample_list_TMB_H_L.Rdata')



## --------------------------------------------
## 柱状图：WES数据中高低TMB患者各有多少
if(F){
  # 数据准备
  df = data.frame(group=c('TMB-H','TMB-L'), number=c(184,383))
  ggplot(df, aes(x=group,y=number,fill=group))+
    geom_bar(stat="identity",position="identity")+
    theme(text=element_text(size=15))+
    theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))+  # 删去背景和网格线
    theme(text=element_text(size=20)) +   # 所有字体增大
    geom_text(aes(label = number, vjust = -0.8, hjust = 0.5))   ## 显示柱条上的数字
  
}


## --------------------------------------------
## TMB的高低与各项临床信息的相关性
## --------------------------------------------
if(F){
  ## step1：准备标注了TMB数值的临床信息谱
  load('Rdata/TCGA_LUAD_wes_TMB.Rdata')
  dim(TMB_LUAD)
  TMB_LUAD$submitter_id = substr(TMB_LUAD$Var1, 1, 12)
  data1 = TMB_LUAD[,c(6,2)]
  load('Rdata/TCGA_LUAD_clinical_from_TCGAbiolinks_for_survival.Rdata')  # 准备临床信息
  data = merge(data1, clinical_luad, by='submitter_id', all=F)  # 合并
  
  # TMB and gender
  if(F){
    num1 = data[data$gender=='male', 'LUAD-mutect']
    num2 = data[data$gender=='female', 'LUAD-mutect']
    diff = t.test(num1,num2)[[3]]
    df1 = data.frame(TMB=num1, group='male')
    df2 = data.frame(TMB=num2, group='female')
    df = rbind(df1,df2)
    p = ggplot(df, aes(x=group,y=TMB, fill=group)) + 
      geom_boxplot()+  #箱线图
      # geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('pvalue', diff))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/boxplot_TMB_with_gender.png', width = 8, height = 6)
    
  }
  
  # TMB with tumor stage
  if(F){
    num1 = data[data$tumor_stage %in% c('stage i', 'stage ia', 'stage ib'), 'LUAD-mutect']
    num2 = data[data$tumor_stage %in% c('stage ii', 'stage iia', 'stage iib'), 'LUAD-mutect']
    num3 = data[data$tumor_stage %in% c('stage iiia', 'stage iiib'), 'LUAD-mutect']
    num4 = data[data$tumor_stage %in% c('stage iv'), 'LUAD-mutect']
    diff12 = t.test(num1,num2)[[3]]
    diff13 = t.test(num1,num3)[[3]]
    diff14 = t.test(num1,num4)[[3]]
    diff23 = t.test(num2,num3)[[3]]
    diff24 = t.test(num2,num4)[[3]]
    diff34 = t.test(num3,num4)[[3]]
    
    df1 = data.frame(TMB=num1, group='stage i')
    df2 = data.frame(TMB=num2, group='stage ii')
    df3 = data.frame(TMB=num3, group='stage iii')
    df4 = data.frame(TMB=num4, group='stage iv')
    
    df = rbind(df1, df2)
    df = rbind(df,  df3)
    df = rbind(df,  df4)
    
    p = ggplot(df, aes(x=group,y=TMB, fill=group)) + 
      geom_boxplot()+  #箱线图
      # geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('2vs4 pvalue', diff24))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/boxplot_TMB_with_tumor_stage_20200120.png', width = 8, height = 6)
    
  }
  
  # TMB with age
  if(F){
    num1 = data[data$age_group=='older', 'LUAD-mutect']
    num2 = data[data$age_group=='younger', 'LUAD-mutect']
    diff = t.test(num1, num2)[[3]]
    df1 = data.frame(TMB=num1, group='older')
    df2 = data.frame(TMB=num2, group='younger')
    df = rbind(df1, df2)
    p = ggplot(df, aes(x=group,y=TMB, fill=group)) + 
      geom_boxplot()+  #箱线图
      # geom_jitter(shape=16, position=position_jitter(0.2))+  # 添加散点
      labs(title=paste('pvalue', diff))+
      theme_bw() + 
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid =element_blank())+  # 删去网格线
      theme(text=element_text(size=20))    # 所有字体增大
    ggsave(p, file='plot/boxplot_TMB_with_age.png', width = 8, height = 6)
    
  }
  
}


## ------------------------------------------------
## 查看TMBH和TMBL两组患者的年龄和性别分布情况
## ------------------------------------------------
if(F){
  load('Rdata/TCGA_LUAD_overlap_patient_in_3omics.Rdata')
  load('Rdata/TCGA_LUAD_clinical_from_TCGAbiolinks_for_survival.Rdata')
  
  group_TMBH_TMBL = c(overlap_tumor_TMBH_3omics, overlap_tumor_TMBL_3omics)
  clin_TMBH_TMBL = clinical_luad[clinical_luad$submitter_id %in% group_TMBH_TMBL,]
  clin_TMBH_TMBL$TMB_group = ifelse(clin_TMBH_TMBL$submitter_id %in% overlap_tumor_TMBH_3omics, 'TMBH', 'TMBL')
  
  # 直方图
  library(ggplot2)
  p = ggplot(clin_TMBH_TMBL, aes(x=age, fill=TMB_group)) + 
    geom_histogram(binwidth=5, position='dodge') + 
    labs(title = 'age in TMBH TMBL')+
    theme(text=element_text(size=30))    # 所有字体增大   推荐！！！
  ggsave(p, file='plot/hist_age_in_TMBH_TMBL.png',width=18, height=14, dpi=300)
  
  
  clin_TMBH = subset(clin_TMBH_TMBL, TMB_group=='TMBH')
  clin_TMBL = subset(clin_TMBH_TMBL, TMB_group=='TMBL')
  table(clin_TMBH$gender)
  table(clin_TMBL$gender)
  
  df = data.frame(group=c('TMBH','TMBH','TMBL','TMBL'), number=c(47,57,194,150), gender=c('female','male','female','male'))
  p = ggplot(df, aes(x=group,y=number,fill=gender))+
    geom_bar(stat="identity",position="dodge")+
    labs(title = 'gender in TMBH TMBL')+
    theme(text=element_text(size=30))
  ggsave(p, file='plot/barplot_gender_in_TMBH_TMBL.png', width=16, height=12, dpi=300)
  
  
}



## TMB的高低对所有患者OS的影响
if(F){
  load('Rdata/TCGA_LUAD_clinical_from_TCGAbiolinks_for_survival.Rdata')
  load('Rdata/TCGA_LUAD_wes_sample_list_TMB_H_L.Rdata')
  
  clinical1 = clinical_luad[clinical_luad$submitter_id %in% c(sample_list_TMB_H,sample_list_TMB_L),]
  head(clinical1)
  
  clinical1$TMB_group = ifelse(clinical1$submitter_id %in% sample_list_TMB_H, 'TMB-H','TMB-L')
  sfit <- survfit(Surv(month, vital_status)~TMB_group, data=clinical1)
  sfit
  # png('plot/KM_TMBH_vs_TMBL_radiation_treated.png', width=850, height=750, res=100)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
  # dev.off()
}


## TMB的高低对训练集中患者OS的影响
if(F){
  load('Rdata/TCGA_LUAD_clinical_from_TCGAbiolinks_for_survival.Rdata')
  load('Rdata/TCGA_LUAD_group1_patient_in_3omics.Rdata')
  
  clinical1 = clinical_luad[clinical_luad$submitter_id %in% c(group1_TMBH, group1_TMBL),]
  head(clinical1)
  
  clinical1$TMB_group = ifelse(clinical1$submitter_id %in% group1_TMBH, 'TMB-H','TMB-L')
  sfit <- survfit(Surv(month, vital_status)~TMB_group, data=clinical1)
  sfit
  png('plot/KM_TMBH_vs_TMBL_radiation_treated.png', width=850, height=750, res=100)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
  dev.off()
}


## chart1的患者中接受或不接受放疗对患者生存期的影响
if(F){
  rm(list=ls())
  gc()
  
  library(survival)
  library(survminer)
  # 载入数据
  load('Rdata/TCGA_LUAD_clinical_from_TCGAbiolinks_for_survival.Rdata')
  load('Rdata/TCGA_LUAD_sample_id_of_TMBH_TMBL_overlap_in_4omics_20200222.Rdata')
  
  clinical1 = clinical_luad[clinical_luad$submitter_id %in% c(list_TMBH, list_TMBL),]
  clinical2 = clinical1[clinical1$treatments_radiation_treatment_or_therapy=='yes',]

  # 绘制生存曲线
  clinical2$TMB = ifelse(clinical2$submitter_id %in% list_TMBH, 'TMBH', 'TMBL')
  sfit <- survfit(Surv(month, vital_status)~TMB, data=clinical2)
  sfit
  png('plot/KM_TMBH_vs_TMBL_radiation_treated.png', width=850, height=750, res=100)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
  dev.off()
}


# 检查每个样品分别来自那个患者
if(F){
  # 获取样品来源的患者id
  sample_id = TMB_LUAD$Var1
  name = substr(sample_id,1,12)
  name2 = unique(name)  # 去掉可能有的重复
  length(name)
  length(name2)   
  all.equal(name,name2)   # name和name2的内容完全一致，567个样品来自567个人
}

# 比较TCGAbiolinks和gdc的临床数据
if(F){
  # 获取TCGAbiolinks来源的522个患者id
  name_bio = clin_bio$submitter_id
  # 获取gdc的患者id
  name_gdc = clin_gdc$ID
  setdiff(name_bio, name_gdc)   # 结果：患者的id都是一样的
}

## 患者只接受化疗  --  total sample
if(F){
  # clin_bio中获取没有接受药物和化疗治疗的患者id
  name_bio_radia = subset(clin_bio, treatments_radiation_treatment_or_therapy=='yes')$submitter_id
  # 获取对应的TMB信息
  TMB_radia = TMB_LUAD[substr(TMB_LUAD$Var1,1,12) %in% name_bio_radia,]
  # 获取对应的临床患者
  clin_bio_radia = clin_bio[clin_bio$submitter_id %in% substr(TMB_radia$Var1,1,12),]
  
  
  # 策略：用不同的TMB作为阈值，分析生存差异的显著性，用最显著的作为阈值
  library(survival)
  library(survminer)
  # 使用LUAD-nutect
  TMB = TMB_radia
  df = data.frame(TMB_cutoff=NA, pvalue=NA)[-1,]  # 创建一个空的data.frame
  n = 0  # 设置循环起始变量
  for (i in 5:30){
    clin_bio_radia$TMB_cutoff=ifelse(as.numeric(TMB[,2])>i,'high','low')   # 将所有的样品分成两组：此基因高表达，此基因低表达
    data_survdiff=survdiff(Surv(month, vital_status)~TMB_cutoff, data=clin_bio_radia)
    pval = 1 - pchisq(data_survdiff$chisq, length(data_survdiff$n) - 1)
    n = n + 1
    df[n,1] = i
    df[n,2] = pval
  }
  df_sign = subset(df, pvalue<0.05)
  df_sign
  ## 结果：TMB等于13,14,15符合要求，选择TMB=14作为cutoff
  
  # 绘制生存曲线
  clin_bio_radia$TMB_cutoff=ifelse(as.numeric(TMB[,2])>14,'high','low')
  sfit <- survfit(Surv(month, vital_status)~TMB_cutoff, data=clin_bio_radia)
  sfit
  png('plot/KM_TMB_H_vs_L_TCGAbiolinks_cutoff14.png', width=750, height=750, res=100)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
  dev.off()
  
  # 将患者分组
  TMB_group = ifelse(as.numeric(TMB[,2]) >= 14,'TMB-H','TMB-L')
  table(TMB_group)
  
  sample_list_TMB_H = as.character((clin_bio_radia[TMB_group=='TMB-H',1]))
  sample_list_TMB_L = as.character((clin_bio_radia[TMB_group=='TMB-L',1]))
  
  save(sample_list_TMB_H, sample_list_TMB_L, file='Rdata/sample_list_2group_TMB_H_L_radiaTreat.Rdata')
}
