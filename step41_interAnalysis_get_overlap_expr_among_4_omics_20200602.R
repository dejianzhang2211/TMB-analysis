## ----------------------------------------
## 获取同时用于4个组学数据的患者
## ----------------------------------------

rm(list=ls())
gc()
options(stringsAsFactors = F)


## ----------------------------------------
## 获取同时用于4个组学数据的患者id
## ----------------------------------------
if(F){
  load('Rdata/TCGA_LUAD_pid_cpg.Rdata')
  load('Rdata/TCGA_LUAD_pid_miR.Rdata')
  load('Rdata/TCGA_LUAD_pid_RNA.Rdata')
  load('Rdata/TCGA_LUAD_pid_WES.Rdata')
  
  pid_4_omics = intersect(LUAD_pid_cpg, LUAD_pid_miR)
  pid_4_omics = intersect(LUAD_pid_RNA, pid_4_omics)
  pid_4_omics = intersect(LUAD_pid_WES, pid_4_omics)
  save(pid_4_omics, file='Rdata/TCGA_LUAD_patients_id_in_4_omics.Rdata')
  
  ## 绘制韦恩图
  library(VennDiagram)
  venn.diagram(list(Methylation=LUAD_pid_cpg, `RNA-seq`=LUAD_pid_RNA, `miRNA-seq`=LUAD_pid_miR, `WES`=LUAD_pid_WES),
               filename='veen_patient_in_4omics_20200222.tiff',
               lwd=1,#圈线粗度
               lty=1, #圈线类型
               col=c('#0099CC','#FF6666','#FFCC99','#0099CC'), #圈线颜色
               fill=c('#0099CC','#FF6666','#FFCC99','#0099CC'), #填充颜色
               cat.col=c('#0099CC','#FF6666','#FFCC99','#0099CC'),#A和B的颜色
               cat.cex = 1.2,# A和B的大小
               rotation.degree = 0,#旋转角度
               # main = "A&B&C&D",#主标题内容
               # main.cex = 2,#主标题大小
               # sub = "plot : example",#亚标题内容
               # sub.cex = 1,#亚标题字大小
               cex=1.5,#里面交集字的大小
               alpha = 0.5,#透明度 
               reverse=TRUE)
  
  
  ## 查看4个组学都有的患者TMB高低的情况
  load('Rdata/TCGA_LUAD_WES_sample_list_TMB_H_L.Rdata')
  pid_TMBH = intersect(sample_list_TMB_H, pid_4_omics)
  pid_TMBL = intersect(sample_list_TMB_L, pid_4_omics)
  save(pid_TMBH, pid_TMBL, file='Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
  
  
  ## 柱状图
  if(F){
    load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
    
    # 数据准备
    df = data.frame(group=c('TMB-H','TMB-L'), number=c(148,292))
    library(ggplot2)
    ggplot(df, aes(x=group,y=number,fill=group))+
      geom_bar(stat="identity",position="identity")+
      theme(text=element_text(size=15))+
      theme(axis.text = element_text(colour='black'))+  # 坐标轴刻度颜色
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"))+  # 删去背景和网格线
      theme(text=element_text(size=20)) +   # 所有字体增大
      geom_text(aes(label = number, vjust = -0.8, hjust = 0.5))   ## 显示柱条上的数字
  }
}


## ----------------------------------------
## 将同时拥有4个组学数据的患者随机分成训练集和测试集（7:3）
## ----------------------------------------
if(F){
  load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
  
  # 将TMBH组拆分成随机的2组
  length(pid_TMBH)
  set.seed(123)
  group1 = sample(1:148, 103, replace = F)  # 随机抽样
  TMBH_group1 = pid_TMBH[group1]  ## group1
  TMBH_group2 = setdiff(pid_TMBH, TMBH_group1)  ## group2
  # 将TMBL组拆分成随机的2组
  length(pid_TMBL)
  set.seed(123)
  group2 = sample(1:292, 204, replace = F)  # 随机抽样
  TMBL_group1 = pid_TMBL[group2]
  TMBL_group2 = setdiff(pid_TMBL, TMBL_group1)
  save(TMBH_group1, TMBL_group1, file='Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
  save(TMBH_group2, TMBL_group2, file='Rdata/TCGA_LUAD_pid_cohort_test.Rdata')
}


## -----------------------------------------------------------
## 为total，training cohort，test cohort构建表达谱
## -----------------------------------------------------------
if(F){
  ## miRNA-seq
  if(F){
    rm(list=ls())
    gc()
    options(stringsAsFactors = F)
    
    ## prepare patients id
    load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
    load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
    load('Rdata/TCGA_LUAD_pid_cohort_test.Rdata')
    
    ## prepare expr
    load('Rdata/TCGA_LUAD_miR_expr_RPM_from_Xena.Rdata')
    expr = expr_miR_RPM_total
    group_list = ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
    table(group_list)
    expr = as.data.frame(expr)
    expr_tumor  = expr[,group_list=='tumor']  ## 仅保留癌症样品
    dim(expr_tumor)
    colnames(expr_tumor) = substr(colnames(expr_tumor),1,12)  ## 去掉患者id之外的其他信息
    
    ## overlap patients
    expr_miR_TMBH_TMBL_overlap = expr_tumor[, c(pid_TMBH, pid_TMBL)]
    group_list_miR_overlap = ifelse(colnames(expr_miR_TMBH_TMBL_overlap) %in% pid_TMBH, 'TMBH', 'TMBL')
    save(expr_miR_TMBH_TMBL_overlap, group_list_miR_overlap, file='Rdata/TCGA_LUAD_miR_expr_RPM_TMBH_TMBL_overlap.Rdata')
    
    ## training cohort
    expr_miR_TMBH_TMBL_training = expr_tumor[, c(TMBH_group1, TMBL_group1)]
    group_list_miR_training = ifelse(colnames(expr_miR_TMBH_TMBL_training) %in% TMBH_group1, 'TMBH', 'TMBL')
    save(expr_miR_TMBH_TMBL_training, group_list_miR_training, file='Rdata/TCGA_LUAD_miR_expr_RPM_TMBH_TMBL_training.Rdata')
    
    ## test cohort
    expr_miR_TMBH_TMBL_test = expr_tumor[, c(TMBH_group2, TMBL_group2)]
    group_list_miR_test = ifelse(colnames(expr_miR_TMBH_TMBL_test) %in% TMBH_group2, 'TMBH', 'TMBL')
    save(expr_miR_TMBH_TMBL_test, group_list_miR_test, file='Rdata/TCGA_LUAD_miR_expr_RPM_TMBH_TMBL_test.Rdata')
  }
  
  
  ## RNA-seq
  if(F){
    rm(list=ls())
    gc()
    options(stringsAsFactors = F)
    
    ## prepare patients id
    load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
    load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
    load('Rdata/TCGA_LUAD_pid_cohort_test.Rdata')
    
    ## prepare expr
    load('Rdata/TCGA_LUAD_RNA_expr_FPKM.Rdata')
    expr = expr_RNA_FPKM
    group_list = ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
    table(group_list)
    expr = as.data.frame(expr)
    expr_tumor  = expr[,group_list=='tumor']  ## 仅保留癌症样品
    dim(expr_tumor)
    colnames(expr_tumor) = substr(colnames(expr_tumor),1,12)  ## 去掉患者id之外的其他信息
    
    ## overlap patients
    expr_RNA_TMBH_TMBL_overlap = expr_tumor[, c(pid_TMBH, pid_TMBL)]
    group_list_RNA_overlap = ifelse(colnames(expr_RNA_TMBH_TMBL_overlap) %in% pid_TMBH, 'TMBH', 'TMBL')
    save(expr_RNA_TMBH_TMBL_overlap, group_list_RNA_overlap, file='Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_overlap.Rdata')
    
    ## training cohort
    expr_RNA_TMBH_TMBL_training = expr_tumor[, c(TMBH_group1, TMBL_group1)]
    group_list_RNA_training = ifelse(colnames(expr_RNA_TMBH_TMBL_training) %in% TMBH_group1, 'TMBH', 'TMBL')
    save(expr_RNA_TMBH_TMBL_training, group_list_RNA_training, file='Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_training.Rdata')
    
    ## test cohort
    expr_RNA_TMBH_TMBL_test = expr_tumor[, c(TMBH_group2, TMBL_group2)]
    group_list_RNA_test = ifelse(colnames(expr_RNA_TMBH_TMBL_test) %in% TMBH_group2, 'TMBH', 'TMBL')
    save(expr_RNA_TMBH_TMBL_test, group_list_RNA_test, file='Rdata/TCGA_LUAD_RNA_expr_FPKM_TMBH_TMBL_test.Rdata')
  }
  
  
  ## DNA methylation
  if(F){
    rm(list=ls())
    gc()
    options(stringsAsFactors = F)
    
    ## prepare patients id
    load('Rdata/TCGA_LUAD_patients_id_in_4_omics_TMBH_TMBL.Rdata')
    load('Rdata/TCGA_LUAD_pid_cohort_training.Rdata')
    load('Rdata/TCGA_LUAD_pid_cohort_test.Rdata')
    
    ## prepare expr
    load('Rdata/TCGA_LUAD_DNA_methylation_450k_expr_clean.Rdata')
    expr = expr_luad
    group_list = ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
    table(group_list)
    expr = as.data.frame(expr)
    expr_tumor  = expr[,group_list=='tumor']  ## 仅保留癌症样品
    dim(expr_tumor)
    colnames(expr_tumor) = substr(colnames(expr_tumor),1,12)  ## 去掉患者id之外的其他信息
    
    ## overlap patients
    expr_cpg_TMBH_TMBL_overlap = expr_tumor[, c(pid_TMBH, pid_TMBL)]
    group_list_cpg_overlap = ifelse(colnames(expr_cpg_TMBH_TMBL_overlap) %in% pid_TMBH, 'TMBH', 'TMBL')
    save(expr_cpg_TMBH_TMBL_overlap, group_list_cpg_overlap, file='Rdata/TCGA_LUAD_cpg_expr_betaValue_TMBH_TMBL_overlap.Rdata')
    
    ## training cohort
    expr_cpg_TMBH_TMBL_training = expr_tumor[, c(TMBH_group1, TMBL_group1)]
    group_list_cpg_training = ifelse(colnames(expr_cpg_TMBH_TMBL_training) %in% TMBH_group1, 'TMBH', 'TMBL')
    save(expr_cpg_TMBH_TMBL_training, group_list_cpg_training, file='Rdata/TCGA_LUAD_cpg_expr_betaValue_TMBH_TMBL_training.Rdata')
    
    ## test cohort
    expr_cpg_TMBH_TMBL_test = expr_tumor[, c(TMBH_group2, TMBL_group2)]
    group_list_cpg_test = ifelse(colnames(expr_cpg_TMBH_TMBL_test) %in% TMBH_group2, 'TMBH', 'TMBL')
    save(expr_cpg_TMBH_TMBL_test, group_list_cpg_test, file='Rdata/TCGA_LUAD_cpg_expr_betaValue_TMBH_TMBL_test.Rdata')
  }
}


