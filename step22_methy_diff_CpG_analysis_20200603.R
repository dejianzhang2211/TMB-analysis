## ---------------------------------------------
## DEG??????DNA methylation
## ---------------------------------------------
rm(list=ls())
gc()
options(stringsAsFactors = F)


## ----------------------------------------------
## TMBH  TMBL
## ----------------------------------------------
library(limma)
## 1. ׼??3??????????
load('Rdata/TCGA_LUAD_cpg_expr_betaValue_TMBH_TMBL_training.Rdata')
# ????Ԥ????
exprset =    expr_cpg_TMBH_TMBL_training
group_list = group_list_cpg_training  # ׼????????Ϣ


# ????????: exprset
dim(exprset)
length(group_list)

# ????????: design
cat('limma analysis begin...','\n')
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprset)
# ?????ȽϾ???: contrast.matrix
contrast.matrix <- makeContrasts("TMBH - TMBL",levels = design)   # ??????tumor????healthy?????бȽ?
## ??ʼ????????
##step1
fit <- lmFit(exprset,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##??һ??????Ҫ?????ҿ??????п???Ч??
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
diff_CpG = na.omit(tempOutput)
diff_CpG <- diff_CpG[order(diff_CpG$adj.P.Val),]

DEG_cpg_TMBH_TMBL = diff_CpG
nrow(subset(DEG_cpg_TMBH_TMBL, adj.P.Val < 0.05 & abs(logFC)>=0.15))

save(DEG_cpg_TMBH_TMBL, file='Rdata/TCGA_LUAD_training_cohort_cpg_DEG_TMBH_TMBL_20200603.Rdata')



## -----------------------------------------
## ???????????????ӻ?: heatmap, volcanoplot
## -----------------------------------------
if(F){
  rm(list=ls())
  gc()
  
  load('Rdata/TCGA_LUAD_training_cohort_cpg_DEG_TMBH_TMBL_20200603.Rdata')
  load('Rdata/TCGA_LUAD_cpg_expr_betaValue_TMBH_TMBL_training.Rdata')
  
  ## ????Ԥ????
  DEG = DEG_cpg_TMBH_TMBL
  
  expr_cpg_TMBH_TMBL_training[is.na(expr_cpg_TMBH_TMBL_training)] = 0
  expr = expr_cpg_TMBH_TMBL_training
  
  group_list = group_list_cpg_training
  
  nrow(subset(DEG, adj.P.Val < 0.05 & abs(logFC)>=0.15))
  p_cutoff = 0.05
  logFC_cutoff = 0.15
  
  
  ## heatmap
  if(F){
    library(pheatmap)
    ## ??ǰ50?????ӵı?????
    choose_gene=head(rownames(DEG),50) ## ѡǰ50????????ͼ
    choose_matrix=expr[choose_gene,]
    choose_matrix[1:4,1:4]
    ## ??׼??
    # choose_matrix=t(scale(t(log2(choose_matrix+1))))  # ??׼??
    range(choose_matrix)
    ## ??????ֵ???????ޣ???ͼ???ÿ?Щ
    # choose_matrix[choose_matrix > 6] = 6  ## ??????ֵҪ?Զ???
    ## ׼??ͼ???еĲ???
    annotation_col = data.frame( group_list=group_list  )
    rownames(annotation_col)=colnames(expr)
    # ??ͼ
    # png(paste0('plot/plot_heatmap_DEG_top50_',n,'.png'),width=750,height=750,res=100)
    pheatmap(choose_matrix,
             show_colnames = F,  # ?Ƿ???ʾ????
             show_rownames = T,  # ?Ƿ???ʾ????
             cluster_rows = T,  # ?Ƿ????н??о???
             cluster_cols = F,  # ?Ƿ????н??о???
             annotation_col = annotation_col)
  }
  
  ## volcano plot
  if(F){
    library(ggplot2)
    # DEG = subset(DEG, abs(logFC) < 20)  ## Ϊ????ͼ?ÿ???ȥ??logFCֵ̫???ķ???
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
