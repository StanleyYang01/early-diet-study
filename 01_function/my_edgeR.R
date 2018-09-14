
## ---- functions to outfiles
# require tidyverse or dplyr package 
## library(tidyverse)
# out_analysis, out_prefix are from global environment
goana_output <- function(go, contrast_i, Treat_FC){
  tr_go <- topGO(go, number = Inf) %>% filter(P.Up<1E-5&P.Down<1E-5) %>% arrange(P.Up|P.Down, Ont)
  if(nrow(tr_go)!=0) {
    out_name <- paste(colnames(con)[contrast_i], collapse="-") 
    write.table(tr_go, paste(out_analysis, out_prefix, "_", out_name, "_Ont_", as.character(Treat_FC), ".txt", sep=""), sep="\t", row=F, quote=F)
  }
}
kegga_output <- function(keg, contrast_i, Treat_FC){
  tr_keg <- topKEGG(keg, truncate=50) %>% filter(P.Up<5E-2|P.Down<5E-2)
  if(nrow(tr_keg)!=0) {
    out_name <- paste(colnames(con)[contrast_i], collapse="-")
    write.table(tr_keg, paste(out_analysis, out_prefix, "_", out_name, "_KEGG_", as.character(Treat_FC), ".txt", sep=""), sep="\t", row=F, quote=F)
  }
}


## ---- my_glmTreat to loop every contrast to glmQLFTest, glmTreat, goana, kegga and output files
my_glmTreat <-function(contrast_i, fit, Treat_FC=1){
  ### , glmTreat + output
  tr <- glmTreat(fit, contrast=con[,contrast_i], lfc=log2(Treat_FC))
  tr_allTags <- topTags(tr, n=Inf)[[1]] %>% mutate(FC=sign(logFC)*2^(abs(logFC)))
  write.table(tr_allTags, paste(out_analysis, out_prefix, "_", paste(colnames(con)[contrast_i], collapse="-"), "_glmTreat_", as.character(Treat_FC), ".txt", sep=""), sep="\t", row=F, quote=F)
  ### goana (based on treatment) + output 
  go <- goana(tr, geneid = tr$genes$Entrezid, species="Mm")
  goana_output(go=go, contrast_i, Treat_FC)
  ### kegga (based on treatment) + output
  keg <- kegga(tr, geneid = tr$genes$Entrezid, species="Mm")
  kegga_output(keg=keg, contrast_i, Treat_FC)  
}


## ----anovaQLFtest--------------------------------------------------------
### Do not bother performing kegga and goana function, there multiple logFC made. It doesn't fit kegga and goana function which only takes one logFC.
### only export FDR<0.05
### need to detach package AnnotationDbi
my_anovaQLF_FDR <-function(contrast_i, fit){
  ### glmQLFTest + output
  res <- glmQLFTest(fit, contrast=con[, contrast_i])
  res_allTags <- topTags(res, n=Inf)[[1]]
  res_allTags_FC <- res_allTags %>% select(contains("logFC")) %>% mutate_all(funs(FC=sign(.)*2^(abs(.)))) %>% select(contains("_FC")) %>% filter(FDR<0.05)
  res_allTags <- cbind(res_allTags, res_allTags_FC)
  write.table(res_allTags, paste(out_analysis, out_prefix, "_", paste(colnames(con)[contrast_i], collapse="-"), "_glmQLFTest", ".txt", sep=""), sep="\t", row=F, quote=F)
}


## ---- plotMD-----
plotMD_All <- function(column, object){
  #object is dge object in edgR
  filename = paste(out_figure,"MD_plot_", as.character(column), "_", colnames(object$count)[column],  ".png", sep="")
  png(filename)
  plotMD(object, column)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
}

## ---- attach gene symbol, Entrezid, and genenames to edgeR object y$gene
symbols <- function(y, keytype="ENSEMBL") {
  library(org.Mm.eg.db)
  y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                           keytype, column="SYMBOL")  # keytype="ENSEMBL", attach gene Symbol from database to DGElist y 
  y$genes$Entrezid <- mapIds(org.Mm.eg.db, rownames(y),
                             keytype, column="ENTREZID")
  y$genes$Genename <- mapIds(org.Mm.eg.db, rownames(y),
                             keytype, column="GENENAME")
  # select(org.Mm.eg.db, keys=rownames(y$genes), columns=c("SYMBOL","GENENAME","ENTREZID"), keytype="ENSEMBL")
  detach("package:org.Mm.eg.db", unload=TRUE)
  detach("package:AnnotationDbi", unload=TRUE)
  return(y)

}

## --

