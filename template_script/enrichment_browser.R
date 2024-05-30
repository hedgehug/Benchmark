options(warn=-1)
library(EnrichmentBrowser)

gmt.file<-file.path('PATHWAY_FILE') 
hsa.gs<-getGenesets(gmt.file)

data.dir <- 'DATA_DIR'
exprs.file<-file.path(data.dir,"TARGET.exprs.tab") 
cdat.file<-file.path(data.dir,"TARGET.colData.tab") 
rdat.file<-file.path(data.dir,"TARGET.rowData.tab") 
se<-readSE(exprs.file,cdat.file,rdat.file)


se<-deAna(se,de.method="DESeq2", filter.by.expr=FALSE)

method_idx <- 1
res_list <- list()

while (method_idx<=length(EnrichmentBrowser::sbeaMethods())) {
  if (method_idx == 5){
    method_idx <- method_idx + 1
    next
  }
  sbea.res<-sbea(method=EnrichmentBrowser::sbeaMethods()[method_idx],se=se,gs=hsa.gs, alpha = 1)
  res_list[[length(res_list)+1]] <- sbea.res
  # [] <- append(res_list, sbea.res)
  method_idx <- method_idx + 1
}

comb.res <- combResults(res_list)
res <- gsRanking(comb.res, signif.only = FALSE)
write.table(res, 'OUTPUT_FILE', sep = '\t', quote = FALSE)
