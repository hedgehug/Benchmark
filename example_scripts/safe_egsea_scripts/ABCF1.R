library('EGSEA')
library('EGSEAdata')
library(limma)

# load an example data first
count_mat <- read.table('../data/count_matrix/ABCF1.txt_data_matrix.txt',
                        header=TRUE, sep = '\t', row.names = 1)
cnt_mat <- data.matrix(count_mat[2:5])

# read the gmt files
pathway_file <- read.table('../pathway_files/K562_DESEq2.gmt',
                           header=FALSE, sep = '\t')
pathway_name <- unlist(pathway_file[1], use.names = FALSE)
pathway_info <- pathway_file[3:ncol(pathway_file)]

# list to df first, then to list
df <- data.frame(matrix(unlist(pathway_info), nrow=length(pathway_info), byrow=TRUE))
pathway_gene_list <- as.list(df)
names(pathway_gene_list) <- pathway_name

# build custom index
gs.annot=buildCustomIdx(geneIDs=rownames(cnt_mat),
                        gsets=pathway_gene_list, species="human", label="custom")
group <- c(1,1,2,2)
group_factor <- factor(group)
genes <- factor(rownames(cnt_mat))
final_genes <- list('FeatureID'=genes, 'Symbols'=genes)
gsa=egsea.cnt(counts=cnt_mat,group=group_factor,
              gs.annots=gs.annot,symbolsMap=NULL, baseGSEAs=egsea.base()[3],
              sort.by="avg.rank", num.threads=4,report=FALSE)
# write results
write.table(gsa$results$custom$test.results$'2vs1', '../example_results/safe_egsea_results/ABCF1.egsea.safe.txt', quote = FALSE, sep = "\t")
