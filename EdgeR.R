library(edgeR)
library(dplyr)


# hg19_gene_length<-read.csv('E:/all_lab_work/reference/hg19/hg19_gene_length_all_exons.txt',sep=' ')
# colnames(hg19_gene_length)<-c('length','gene_id')
# 
# hg38_gene_length<-read.csv('E:/all_lab_work/daily_script_useful_part/extract_exons_gene_length_for_FPKM/hg38_gene_length_all_exons.txt',sep=' ')
# colnames(hg38_gene_length)<-c('length','gene_id')



#' creates a edgeR object and filters for genes which are expressed more than 0.5 times per million in more than 1 sample
#' @param raw_count A genes x samples dataframe, genes are rows, cols are samples. Can accept other inputs per the edgerDGE list function  
#' @param genes a list of gene names if not using standard dataframe for raw_count. defaults to rows in raw_count
#' @param meta_df the metadata
#' @param group_name a string representing the comparison group in the meta_df
#' @param group_factor a factor representing the group
#' 
edgeRstep<-function(raw_count,genes = row.names(raw_count)){
  if (hasArg(group_name)){
    group_factor <- factor(meta_df[[group_name]])
  }
  DGEList.obj <- DGEList(counts=raw_count,genes = genes)
  isexpr <- rowSums(cpm(DGEList.obj)>=0.5) >=1
  DGEList.obj <- DGEList.obj[isexpr,]
  DGEList.obj$samples$lib.size <- colSums(DGEList.obj$counts)
  DGEList.obj <- calcNormFactors(DGEList.obj, method="upperquartile")
  genes_kept<-DGEList.obj$genes
  return(DGEList.obj)
}


# a helper function to create design matrix.  Converts all string cols to categoricals
# note this doesn't work for the e bin category because the levels are 1 and 0, which are ints

make_design<-function(metadata, group){
  metadata <- metadata%>%mutate_if(is.character, as.factor)
  formula <- "~"
  for (covar in colnames(meta_data)){
    if (covar != group)
    formula <- paste(formula,covar, "+")
  }
  formula <- paste(fomula, group)
  return(model.matrix(formula, metadata))
}



isEntrez<-function(gene_ids){
  # figure out if the gene labels are entrez ids
  if (((substr(gene_ids[1], start = 1, stop = 3) == "ENS") & (substr(gene_ids[2], start = 1, stop = 3) == "ENS")) 
    | ((strtoi(gene_ids[1]) == gene_ids[1]) &  (strtoi(gene_ids[2]) == gene_ids[2]))){
    return(TRUE)
  }
  return(FALSE)
}



edgeR_differential_study<-function(counts, metadata, group, plots = TRUE, GO = TRUE){
  if (GO & !row.names(counts)){
    print(rename)
    
  }
  edge_R_object <- edgeRstep(counts)
  design <- make_formula(metadata, group)
  y <- estimateDisp(edgeR_object,design,tagwise=TRUE, robust=TRUE)
  if (plots){
    plotBCV(y)
  }
  fit <- glmQLFit(y, design)
  if (plots) {
    plotQLDisp(fit)
  }
  qlf <- glmQLFTest(fit)
  qlf_ranked <- topTags(qlf)
  if (GO){
    return(GOterms)
  }
  return(qlf_ranked)
}