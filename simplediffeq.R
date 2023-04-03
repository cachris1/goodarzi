library(edgeR)
library(org.Hs.eg.db)
# a function I made to do differential expression testing reproducably.  The paramaters are different 
# covars and tests. all_int and all_bin are wether to subset, and int, ebin, and CCR are things to test.
diffexpr <- function (all_int, all_bin, int, ebin, CCR){
  
  timer2 <- read.csv("./timer2.csv", row.names = 1)
  if (all_int == FALSE){
    timer2 <- timer2[timer2['HPVint'] == 'pos',]
  }
  if (all_bin == FALSE){
    timer2 <- timer2[timer2['e.bin'] == 1,]
  }
  counts <- (read.csv("./data/both_counts.tsv", sep = "\t", row.names = 1))
  counts <- counts[,rownames(timer2)]
  timer2 <- timer2[colnames(counts),]
  
  object <- DGEList(counts=counts )
  isexpr <- rowSums(cpm(object)>=0.5) >=1
  object <- object[isexpr,]
  object$samples$lib.size <- colSums(object$counts)
  object <- calcNormFactors(object, method="upperquartile")
  timer2$e.bin <- as.factor(timer2$e.bin)
  batch <- as.factor(timer2$Cohort)
  bin <- factor(timer2$e.bin)
  subtype <- as.factor(timer2$subtype.CCR)
  site <- factor(timer2$Site.combined)
  stage <- factor(timer2$Stage.combined)
  if (int){
    design <- model.matrix(~batch +stage+timer2$Age +site+factor(timer2$HPVint))
  }
  if (CCR){
    design <- model.matrix(~batch+stage+timer2$Age +site+factor(timer2$HPVint)+subtype)
  }
  if(ebin){
    design <- model.matrix(~batch+stage+timer2$Age +site+bin)
  }
  print(table(timer2$e.bin))
  object <- estimateDisp(object, design, robust=TRUE)
  fit <- glmQLFit(object, design, robust = TRUE)
  qlf <- glmQLFTest(fit)
  return(qlf)
  #go <- goana(qlf, geneid = 'mapIds.org.Hs.eg.db..row.names.counts....ENTREZID....SYMBOL...')
}
o <- data.frame(mapIds(org.Hs.eg.db, row.names(counts), 'ENTREZID', 'SYMBOL', multiVals = 'first'))
qlf <- diffexpr(all_int = FALSE, all_bin = TRUE, int = FALSE, ebin = TRUE, CCR = FALSE )
lrpath <- data.frame(topTags(qlf, n = nrow(qlf)))
o <- plyr::rename(o, c("mapIds.org.Hs.eg.db..row.names.counts....ENTREZID....SYMBOL..." = "ENTREZ_ID" ))
lrpath <- merge(o, lrpath, by = 0)
lrpath<- lrpath %>% drop_na()
lrpath <- lrpath[,c('ENTREZ_ID', "PValue", "logFC")]
write_tsv(lrpath, "/mnt/c/Users/Chris/Downloads/int_all_diffexpr_rawPval.txt")


