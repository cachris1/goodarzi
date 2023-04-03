library("stringr")  
substr_func <- function(x)  {
  if (substr(x, 0, 2) == "LC"){
    return ("LCR");
  }
  return (substr(x, 0, 2));
}
library(cowplot)


data <- read.csv("../HPV_genomes/methylation_per_gene.csv")
meth <- data[,c(grep("methylation", colnames(data)), 57, 25)]
cov <- data[,c(grep("reads_per_cpg", colnames(data)), 57, 25)]
rna <- timer2[,c(grep("e*.cpm", colnames(timer2)))]
rna <- rna[-1]
rna$HPVint <- timer2$HPVint
rna$ID <- timer2$Sample.ID
meth2 <- meth %>% pivot_longer(cols = 1:9, names_to = "gene", values_to = "methylation")
cov2 <- cov %>% pivot_longer(cols = 1:9, names_to = "gene", values_to = "coverage")
rna2 <- rna %>% pivot_longer(cols = 1:8, names_to = "gene", values_to = "expression")
meth2$gene <- sapply(meth2$gene, substr_func)
rna2$gene <- sapply(str_to_upper(rna2$gene), substr_func)
cov2$gene <- sapply(str_to_upper(cov2$gene), substr_func)
p1 <- ggplot(meth2, aes(x = gene, y = methylation, color = HPV.integration,  group = Seq.ID)) + geom_jitter(width = 0.1) + geom_line() + xlab("HPV16 gene") + ylab("Percent methylated") + ggtitle("Average methylation per HPV gene in 13 UM samples")
p2 <- ggplot(cov2, aes(x = gene, y = log(coverage), color = HPV.integration, group = Seq.ID)) + geom_jitter(width = 0.1)+ geom_line() + xlab("HPV16 gene") + ylab("number of reads/number of CpG") + ggtitle("Coverage per HPV gene in 13 UM samples")
p3 <- ggplot(rna2, aes(x = gene, y = log(expression), color = HPVint, group = ID)) + geom_jitter(width = 0.1) + geom_line() + xlab("HPV16 gene") + ylab("Expression in log(cpm)") + ggtitle("Expression of HPV genes in 13 UM samples and 66 TCGA samples")
f1 <- cowplot::plot_grid(p1, p2, p3, labels = c("A", "B", "C"))


genes <- c("E1", "E2", "E4", "E5", "E6", "E7")
plots <- vector('list', 6)

count <- 1
plots_cov <- vector('list', 6)
plots_meth <- vector('list', 6)
plots_multi <- vector('list', 6)
plots_multi_neg <-vector('list', 6)
plots_multi_pos <- vector('list', 6)



for (i in genes){
  meth_col <- paste(i, "methylation", sep = "")
  cov_col <-  paste(i, "reads_per_cpg", sep = "")
  
  cpm_col <- paste(str_to_lower(i), ".cpm", sep = "")
  
plots[[count]] <-  ggplot(data, aes_(x = as.name(meth_col), y = as.name(cpm_col), color = as.name("HPV.integration"))) + 
  geom_jitter() + geom_smooth(method = "lm", se = FALSE) + ylab(paste(i, "cpm")) 

plots_cov[[count]] <-  ggplot(data, aes_(x = as.name(cov_col), y = as.name(cpm_col), color = as.name("HPV.integration"))) + 
  geom_jitter() + geom_smooth(method = "lm", se = TRUE) 

plots_meth[[count]] <-  ggplot(data, aes_(x = as.name(cov_col), y = as.name(meth_col), color = as.name("HPV.integration"))) + 
  geom_jitter() + geom_smooth(method = "lm", se = TRUE) 

  fm <- as.formula(paste(cpm_col, "~", cov_col, "+", meth_col, "+", "HPV.integration", sep=""))

  m <- lm(fm, data)
  print(as.name(meth_col))
  pred <- data.frame(prediction = predict(m, data), methyl = data[,(cov_col)])
  print(as.name(meth_col))
  plots_multi[[count]] <- ggplot(data, aes_(x = as.name(cov_col), y = as.name(cpm_col))) + 
    geom_jitter() + geom_line(data = pred, aes(y = prediction, x = methyl))
  
  fm2 <- as.formula(paste(cpm_col, "~", cov_col, "+", meth_col,  sep=""))
  mp <- lm(fm2, data[data$HPV.integration == "pos",])
  print(as.name(meth_col))
  pred_pos <- data.frame(prediction = predict(m, data[data$HPV.integration == "pos",]), methyl = data[data$HPV.integration == "pos",(cov_col)])
  print(as.name(meth_col))
  plots_multi_pos[[count]] <- ggplot(data[data$HPV.integration == "pos",], aes_(x = as.name(cov_col), y = as.name(cpm_col))) + 
    geom_jitter() + geom_line(data = pred, aes(y = prediction, x = methyl))
  
  mp <- lm(fm2, data[data$HPV.integration == "neg",])
  print(as.name(meth_col))
  pred_pos <- data.frame(prediction = predict(m, data[data$HPV.integration == "neg",]), methyl = data[data$HPV.integration == "neg",(cov_col)])
  print(as.name(meth_col))
  plots_multi_neg[[count]] <- ggplot(data[data$HPV.integration == "neg",], aes_(x = as.name(cov_col), y = as.name(cpm_col))) + 
    geom_jitter() + geom_line(data = pred, aes(y = prediction, x = methyl))
  
count = count + 1
}

f3 <- multiplot(plotlist = plots_cov, ncol = 3)
f4<- multiplot(plotlist = plots_meth, ncol = 3)

f2 <- multiplot(plotlist = plots, ncol = 3)


substr_func <- function(x)  {
  if (substr(x, 0, 2) == "LC"){
    return ("LCR");
  }
  return (substr(x, 0, 2));
}



