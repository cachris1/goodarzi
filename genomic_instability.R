metadata <- read.csv("../../HPV_genomes/Master R01 HNSCC file with Leo immune scores and HPV subtypes.csv")
metadata$HPV.ID.2 <- gsub("\\+", "pos", metadata$HPV.ID)
metadata$HPV.ID.2 <- gsub("\\-", "neg", metadata$HPV.ID.2)
rownames(metadata) <- metadata$HPV.ID.2
omics <- read.csv("./genomic_instatbility_score_and_DNA_methylation_score.csv", row.names = 2)
omics <- merge(omics, metadata[,c("HPV.ID.2", "Seq.ID")], by = 0)
edata <- read.csv("./edata.csv")
omics$Seq.ID <- paste("UM_", omics$Seq.ID, sep = "")
omics <- merge(omics, edata, by.x = "Seq.ID", by.y = "X")
plot(patient_methylation_score ~ as.factor(e.bin), omics)
wilcox.test(patient_methylation_score ~ as.factor(e.bin), omics)


omics <- read.csv("./genomic_instatbility_score_and_DNA_methylation_score.csv", row.names = 2)
rownames(omics) <- gsub("neg", "\\-", rownames(omics))
rownames(omics) <- gsub("pos", "\\+", rownames(omics))


UMdata <- read.csv("../../HPV_genomes/methylation_per_gene.csv")

rownames(UMdata) <- UMdata$X
rownames(metadata) <- metadata$HPV.ID

meth <- UMdata[,c(grep("methylation", colnames(UMdata)))]
UMdata$sum <- rowSums(meth)

meth<- merge(UMdata, omics, by = 0)

plots2 <- vector('list')
plots3 <- vector('list')
omics <- omics[-1]
omics <- omics[-1]
count <- 1
meth$Subgroup <- as.factor(meth$Subgroup)
res <- data.frame()
for (i in colnames(omics)){
  print(i)
  if (i != "X") {
  f <-  as.formula(paste(i, "~", "sum"))
  m <-  lm(f, meth)
  m
  m2 <- lm(f, meth[meth$HPV.integration == "pos",])
  mneg <- lm(f, meth[meth$HPV.integration == "neg",])
  mimu <- lm(f, meth[meth$Subgroup == "1",])
  mkrt <- lm(f, meth[meth$Subgroup == "2",])
  n <- summary(m2)$coefficients[2,4] 
  n <- as.character(n)
  n <- substr(n, 0, 4)
  res[i,"none"] <- summary(m)$coefficients[2,4]
  res[i,"intpos"] <- summary(m2)$coefficients[2,4]
  res[i,"intneg"] <- summary(mneg)$coefficients[2,4]
  res[i,"imu"] <- summary(mimu)$coefficients[2,4]
  res[i,"krt"] <- summary(mimu)$coefficients[2,4]
  
  
  
  
  plots2[[count]]  <- ggplot(meth, aes_(x = as.name("sum"), y = as.name(i), color = (as.name("Subgroup")))) + geom_jitter() + geom_smooth(method = "lm")
  
  plots3[[count]]  <- ggplot(meth, aes_(x = as.name("sum"), y = as.name(i), color = (as.name("HPV.integration")))) + geom_jitter() + geom_smooth(method = "lm") + 
    ggtitle(paste("pval of int neg line:", n))
  count = count + 1
  }
  }
