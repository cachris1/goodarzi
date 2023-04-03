library(edgeR)
source("./EdgeR.R")

timer2 <- read.csv("../timer2.csv", row.names = 1)


counts <- (read.csv("./both_counts.txt", sep = "\t", row.names = 1))
counts <- counts[,rownames(timer2)]
timer2 <- timer2[colnames(counts),]
object <- DGEList(counts=counts)
isexpr <- rowSums(cpm(object)>=0.5) >=1
object <- object[isexpr,]
object$samples$lib.size <- colSums(object$counts)
object <- calcNormFactors(object, method="upperquartile")
timer2$HPVint <- as.factor(timer2$HPVint)
batch <- as.factor(timer2$Cohort)
int <- factor(timer2$HPVint)
site <- factor(timer2$Site.combined)
stage <- factor(timer2$Stage.combined)  
design <- model.matrix(~batch +stage+timer2$Age+site+int)
object <- estimateDisp(object, design, robust=TRUE)
fit <- glmQLFit(object, design, robust = TRUE)
qlfint <- glmQLFTest(fit)
print(topTags(qlfint))