library(DOSE)
library(enrichplot)
library(multienrichjam)
library(tidyr)
library(cowplot)
library(MASS)

gos <- read.delim("./lrpath_ebin_noCCR_raw_pval2.tsv")

bp <- gos[gos["ConceptType"] == "GOBP",]

bp$OddsRatio <- fractions(bp$OddsRatio)
bp$qScore <- -log10(bp$FDR)
down <- which(bp$Direction == "down")
bp[down,]$qScore <-bp[down,]$qScore * -1
bar <- ggplot(head(bp, 20), aes(y = fct_rev(fct_inorder(Name)), x = qScore, fill = str_count(SigGenes, ","))) +
  geom_bar(stat = "identity") + labs(y = "Process", color = "SigGenes")

bp_enrich <-  multienrichjam::enrichDF2enrichResult(
  bp,
  keyColname = ("Id"),
  pathGenes = "X.Genes",
  geneColname = "SigGenes",
  geneRatioColname = "OddsRatio",
      
  pvalueColname = "FDR",
  verbose = FALSE
)
bp_enrich <- setReadable(bp_enrich, 'org.Hs.eg.db', 'ENTREZID')

bp_enrich <- pairwise_termsim(bp_enrich)
p3 <- treeplot(bp_enrich, ncluster = 8)
p4 <- upsetplot(bp_enrich)
# top is the edgeR output
logfc <- tags2[,c("geneid", "logFC")]
logfc <- logfc %>% drop_na()
p1 <- cnetplot(bp_enrich, foldChange = tibble::deframe(test), color_edge = TRUE, showCategory = 40, node_label = "category")
p2<- emapplot(bp_enrich)
p5 <- emapplot(bp_enrich, showCategory = 40, node_label = 'group')

cowplot::plot_grid(p3, p4, bar, ncol=2, labels=LETTERS[1:5])
