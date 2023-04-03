## a script to make charts shwoing difference in results from lrpath outputs, although could also be used for gsea, david, etc
#bp1 and bp2 are the two results, downloaded from the LRpath, where the columns are name, concept type, genes, coeff, oddsratio, pval, FDR, direction, siggenes
bp1 <- read.delim("./earlylate_regions/lrpathCCR_all1.txt")
bp1 <- bp1[bp1["ConceptType"] == "GOBP",]
bp1$qScore <- -log10(bp1$FDR)
down <- which(bp1$Direction == "down")
bp1[down,]$qScore <-bp1[down,]$qScore * -1
bp2 <- read.delim("./earlylate_regions/lrpathCCR_noearlylate1.txt")
bp2 <- bp2[bp2["ConceptType"] == "GOBP",]
bp2$qScore <- -log10(bp2$FDR)
down2 <- which(bp2$Direction == "down")
bp2[down2,]$qScore <-bp2[down2,]$qScore * -1
bp <- merge(bp2, bp1, suffixes = c(".not_late", ".all"), by.x = "Name", by.y = "Name")
bp$diff <- (bp$Coeff.not_late - bp$Coeff.all)  
bp <- bp[order(-(abs(bp$diff))),]
colors <- c("without IMU cofactor" = "green", "with IMU cofactor" = "red")
diff <- ggplot(head(bp, 50), aes(y =fct_rev(fct_inorder(Name)))) + 
  geom_bar(stat = "identity", aes(x = qScore.all, fill = "all"),  alpha = 0.5) +
  geom_bar(stat = "identity", aes(x = qScore.not_late, fill = "excluding late genes"  ),  width = 0.4, alpha = 0.5) + xlab("qScores, pos = enriched in e.bin 1, down = enriched in e.bin 0") +
  ylab("Biological Process")  + labs(fill = "Legend") + scale_color_manual(values = colors) +
  ggtitle("In order of most biggest difference between groups")
bp <- bp[order(-abs(bp$qScore.all)),]
old <- ggplot(head(bp, 20), aes(y =fct_rev(fct_inorder(Name)))) + 
  geom_bar(stat = "identity", aes(x = qScore.all, fill = "all"),  alpha = 0.5) +
  geom_bar(stat = "identity", aes(x = qScore.not_late, fill = "excluding early and late genes"  ),  width = 0.4, alpha = 0.5) + xlab("qScores, pos = enriched in IMU") +
  ylab("Biological Process")  + labs(fill = "Legend") + scale_color_manual(values = colors) +
  ggtitle("In order of most significant in control")
bp <- bp[order(-abs(bp$qScore.not_late)),]  
new  <- ggplot(head(bp, 20), aes(y =fct_rev(fct_inorder(Name)))) + 
  geom_bar(stat = "identity", aes(x = qScore.pos_vs_neg_bin1, fill = "pos_vs_neg_bin1"),  alpha = 0.5) +
  geom_bar(stat = "identity", aes(x = qScore.bin1_vs_0_intpos, fill = "bin1_vs_0_intpos"  ),  width = 0.4, alpha = 0.5) + xlab("qScores, pos = enriched in e.bin 1, down = enriched in e.bin 0") +
  ylab("Biological Process")  + labs(fill = "Legend") + scale_color_manual(values = colors) +
  ggtitle("In order of most biggest difference between groups")


ggplot(bp[bp$Name %in% intresting_terms,], aes(y =fct_rev(fct_inorder(Name)))) + 
  geom_bar(stat = "identity", aes(x = OddsRatio.all, fill = "all"),  alpha = 0.5) +
  geom_bar(stat = "identity", aes(x = OddsRatio.not_late, fill = "excluding late genes"  ),  width = 0.4, alpha = 0.5) + xlab("qScores, pos = enriched in e.bin 1, down = enriched in e.bin 0") +
  ylab("Biological Process")  + labs(fill = "Legend") + scale_color_manual(values = colors) +
  ggtitle("In order of most biggest difference between groups")     



