DMR_genes <- read.csv("./annotSigDMRS_IMU_v_KRT_genes_id.csv")
DMR_genes <- DMR_genes %>% drop_na() %>% distinct()
late_genes <- read.csv("earlylate replicated genes.csv")
late_genes <- late_genes$gene.LATE
late_genes <- Filter(function(x) x != "", late_genes)
DMR_genelist <- as.character(DMR_genes$annot.gene_id)
keys <- data.frame(mapIds(org.Hs.eg.db, late_genes,  keytype = 'SYMBOL', column = 'ENTREZID', multiVals = 'first'))
keys_DMR <- data.frame(mapIds(org.Hs.eg.db, DMR_genelist,  keytype = 'ENTREZID', column = 'SYMBOL', multiVals = 'first'))
late_genes2 <-intersect(DMR_genes$annot.gene_id, keys$mapIds.org.Hs.eg.db..late_genes..keytype....SYMBOL...column....ENTREZID...)
#number of total genes = ~ 20,000
t_genes <- 20000
fish <- data.frame(c(length(late_genes2), length(late_genes) - length(late_genes2)), c(nrow(DMR_genes) - length(late_genes2), t_genes))
fisher.test(fish)
cpm <- cpm(read.csv("./both_counts.txt", sep = "\t", row.names = 1))
earlylate <- read.csv("earlylate replicated genes.csv")

early_rows <- intersect(rownames(cpm,), earlylate$gene.EARLY[earlylate$gene.EARLY != ""])
late_rows <- intersect(rownames(cpm), earlylate$gene.LATE[earlylate$gene.LATE != ""])
for_plot <- data.frame()
ggplot()
ggplot()

