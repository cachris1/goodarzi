gos <- read.delim("./file573847206.txt")
bp <- gos[gos["ConceptType"] == "GOBP",]
bp_sig <- bp[bp["FDR"] < 0.05,]
sperms <- filter(gos, grepl("sperm", Name))
immune <- filter(gos, grepl("immune", Name))
virus <- filter(gos, grepl("vir", Name))
print(nrow(bp_sig))
