lada.csv <- read.csv("./DE_genes_by_integration-status_fromLada.csv")
lada.csv <- lada.csv %>% drop_na()
row.names(lada.csv) <- lada.csv[,"GeneID"]
lada.slack <- readxl::read_xlsx("./DE_genes_by_integration-status_fromLada.xlsx") 
lada.slack <- lada.slack[lada.slack$GeneID != "NA",]
lada.slack <- lada.slack %>% tibble::column_to_rownames("GeneID")
lad.slack <- data.frame(lada.slack)
#me.qlf <- diffexpr(all_int = TRUE, all_bin = TRUE, int = TRUE, ebin = FALSE, CCR = FALSE)
me.tags <- data.frame(topTags(me.qlf, n = nrow(me.qlf)))
me.tags <- merge(me.tags, o, by = 0)
me.tags <- me.tags %>% drop_na()
row.names(me.tags) <- me.tags$ENTREZ_ID
both <- merge(lad.slack, me.tags, by = 0, suffixes = c("lada.slack", "me"))