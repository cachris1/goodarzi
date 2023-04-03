cpms <- read.csv("~/Bioinformatics/sartor/sartor/Integration/both_counts_logcpm.csv")
timer2 <- read.csv("../timer2.csv", row.names = 1)
timer2 <- timer2[timer2['HPVint'] == 'pos',]
cpms <- data_frame(cpms)
cpms <- t(cpms)
colnames(cpms) <- cpms[1,]
cpms <- cpms[-1,]
cpms <- cpms[rownames(timer2),]
cpms["DMBT1"] <- as.numeric(cpms["DMBT1"])
cpms <- cpms["DMBT1"]
timer2 <- merge(timer2, cpms, by = 0, all.y = TRUE)
boxplot(e.bin ~ DMBT1, timer2)