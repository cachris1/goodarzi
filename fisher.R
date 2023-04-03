timer2 <- read.csv("../timer2.csv", row.names = 1)
#timer2 <- timer2[timer2$HPVint == "pos",]
fish <- data.frame(table(timer2[timer2$e.bin == 1,]$Cohort))
fish[,3] <- (table(timer2[timer2$e.bin == 0,]$Cohort))
fish <- fish[-1]