library(ggbiplot)
library(factoextra)
library(tidyverse)  # data manipulation
library(cluster) 
cibersortx <- read.csv("../cibersortx_data.csv", row.names = "X")
ciber.pca <- prcomp(cibersortx[1:10])
g <- ggbiplot(ciber.pca, groups = as.factor(cibersortx$e.bin))
timer2 <- read.csv("../timer2.csv", row.names = 1)
timer2.pca <- prcomp(timer2[1:8])
gt <- ggbiplot(timer2.pca, groups = as.factor(timer2$e.bin))

fviz_nbclust(cibersortx[1:10], kmeans, method = "wss")

ciber.k <- kmeans(cibersortx[1:10], centers = 2  )
fviz_cluster(ciber.k, data = cibersortx[1:10])
