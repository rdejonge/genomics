##
## analyse pidents 
##

library(lattice)
library(ggplot2)
library(plotly)
library(reshape2)

setwd("C:/Users/Ronnie de Jonge/Documents/MEGAsync/Sync/work/current/Cercospora/CTB_Aurofusarin/")

df <- read.table("cbet_best_allspeciesheadersX", header = T, fill = F)
df_col_means <- colMeans(df, na.rm = T)
df <- df[,order(df_col_means)]

df_2 <- read.table("cbet_best_allspeciesheaders_ctbX", header = T, fill = F)
df_2 <- df_2[,order(df_col_means)]

df_3 <- read.table("cbet_best_allspeciesheaders_core", header = T, fill = F)
df_3 <- df_3[,order(df_col_means)]

    
boxplot(df_3, las = 2)
stripchart(df_2, vertical = T, data = df_2, method = "jitter", add = T, pch = 20, col = "red")
stripchart(df_3, vertical = T, data = df_3, method = "jitter", add = T, pch = 20, col = "blue")

