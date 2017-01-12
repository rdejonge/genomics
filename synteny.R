###plotting synteny score data for clusters
###Ronnie de Jonge r.dejonge@uu.nl January 9, 2017

library(ggplot2)
library(Hmisc)
library(dplyr)
library(reshape)
library(lme4)
library(nlme)
library(labeling)

setwd("C:/Users/Ronnie de Jonge/Documents/MEGAsync/Sync/work/current/Cercospora/D_Cercospora_sync/Manuscript/Drafts/draft6/manuscript_/")

data <- read.table("synteny_scores.txt", header = T)
data_melted <- melt(data = data, measure.vars = c("Cg", "Mo"))
##prepare levels for sorting by chromosome size and apply
f=c('Chr1','Chr3','Chr2','Chr6','Chr9','Chr5','Chr7','Chr4','Chr8','Chr10')

data_melted$scaffold = factor(data_melted$scaffold, levels=f)

d2 <- data_melted %>%
  group_by(variable) %>%
  summarize(upper = log10(quantile(value, probs = .99)))

##renaming labeller
chromosomes <- c(
  `Chr1` = "1",
  `Chr3` = "2",
  `Chr2` = "3",
  `Chr6` = "4",
  `Chr9` = "5",
  `Chr5` = "6",
  `Chr7` = "7",
  `Chr4` = "8",
  `Chr8` = "9",
  `Chr10` = "10"
)

##plot the data by facetting
p <- ggplot(data = data_melted, aes(x=start_loc, y=log10(value))) +
  geom_point(aes(color=variable))
p + facet_grid(~scaffold, scales = "free", space = "free", labeller = as_labeller(chromosomes)) + 
  scale_x_continuous(breaks = seq(1000000, 7500000, by = 1000000), labels = NULL) +
  geom_hline(data = d2, aes(yintercept = upper, colour = variable), linetype="dashed") +
  labs(x = "")
