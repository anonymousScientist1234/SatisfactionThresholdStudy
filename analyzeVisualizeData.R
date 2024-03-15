#Clear Work Space
rm(list=ls())

library(tidyverse)
library(ggpubr)
library(ggthemes)
library(EnvStats)
library(FSA)
library(forcats)
library(ggsignif)
library(readr)

setwd("~/Satisfaction Threshold/Simple Model/final model/from laptop/final final final")

data = read.csv("thresholdComparisons.csv")

f <- function(x) {
  r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#Data analysis. Do kruskal wallis tests, dunn tests, and correlation tests, showing what factors influence performance metrics 

#First we look into DOL

dataWilcox = data[order(data$CueType),]
dataTD = subset(data, CueType == "Task Demand")
dataTC = subset(data, CueType == "Task Completion")

kruskal.test(DOL~ThresholdType, data=data)
dunnTest(DOL~ThresholdType, data=data)
wilcox.test(DOL ~ CueType, data = dataWilcox)
cor.test(as.numeric(data$T), data$DOL, method="spearman", na.omit=TRUE)
aggregate(DOL~ThresholdType, data=data, FUN = median)
aggregate(DOL~CueType, data=data, FUN = median)
aggregate(DOL~T, data=data, FUN = median)

kruskal.test(DOL~ThresholdType, data=dataTD)
dunnTest(DOL~ThresholdType, data=dataTD)
kruskal.test(DOL~ThresholdType, data=dataTC)
dunnTest(DOL~ThresholdType, data=dataTC)

#Task switching

kruskal.test(TaskSwitches~ThresholdType, data=data)
dunnTest(TaskSwitches~ThresholdType, data=data)
wilcox.test(TaskSwitches ~ CueType, data = dataWilcox, paired = TRUE)
cor.test(as.numeric(data$T), data$TaskSwitches, method="spearman")
aggregate(TaskSwitches~ThresholdType, data=data, FUN = median)
aggregate(TaskSwitches~CueType, data=data, FUN = median)
aggregate(TaskSwitches~T, data=data, FUN = median)

kruskal.test(TaskSwitches~ThresholdType, data=dataTD)
dunnTest(TaskSwitches~ThresholdType, data=dataTD)
kruskal.test(TaskSwitches~ThresholdType, data=dataTC)
dunnTest(TaskSwitches~ThresholdType, data=dataTC)

#Loading

kruskal.test(Loading~ThresholdType, data=data)
dunnTest(Loading~ThresholdType, data=data)
wilcox.test(Loading ~ CueType, data = dataWilcox, paired = TRUE)
cor.test(as.numeric(data$T), data$Loading, method="spearman")
aggregate(Loading~ThresholdType, data=data, FUN = median)
aggregate(Loading~CueType, data=data, FUN = median)
aggregate(Loading~T, data=data, FUN = median)

kruskal.test(Loading~ThresholdType, data=dataTD)
dunnTest(Loading~ThresholdType, data=dataTD)
kruskal.test(Loading~ThresholdType, data=dataTC)
dunnTest(Loading~ThresholdType, data=dataTC)

#Timesteps to EQ

kruskal.test(TimetoEq~ThresholdType, data=data)
dunnTest(TimetoEq~ThresholdType, data=data)
wilcox.test(TimetoEq ~ CueType, data = dataWilcox, paired = TRUE)
cor.test(as.numeric(data$T), data$TimetoEq, method="spearman")
aggregate(TimetoEq~ThresholdType, data=data, FUN = median)
aggregate(TimetoEq~CueType, data=data, FUN = median)
aggregate(TimetoEq~T, data=data, FUN = median)

kruskal.test(TimetoEq~ThresholdType, data=dataTD)
dunnTest(TimetoEq~ThresholdType, data=dataTD)
kruskal.test(TimetoEq~ThresholdType, data=dataTC)
dunnTest(TimetoEq~ThresholdType, data=dataTC)

### Figure 1

cols = c('#009E73','#CC79A7')

dataWilcox$ThresholdType = fct_relevel(dataWilcox$ThresholdType, "Response Threshold", "Satisfaction Threshold", "Composite Threshold", "Random Choice")

data$ThresholdType = fct_relevel(data$ThresholdType, "Response Threshold", "Satisfaction Threshold", "Composite Threshold", "Random Choice")

anno_df = compare_means(DOL ~ CueType, group.by = "ThresholdType", data = dataWilcox, paired = TRUE) %>%
  mutate(y_pos = .95)
anno_df$p.adj.char = paste("p =", anno_df$p.adj)
dataWilcox$DOL[dataWilcox$DOL>1] = NA

p1 = ggplot(dataWilcox, aes(x = CueType, y = DOL))  + facet_wrap(~ThresholdType) + geom_jitter(aes(color = CueType), width = .15, alpha = .8) + stat_summary(fun.data=f, geom="boxplot", aes(color = CueType), alpha = .25, size = 1.2)+ theme_bw() + guides(color="none") + scale_color_manual(values=cols) + ylab("Division of Labor Index") + xlab("Cue Type") + theme(text = element_text(size=14)) + scale_x_discrete(guide = guide_axis(n.dodge = 2)) + scale_y_continuous(labels = function(x) {
  format(x, nsmall = 2)}, limits = c(0,1.1))  + ggtitle("A)") +
  ggsignif::geom_signif(
    data=anno_df,
    aes(xmin=group1, xmax=group2, annotations=p.adj.char, y_position=y_pos),
    manual=TRUE, size = 0)

anno_df = compare_means(TaskSwitches ~ CueType, group.by = "ThresholdType", data = dataWilcox, paired = TRUE) %>% mutate(y_pos = 5200000)
anno_df$p.adj.char = paste("p =", anno_df$p.adj)

p2 = ggplot(dataWilcox, aes(x = CueType, y = TaskSwitches))  + facet_wrap(~ThresholdType) + geom_jitter(aes(color = CueType), width = .15, alpha = .8) + stat_summary(fun.data=f, geom="boxplot", aes(color = CueType), alpha = .25, size = 1.2)+ theme_bw() + guides(color="none") + scale_color_manual(values=cols) + ylab(expression("#"~Task~Switches~(x10^4))) + xlab("Cue Type") + theme(text = element_text(size=14)) + scale_x_discrete(guide = guide_axis(n.dodge = 2))+ scale_y_continuous(labels = function(x) x/10000, limits = c(0,6000000)) + ggtitle("D)") +
  ggsignif::geom_signif(
    data=anno_df,
    aes(xmin=group1, xmax=group2, annotations=p.adj.char, y_position=y_pos),
    manual=TRUE, size = 0)

anno_df = compare_means(Loading ~ CueType, group.by = "ThresholdType", data = dataWilcox, paired = TRUE) %>% mutate(y_pos = 3100000)
anno_df$p.adj.char = paste("p =", anno_df$p.adj)

p3 = ggplot(dataWilcox, aes(x = CueType, y = Loading))  + facet_wrap(~ThresholdType) + geom_jitter(aes(color = CueType), width = .15, alpha = .8) + stat_summary(fun.data=f, geom="boxplot", aes(color = CueType), alpha = .25, size = 1.2)+ theme_bw() + guides(color="none") + scale_color_manual(values=cols) + ylab(expression(Loading~(Worker-Timesteps~x10^4)))+ xlab("Cue Type") + theme(text = element_text(size=14)) + scale_x_discrete(guide = guide_axis(n.dodge = 2))+ scale_y_continuous(labels = function(x) x/10000, limits = c(-5000000, 4000000)) + ggtitle("B)") +
  ggsignif::geom_signif(
    data=anno_df,
    aes(xmin=group1, xmax=group2, annotations=p.adj.char, y_position=y_pos),
    manual=TRUE, size = 0) + geom_hline(yintercept = 0, size = 1.2, color = "darkgray", linetype="dashed")

anno_df = compare_means(TimetoEq ~ CueType, group.by = "ThresholdType", data = dataWilcox, paired = TRUE) %>% mutate(y_pos = 5200)
anno_df$p.adj.char = paste("p =", round(anno_df$p, 3))

p4 = ggplot(dataWilcox, aes(x = CueType, y = TimetoEq))  + facet_wrap(~ThresholdType) + geom_jitter(aes(color = CueType), width = .15, alpha = .8) + stat_summary(fun.data=f, geom="boxplot", aes(color = CueType), alpha = .25, size = 1.2)+ theme_bw() + guides(color="none") + scale_color_manual(values=cols) + ylab(expression(Timesteps~to~Equilibrium~(x10^1))) + xlab("Cue Type") + theme(text = element_text(size=14)) + scale_x_discrete(guide = guide_axis(n.dodge = 2))+ scale_y_continuous(labels = function(x) x/10, limits = c(2,6000)) + ggtitle("C)") +
  ggsignif::geom_signif(
    data=anno_df,
    aes(xmin=group1, xmax=group2, annotations=p.adj.char, y_position=y_pos),
    manual=TRUE, size = 0)

ggarrange(p1, p3, p4, p2, nrow= 2, ncol = 2)

#figure 2: colony size vs performance

cols = c("Response Threshold" = '#66ccee', "Satisfaction Threshold" = '#009988', "Composite Threshold" = "#cc3311", "Random Choice" = "#bbbbbb")

p1 = ggplot(data, aes(x = N, y = DOL, color = ThresholdType)) + geom_point(size = .25) + geom_smooth(size = 1.2, se = T) + theme_bw() + labs(color = "Decision Rule") + theme(text = element_text(size=14)) + ggtitle("A)") + scale_color_manual(values=cols) + scale_y_continuous(labels = function(x) {format(x, nsmall = 2)}, limits = c(0,1)) + ylab("Division of Labor Index")

p2 = ggplot(data, aes(x = N, y = TaskSwitches/N, color = ThresholdType)) + geom_point(size = .25) + geom_smooth(size = 1.2, se = T) + theme_bw() + labs(color = "Decision Rule") + theme(text = element_text(size=14)) + ggtitle("D)") + scale_color_manual(values=cols) + scale_y_continuous(labels = function(x) x/10, limits = c(0,5000)) + ylab(expression("#"~Task~Switches~per~Capita~(x10^1)))

p3 = ggplot(data, aes(x = N, y = Loading, color = ThresholdType)) + geom_point(size = .25) + geom_smooth(size = 1.2, se = T) + theme_bw() + labs(color = "Decision Rule") + theme(text = element_text(size=14)) + ggtitle("B)") + scale_color_manual(values=cols) + scale_y_continuous(labels = function(x) x/10000, limits = c(-4500000, 2500000)) + ylab(expression(Loading~(Worker-Timesteps~x10^4)))

p4 = ggplot(data, aes(x = N, y = TimetoEq, color = ThresholdType)) + geom_point(size = .25) + geom_smooth(size = 1.2, se = T) + theme_bw() + labs(color = "Decision Rule") + theme(text = element_text(size=14)) + ggtitle("C)") + scale_color_manual(values=cols) + scale_y_continuous(labels = function(x) x/10, limits = c(2,5000)) + ylab(expression(Timesteps~to~Equilibrium~(x10^1)))

ggarrange(p1, p3, p4, p2, nrow= 2, ncol = 2, common.legend = TRUE)

#figure 3: colony variance vs performance 

p1 = ggplot(data, aes(x = Sigma, y = DOL, color = ThresholdType)) + geom_point(size = .25) + geom_smooth(size = 1.2, se = T) + theme_bw() + labs(color = "Decision Rule") + theme(text = element_text(size=14)) + ggtitle("A)") + scale_color_manual(values=cols) + scale_y_continuous(labels = function(x) {format(x, nsmall = 2)}, limits = c(0,1)) + ylab("Division of Labor Index") + xlab(expression(sigma))

p2 = ggplot(data, aes(x = Sigma, y = TaskSwitches/N, color = ThresholdType)) + geom_point(size = .25) + geom_smooth(size = 1.2, se = T) + theme_bw() + labs(color = "Decision Rule") + theme(text = element_text(size=14)) + ggtitle("D)") + scale_color_manual(values=cols) + scale_y_continuous(labels = function(x) x/10, limits = c(0,5000)) + ylab(expression("#"~Task~Switches~per~Capita~(x10^1))) + xlab(expression(sigma))

p3 = ggplot(data, aes(x = Sigma, y = Loading, color = ThresholdType)) + geom_point(size = .25) + geom_smooth(size = 1.2, se = T) + theme_bw() + labs(color = "Decision Rule") + theme(text = element_text(size=14)) + ggtitle("B)") + scale_color_manual(values=cols) + scale_y_continuous(labels = function(x) x/10000, limits = c(-4500000, 2500000)) + ylab(expression(Loading~(Worker-Timesteps~x10^4))) + xlab(expression(sigma))

p4 = ggplot(data, aes(x = Sigma, y = TimetoEq, color = ThresholdType)) + geom_point(size = .25) + geom_smooth(size = 1.2, se = T) + theme_bw() + labs(color = "Decision Rule") + theme(text = element_text(size=14)) + ggtitle("C)") + scale_color_manual(values=cols) + scale_y_continuous(labels = function(x) x/10, limits = c(2,5000)) + ylab(expression(Timesteps~to~Equilibrium~(x10^1))) + xlab(expression(sigma))

ggarrange(p1, p3, p4, p2, nrow= 2, ncol = 2, common.legend = TRUE)

#figure 4: effect of task number

cols = c("1" = '#fc8961', "2" = '#b73779', "4" = "#51127c", "8" = "#000004")
data$T = as.factor(data$T)
data$DOL[data$DOL>1] = NA

p1 = ggplot(data, aes(x = as.numeric(T), y = DOL)) + facet_wrap(~ThresholdType) + geom_jitter(aes(color = T), width = .15, alpha = .8) + stat_summary(fun.data=f, geom="boxplot", aes(color = T), alpha = .25, size = 1.2) + theme_bw() + guides(color="none") + ylab("Division of Labor Index") + xlab("# Tasks (T)") + theme(text = element_text(size=14)) + ggtitle("A)") + stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p.adjust(readr::parse_number(..p.label..), n = 4))), method = "spearman") + scale_x_continuous(breaks = 1:4, labels = c(1, 2, 4, 8))  + scale_color_manual(values=cols)  + scale_y_continuous(labels = function(x) {format(x, nsmall = 2)}, limits = c(0,1.1))

p2 = ggplot(data, aes(x = as.numeric(T), y = TaskSwitches)) + facet_wrap(~ThresholdType) + geom_jitter(aes(color = T), width = .15, alpha = .8) + stat_summary(fun.data=f, geom="boxplot", aes(color = T), alpha = .25, size = 1.2) + theme_bw() + guides(color="none") + scale_color_manual(values=cols) + ylab(expression("#"~Task~Switches~(x10^4)))  + xlab("# Tasks (T)") + theme(text = element_text(size=14)) + ggtitle("D)") + scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 5100000))  + stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p.adjust(readr::parse_number(..p.label..), n = 4))), method = "spearman") + scale_x_continuous(breaks = 1:4, labels = c(1, 2, 4, 8)) + scale_y_continuous(labels = function(x) x/10000, limits = c(0,5550000))

p3 = ggplot(data, aes(x = as.numeric(T), y = Loading)) + facet_wrap(~ThresholdType) + geom_jitter(aes(color = T), width = .15, alpha = .8) + stat_summary(fun.data=f, geom="boxplot", aes(color = T), alpha = .25, size = 1.2) + theme_bw() + guides(color="none") + scale_color_manual(values=cols) + ylab(expression(Loading~(Worker-Timesteps~x10^4)))+ xlab("# Tasks (T)") + theme(text = element_text(size=14)) + ggtitle("B)") + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))  + stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p.adjust(readr::parse_number(..p.label..), n = 4))), method = "spearman") + scale_x_continuous(breaks = 1:4, labels = c(1, 2, 4, 8)) + scale_y_continuous(labels = function(x) x/10000, limits = c(-4500000, 3900000))

p4 = ggplot(data, aes(x = as.numeric(T), y = TimetoEq)) + facet_wrap(~ThresholdType) + geom_jitter(aes(color = T), width = .15, alpha = .8) + stat_summary(fun.data=f, geom="boxplot", aes(color = T), alpha = .25, size = 1.2) + theme_bw() + guides(color="none") + scale_color_manual(values=cols) + ylab(expression(Timesteps~to~Equilibrium~(x10^1))) + xlab("# Tasks (T)") + theme(text = element_text(size=14)) + ggtitle("C)") + scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0,6000))  + stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p.adjust(readr::parse_number(..p.label..), n = 4))), method = "spearman") + scale_x_continuous(breaks = 1:4, labels = c(1, 2, 4, 8)) + scale_y_continuous(labels = function(x) x/10, limits = c(2,6000))

ggarrange(p1, p3, p4, p2, nrow= 2, ncol = 2)

### Fig S1: Truncated Normal Distribution

x = seq(0,1, length.out = 1000)
y1 = dnormTrunc(x, mean = 0, sd = .1, min = 0, max = 1)
y2 = dnormTrunc(x, mean = .5, sd = .1, min = 0, max = 1)
y3 = dnormTrunc(x, mean = 1, sd = .1, min = 0, max = 1)

y4 = dnormTrunc(x, mean = 0, sd = .5, min = 0, max = 1)
y5 = dnormTrunc(x, mean = .5, sd = .5, min = 0, max = 1)
y6 = dnormTrunc(x, mean = 1, sd = .5, min = 0, max = 1)

y7 = dnormTrunc(x, mean = 0, sd = 1, min = 0, max = 1)
y8 = dnormTrunc(x, mean = .5, sd = 1, min = 0, max = 1)
y9 = dnormTrunc(x, mean = 1, sd = 1, min = 0, max = 1)

dataGraph = data.frame(x = x, y = c(y1, y2, y3, y4, y5, y6, y7, y8, y9), Mu = c(rep(0, 1000), rep(.5, 1000), rep(1, 1000), rep(0, 1000), rep(.5, 1000), rep(1, 1000), rep(0, 1000), rep(.5, 1000), rep(1, 1000)), Sigma = c(rep(.1, 3000), rep(.5, 3000),rep(1, 3000)))
dataGraph$Sigma = as.character(dataGraph$Sigma)

dataGraph$Sigma = factor(dataGraph$Sigma, levels = c("0.1", "0.5", "1"), labels = c(expression(paste(sigma, " = 0.1")), expression(paste(sigma, " = 0.5")), expression(paste(sigma, " = 1"))))

ggplot(dataGraph, aes(x = x, y = y, color = as.factor(Mu))) + geom_line(size=1.2) + facet_wrap(~Sigma, labeller = label_parsed) + theme_bw() + guides(color=guide_legend(title=expression(mu)))  + ylab("PDF") + xlab(expression(theta)) + theme(text = element_text(size=14)) + geom_hline(yintercept=1, linetype="dashed", color = "gray", size=1.2)

