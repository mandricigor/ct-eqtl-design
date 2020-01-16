
library(reshape2)

setwd("/Users/imandric/workspace/projects/low-coverage/natcomm-revision")

# cell counts
#46156 B cells
#98623 CD14+ Monocytes
#124827 CD4 T cells
#40579 CD8 T cells
#3578 Dendritic cells
#16670 FCGR3A+ Monocytes
#7666 Megakaryocytes
#21901 NK cells

# dendritic -> megakaryocytes -> fcgr3a -> nk -> cd8 -> b -> cd14 -> cd4

labels <- c("Dendritic cells", "Megakaryocytes", "FCGR3A+ Monocytes", "NK cells",
            "CD8 T cells", "B cells", "CD14+ Monocytes", "CD4 T cells")



resultfile <- "ALL_CELLS.CT_DMX.csv"
results <- read.csv("ALL_CELLS.CT_DMX.csv")

celltypes <- c("dendritic", "megakaryocytes", "fcgr3a", "nk", "cd8", "b", "cd14", "cd4")

get_best_coverage <- function(ct, dollars) {
  tt <- results[results$celltype == ct & results$money == dollars,]
  mean(tail(tt[order(tt$ess),]$reads, n = 3)) / 1000  
}

k15 <- sapply(celltypes, get_best_coverage, 15000)
k20 <- sapply(celltypes, get_best_coverage, 20000)
k25 <- sapply(celltypes, get_best_coverage, 25000)
k30 <- sapply(celltypes, get_best_coverage, 30000)
k35 <- sapply(celltypes, get_best_coverage, 35000)
k40 <- sapply(celltypes, get_best_coverage, 40000)
k45 <- sapply(celltypes, get_best_coverage, 45000)
k50 <- sapply(celltypes, get_best_coverage, 50000)

kk <- data.frame(k15, k20, k25, k30, k35, k40, k45, k50)

colnames(kk) <- c(15, 20, 25, 30, 35, 40, 45, 50)

kk$celltype <- rownames(kk)
rownames(kk) <- NULL

kkmelted <- melt(kk, id.vars = "celltype")

p_cov <- ggplot(data=kkmelted, aes(x=variable, y=value, group=celltype, colour = celltype)) + 
  geom_line(linetype = "solid", size = 1) +
  geom_point(size = 2.5) +
  ylim(0, 17) +
  ylab(expression("Coverage (10"^3~"reads per cell)")) +
  xlab("Budget (thousand dollars)") +
  scale_color_discrete(name = "Cell type", labels = labels) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=30),
        plot.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20)) +
  geom_segment(aes(x=0.75,xend=8.5,y=10,yend=10), colour = "black", size = 1) +
  geom_segment(aes(x=0.75,xend=8.5,y=7.5,yend=7.5), colour = "black", linetype = "dashed", size = 1) +
  geom_segment(aes(x=0.75,xend=8.5,y=12.5,yend=12.5), colour = "black", linetype = "dashed", size = 1)
  
p_cov


onecase1 <- results[(results$individuals == 96 & results$cells == 2000 & results$money == 35000),]
onecase2 <- results[(results$individuals == 80 & results$cells == 2000 & results$money == 35000),]
onecase3 <- results[(results$individuals == 64 & results$cells == 2000 & results$money == 35000),]
onecase4 <- results[(results$individuals == 48 & results$cells == 2000 & results$money == 35000),]

rownames(onecase1) <- onecase1$celltype
onecase1 <- onecase1[celltypes,]
onecase1$labels <- labels
onecase1$experiment <- "A"

rownames(onecase2) <- onecase2$celltype
onecase2 <- onecase2[celltypes,]
onecase2$labels <- labels
onecase2$experiment <- "B"

rownames(onecase3) <- onecase3$celltype
onecase3 <- onecase3[celltypes,]
onecase3$labels <- labels
onecase3$experiment <- "C"

rownames(onecase4) <- onecase4$celltype
onecase4 <- onecase4[celltypes,]
onecase4$labels <- labels
onecase4$experiment <- "D"

onecase <- rbind(onecase1[c("labels", "experiment", "ess")], 
                 onecase2[c("labels", "experiment", "ess")],
                 onecase3[c("labels", "experiment", "ess")],
                 onecase4[c("labels", "experiment", "ess")])

onecase$labels <- factor(onecase$labels, levels = labels)

p_ess <- ggplot(data = onecase, aes(x = labels, y = ess, group = experiment, colour = experiment)) + geom_line() +
  xlab("Cell types") + 
  geom_point(size = 3) +
  scale_color_discrete(name = "Sample size", labels = c("96", "80", "64", "48")) + 
  ylab("Effective sample size") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(axis.text = element_text(size=20, angle = 90),
        axis.title = element_text(size=30),
        plot.title = element_text(size=20),
        legend.text = element_text(size=25),
        legend.title = element_text(size=20))
p_ess

plot_grid(p_ess, p_cov, scale=0.9, labels=c("A", "B"), label_size = 30, align = "v", ncol = 1)





onecase1$labels <- factor(onecase1$labels, levels = labels)
p_trend <- ggplot(data = onecase1, aes(x = labels, y = ess, group = experiment)) + geom_line() +
  xlab("Cell types") + 
  geom_smooth(method='lm', formula=y~x, se = F) +
  geom_point(size = 3) +
  scale_color_discrete(name = "Sample size", labels = c("96", "80", "64", "48")) + 
  ylab("Effective sample size") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_text(x = 4, y = 53, label = expression("R"^"2"~" = 0.67"), size = 10) +
  theme(axis.text = element_text(size=20, angle = 90),
        axis.title = element_text(size=30),
        plot.title = element_text(size=20),
        legend.text = element_text(size=25),
        legend.title = element_text(size=20))
p_trend

