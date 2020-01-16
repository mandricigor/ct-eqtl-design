
library(ggplot2)
setwd("/Users/imandric/workspace/projects/low-coverage/natcomm-revision")


esstable <- read.csv("ALL_CELLS.CT_DMX.csv")

esstable <- esstable[esstable$celltype == "dendritic",]
esstable$cost <- paste0("$", format(as.integer(esstable$money), nsmall=1, big.mark=","))


p11 <- ggplot(esstable, aes(x=reads / 1000, y=ess)) + geom_point(size = 1) + 
  facet_wrap(~cost, nrow = 4, scales = "free") +
  xlab("Coverage (thousand reads per cell)") +
  ylab("Effective sample size") +
  ylim(0, 90) +
  geom_vline(xintercept = 10, linetype = "dashed", colour = "blue") +
  geom_rect(aes(xmin=7.5, xmax=12.5, ymin=0, ymax=Inf), alpha = 0.005, colour = "red") +
  theme(strip.text.x = element_text(size = 15, colour = "black", face="bold")) +
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
p11


optimal_esstable <- esstable[esstable$reads <= 12500 & esstable$reads >= 7500,]
best_mean_ess <- aggregate(optimal_esstable$ess, list(optimal_esstable$money), mean)
colnames(best_mean_ess) <- c("cost", "ess")
best_mean_ess$cost <- best_mean_ess$cost / 1000

plot(best_mean_ess$cost, best_mean_ess$ess)

lm(formula = best_mean_ess$ess ~ best_mean_ess$cost)

