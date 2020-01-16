library(ggplot2)
library(cowplot)

setwd("/Users/imandric/workspace/projects/low-coverage/natcomm-revision")

essfile <- "ALL_CELLS.CT_DMX.csv"

ess <- read.csv(essfile, header=T)
rnames <- paste(ess$celltype, ess$cells, sep = "-")
rnames <- paste(rnames, ess$money, sep = "-")
rnames <- paste(rnames, ess$individuals, sep = "-")
rownames(ess) <- rnames

eqtlfile <- "eqtl.ALL.35000.txt"
table <- read.csv(eqtlfile, header=T, sep=" ")
rnames <- paste(table$celltype, table$cells, sep = "-")
rnames <- paste(rnames, "35000", sep = "-")
rnames <- paste(rnames, table$indiv, sep = "-")
rownames(table) <- rnames

table$ess <- ess[rownames(table),]$ess

table <- table[table$celltype == "cd4",]

table$Coverage <- ifelse(table$reads < 50000, "< 50,000", "> 50,000")

table$indiv2 <- factor(paste0(table$indiv, " samples"))
levels(table$indiv2) <- c("40 samples", "48 samples", "56 samples", "64 samples",
                          "72 samples", "80 samples", "88 samples", "96 samples",
                          "104 samples", "112 samples", "120 samples")


p11 <- ggplot(table, aes(x=sens, y=ppv)) + geom_point(size = 2) +
  geom_smooth(method='lm', formula=y~x, se = F) + facet_wrap(~indiv2, nrow = 4) +
  xlab("Recall (power estimate)") +
  ylab("Precision") +
  theme(strip.text.x = element_text(size = 15, colour = "black")) +
  theme(plot.title = element_text(size = 15), axis.text = element_text(size=20),
        axis.title = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        panel.grid.minor = element_line(colour = "grey50"))



p33 <- ggplot(table, aes(x=ess, y=sens)) + geom_point(size = 2) +
  geom_smooth(method='lm', formula=y~x, se = F) + facet_wrap(~indiv2, nrow = 4) +
  xlab("Recall (power estimate)") +
  ylab("Effective sample size") +
  theme(strip.text.x = element_text(size = 15, colour = "black")) +
  theme(plot.title = element_text(size = 15), axis.text = element_text(size=20),
        axis.title = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        panel.grid.minor = element_line(colour = "grey50"))



p22 <- ggplot(table, aes(x=reads / 1000, y=sens, colour = indiv, shape = Coverage)) + 
  geom_point(size = 3) +
  guides(colour = guide_legend(override.aes = list(size=5)), 
         shape = guide_legend(override.aes = list(size=5))) +
  ylab("Recall (power estimate)") +
  xlab(expression("Coverage (10"^3~"reads per cell)")) +
  scale_colour_gradient(low = "blue4", high = "orange", name = "Sample size") + 
  theme(plot.title = element_text(size = 15, face = "bold"), axis.text = element_text(size=20),
        axis.title = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.position = "none")



plot_grid(p33, p11, scale=0.9, labels=c("A", "B"), label_size = 30, align = "v", ncol = 1)

p22



# plot common with variance eqtl plot (see eqtl-plot-variance.R)
prow <- plot_grid(p22, p222 + theme(legend.position = "none"), scale=0.9, labels=c("A", "B"), label_size = 30, ncol = 2)
legend <- get_legend(
  # create some space to the left of the legend
  p222 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plot_grid(prow, legend, rel_widths = c(3, 0.4))



