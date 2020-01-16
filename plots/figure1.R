

library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)


setwd("/Users/imandric/workspace/projects/low-coverage/natcomm-revision")
gexp_ <- "final_data_24K/final_data/final_data.down_100.ct_0.csv"
cors_ <- "final_data_24K/r2s/r2.down_90.ct_0.csv" # this is in fact 10 percent downsampling
gexp <- read.csv(gexp_)
rownames(gexp) <- gexp$X
gexp <- gexp[c("AZ","HP1502401","HP1504101T2D","HP1504901","HP1506401","HP1507101","HP1508501T2D","HP1509101","HP1525301T2D","HP1526901T2D")]
cors <- read.csv(cors_)
means <- rowMeans(gexp)
cors$means <- means
nonzero <- sapply(1:24181, function(x) length(gexp[x,][gexp[x,] > 0]))
cors$nonzero <- nonzero
cors2 <- cors[complete.cases(cors),]
cors2$nonzero <- factor(cors2$nonzero)
cormean <- function(x) {mean(cors2[cors2$nonzero == x,]$r2)}
corsd1 <- function(x) {quantile(cors2[cors2$nonzero == x,]$r2)[2]}
corsd2 <- function(x) {quantile(cors2[cors2$nonzero == x,]$r2)[4]}
cms <- sapply(1:8, cormean)
csd1 <- sapply(1:8, corsd1) 
csd2 <- sapply(1:8, corsd2) 
cdata <- data.frame(x=1:8, cs1=csd1, cs2=csd2, cms=cms)

p1 <- ggplot(cdata, aes(x=x, y=cms)) +
  geom_errorbar(aes(ymin=csd1, ymax=csd2), width=.2, size=1) + 
  geom_line(size=1) +
  geom_point(size=5) +
  ylab(expression("Pearson R"^"2")) +
  scale_x_continuous("# individuals", labels = 1:8, breaks = 1:8) +
  theme(legend.position = "none",
        axis.text = element_text(size=25),
        axis.title = element_text(size=25),
        plot.title = element_text(size=25))

p1




dataStratify <- data.frame(x=1:8, y=as.vector(table(cors2$nonzero) / 10000))
p2 <- ggplot(dataStratify, aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  ylab(expression("# genes, 10"^"3")) +
  scale_x_continuous("# individuals", labels = 1:8, breaks = 1:8) +
  theme(legend.position = "none",
        axis.text = element_text(size=30),
        axis.title = element_text(size=30),
        plot.title = element_text(size=30))
p2







expr2 <- function(x) {
  print (length(cors2[cors2$means > 0.5 * x & cors2$means < 0.5 * (x + 1),]$r2))
  mean(cors2[cors2$means > 0.5 * x & cors2$means < 0.5 * (x + 1),]$r2) 
}

expr3 <- function(x) {
  sd(cors2[cors2$means > 0.5 * x & cors2$means < 0.5 * (x + 1),]$r2) 
}

e2 <- sapply(1:15, expr2)
e3 <- sapply(1:15, expr3)
ee2 <- data.frame(x=1:15, y=e2)
ee3 <- data.frame(x=1:15, y=e3)
p3 <- ggplot(ee2, aes(x=x, y=y)) +
  geom_line(size=2) +
  geom_point(shape=21, fill="white", size=3, stroke=2) + 
  #scale_x_continuous("Mean expression", labels=seq(0.5, 7.5, 1), breaks=seq(1, 15, 2)) +
  scale_x_continuous("Mean expression", 
                     labels=sapply(seq(0.5, 7.5, 1), function(x) sprintf("%.1f-%s", x - 0.5, x)), 
                     breaks=seq(1, 15, 2)) +
  geom_ribbon(aes(ymin=ee2$y - ee3$y, ymax=pmin(rep(1, length(ee2$y)), ee2$y + ee3$y)), linetype=3, alpha=0.2) +
  ylab(expression("Pearson R"^"2")) +
  ylim(0, 1) +
  theme(legend.position = "none",
        axis.text = element_text(size=30),
        axis.text.x = element_text(size=13, angle = 45, hjust=1, face = "bold"),
        axis.title = element_text(size=30),
        plot.title = element_text(size=30))
p3








#setwd("/home/igorm/workspace/PAPER_REAL_PLOTS/smartseq/")
scaleFUN <- function(x) sprintf("%.2f", x)
results <- list()
sdres <-list()
for (ct in 0) {
  r2s <- numeric()
  s2 <- numeric()
  for (down in c(90, 80, 70, 60, 50, 40, 30, 20, 10, 100)) {
    table <- read.csv(sprintf("final_data_24K/r2s/r2.down_%s.ct_%s.csv", down, ct), header=T, row.names=1)
    meanr2 <- mean(table$r2, na.rm=T)
    r2s <- c(r2s, meanr2)
    ss2 <- sd(table$r2, na.rm=T)
    s2 <- c(s2, ss2)
  }
  results[[sprintf("celltype_%s", ct)]] <- r2s
  sdres[[sprintf("celltype_%s", ct)]] <- s2
}



results <- data.frame(results)
results$x <- seq(75, 750, 75)
p4 <- ggplot(results, aes(x=x)) +
  geom_line(aes(y = celltype_0), size=2) +
  geom_point(aes(y = celltype_0), shape=21, fill="white", size=3, stroke=2) +
ylim(0, 1) +
  ylab(expression("Pearson R"^"2")) +
  #scale_y_continuous(labels=scaleFUN) +
  scale_x_continuous(expression(paste("Coverage per cell,10"^"3", " reads")), labels=c(75, 300, 525, 750), breaks=c(75, 300, 525, 750)) +
  #scale_x_continuous(expression(paste("Coverage,10"^"3", " reads")), labels=results$x, breaks=results$x) +
  geom_ribbon(aes(ymin=results$celltype_0 - sdres$celltype_0, ymax=pmin(rep(1, length(results$celltype_0)), results$celltype_0 + sdres$celltype_0)), linetype=3, alpha=0.2) + 
  theme(axis.text.y = element_text(size=30),
        axis.text.x = element_text(size=30),
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        legend.position = "none")
p4





plot_grid(p4, p3, scale=0.9, labels=c("A", "B"), label_size = 30)

plot_grid(p2, p1, scale=0.9, labels=c("A", "B"), label_size = 30)

