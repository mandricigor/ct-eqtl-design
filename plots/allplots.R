
# Figure S9 B)

library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)


setwd("/Users/imandric/workspace/projects/low-coverage/natcomm-revision")

gexp_ <- "final_data_24K/final_data/final_data.down_100.ct_0.csv"
cors_ <- "final_data_24K/r2s/r2.down_90.ct_0.csv"
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


# Figure S9 A)

dataStratify <- data.frame(x=1:8, y=as.vector(table(cors2$nonzero) / 10000))
p2 <- ggplot(dataStratify, aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  ylab(expression("# genes, 10"^"3")) +
  scale_x_continuous("# individuals", labels = 1:8, breaks = 1:8) +
  theme(legend.position = "none",
        axis.text = element_text(size=30),
        axis.title = element_text(size=30),
        plot.title = element_text(size=30))




# Figure 1 B)

expr2 <- function(x) {
  print (length(cors2[cors2$means > 0.5 * x & cors2$means < 0.5 * (x + 1),]$r2))
  mean(cors2[cors2$means > 0.5 * x & cors2$means < 0.5 * (x + 1),]$r2)
}

e2 <- sapply(1:15, expr2)

ee2 <- data.frame(x=1:15, y=e2)
p3 <- ggplot(ee2, aes(x=x, y=y)) +
  geom_line(size=2) +
  geom_point(shape=21, fill="white", size=4, stroke=2) +
  scale_x_continuous("Mean expression", labels=seq(0.5, 7.5, 1), breaks=seq(1, 15, 2)) +
  ylab(expression("Peasron R"^"2")) +
  ylim(0, 1) +
  theme(legend.position = "none",
        axis.text.y = element_text(size=30),
        axis.text.x = element_text(size=20),
        axis.title = element_text(size=30),
        plot.title = element_text(size=30))



# Figure 2

aa <- read.csv("B_CELLS.CT.csv", header=T)
#a <- read.csv("/home/igorm/workspace/PAPER_SINGLE_CELL_PLOTS/summaries2/B_CELLS", header=F)

colnames(aa) <- c("ss", "nc", "cov", "ess", "budget")
aa$nc <- round(aa$nc / 50) * 50
aa$ss_nc <- paste(as.character(aa$ss), as.character(aa$nc), sep=",")
aa$experiment <- paste(as.character(aa$ss), as.character(aa$nc), as.character(aa$cov), sep=",")




important <- aa[(aa$ess == 45 & aa$ss == 56 & aa$nc == 2750 & aa$cov == 33) |
                 (aa$ess == 45 & aa$ss == 96 & aa$nc == 1250 & aa$cov == 39) |
                 (aa$ess == 45 & aa$ss == 120 & aa$nc == 1000 & aa$cov == 25) |
                 (aa$ess == 45 & aa$ss == 64 & aa$nc == 2250 & aa$cov == 16) |
                 (aa$ess == 45 & aa$ss == 96 & aa$nc == 2500 & aa$cov == 1.5),]
important$name <- c("BEST DESIGN", "LARGE # CELLS", "LARGE SAMPLE SIZE", "HIGH COVERAGE", "RECOMMENDED DESIGN")


p <- ggbarplot(important, "name", "budget", orientation = "vert", fill="green", xlab="Experimental design", ylab="Budget ($1000)", order=c(rev(important$name)))
p + scale_x_discrete(labels=function(x) sub(" ","\n",x,fixed=TRUE)) +
  ggtitle(expression("Effective sample size"~N[eff]%~~%"45")) +
  #ggtitle("Effective sample size\n- Known cell types: N=50\n- Label transfer: N=45") +
  theme(axis.title = element_text(size=18, face="bold"),
        axis.text.y = element_text(color="blue", size=25),
        axis.text.x = element_text(color="blue", size=10),
        plot.title = element_text(face="bold", hjust=0.5, size=25)) +
  geom_text(label=c("N=56\nM=2750\nr=33000",
                    "N=96\nM=1250\nr=39000",
                    "N=120\nM=1000\nr=25000",
                    "N=64\nM=2250\nr=16000",
                    "N=96\nM=2500\nr=1500"),
            size=4, position=position_stack(vjust=0.5)) +
  theme(plot.subtitle=element_text(hjust=0.5, face="italic", color="black"))


