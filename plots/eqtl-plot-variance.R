

setwd("/Users/imandric/workspace/projects/low-coverage/natcomm-revision")

#essfile <- "ALL_CELLS.CT_DMX.csv"

#ess <- read.csv(essfile, header=T)
#rnames <- paste(ess$celltype, ess$cells, sep = "-")
#rnames <- paste(rnames, ess$money, sep = "-")
#rnames <- paste(rnames, ess$individuals, sep = "-")
#rownames(ess) <- rnames

eqtlfile <- "ALL.varEQTL.cd14.35000.txt"
table <- read.csv(eqtlfile, header=T, sep=" ")
rnames <- paste(table$celltype, table$cells, sep = "-")
rnames <- paste(rnames, "35000", sep = "-")
rnames <- paste(rnames, table$indiv, sep = "-")
rownames(table) <- rnames

#table$ess <- ess[rownames(table),]$ess

#table <- table[table$celltype == "cd4",]

table$Coverage <- ifelse(table$reads < 50000, "< 50,000", "> 50,000")








p222 <- ggplot(table, aes(x=reads / 1000, y=sens, colour = indiv, shape = Coverage)) + 
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
        legend.title = element_text(size=15))



p222

