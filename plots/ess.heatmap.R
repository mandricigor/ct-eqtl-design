library(gplots)


args <- commandArgs(trailingOnly=T)

celltype <- args[1]
cost <- args[2]


files <- list.files()[grepl(paste0(celltype, "_", cost), list.files())]


r2 <- function(x) read.table(x, header=F, stringsAsFactors=F, row.names=NULL)
r2s <- sapply(files, r2)

colnames(r2s) <- NULL
rownames(r2s) <- c("celltype", "cost", "nr_individuals", "nr_cells", "nr_reads", "r2", "ess")
r2s <- as.data.frame(t(r2s))


#ss <- as.integer(sapply(strsplit(names(r2s), split="_"), function(x) x[4]))



#cols_ <- as.integer(sapply(strsplit(names(r2s), split="_"), function(x) x[5])) + as.integer(sapply(strsplit(names(r2s), split="_"), function(x) x[6]))



#r2s$nr_cells[r2s$nr_cells == 1950] <- 2000
#r2s$nr_cells[r2s$nr_cells == 2200] <- 2250
#r2s$nr_cells[r2s$nr_cells == 2450] <- 2500
#r2s$nr_cells[r2s$nr_cells == 2700] <- 2750

cols_ <- as.integer(r2s$nr_cells) / 250
cols <- as.integer(cols_) - 1


#cols_ <- round(cols_ / 50) * 50
#cols_ <- cols_ / 250
#cols <- cols_ - 1


#rows <- as.integer(sapply(strsplit(names(r2s), split="_"), function(x) x[4])) / 8 - 4


rows <- as.integer(r2s$nr_individuals) / 8 - 4




asmat <- data.frame(rows = rows, cols = cols, data = as.integer(r2s$ess))






mymat <- matrix(rep(0, 110), byrow=T, nrow=11)

for (i in 1:110) {mymat[asmat[i,]$rows, asmat[i,]$cols] = asmat[i,]$data}

colorfun <- function(x) {
    my_palette <- colorRampPalette(c("white", "yellow", "orange", "red"))(n = 10)
    my_palette[floor(x / 10) + 1]
}

mycols <- apply(mymat, 1:2, colorfun)

pdf(sprintf("heatmap_ess.%s_cells.cost_%s.pdf", celltype, cost))
heatmap.2(400 - mymat, key=F, trace="none", dendrogram="none", Rowv=F, Colv=F, labCol=c(500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750), labRow=c(40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120), cellnote=as.matrix(mymat), srtCol=45, notecol="black", ylab="# individuals", xlab="# cells", hline=T)
dev.off()



