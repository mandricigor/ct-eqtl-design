
library(ggpubr)

setwd("/Users/imandric/workspace/projects/low-coverage/natcomm-revision")

aa <- read.csv("ALL_CELLS.CT_DMX.csv", header=T)


#colnames(aa) <- c("ss", "nc", "cov", "ess", "budget")
#aa$nc <- round(aa$nc / 50) * 50

aa$ss_nc <- paste(as.character(aa$individuals), as.character(aa$cells), sep=",")
aa$reads <- aa$reads / 1000
aa$money <- aa$money / 1000

aa$experiment <- paste(as.character(aa$individuals), as.character(aa$cells), 
                       as.character(aa$reads), sep=",")



exp <- "88,2250"
par(mar=c(5,5,2,5))
dd <- aa[(aa$ss_nc == exp) & (aa$celltype == "cd4"),]
with(dd, plot(reads, money, type="b", col="red3",lwd=2.5, ylim=c(0, 50), lty=1,
    xlab=expression("Coverage (10"^3~"reads per cell)"), 
    cex.lab=2,
    col.axis = "red3",
    ylab=NA,
    axes=F,
    main="N = 88, M = 2250",
    cex.main=2.5
))
axis(side=1, cex.axis=2, lwd.ticks = 3)
axis(side=2, cex.axis=2, col.axis="red3", lwd.ticks = 3)
par(new=T)
par(new=T)
with(dd, plot(reads, ess, type="b", axes=F, xlab=NA, ylab=NA, lwd=4, pch=1, lty=2,
    ylim=c(0, 105),
    cex.lab=2, cex.axis=2,
    col="blue"
))
box(which = "plot", lty = "solid")
axis(side = 4, cex.axis=2, col.axis="blue4", lwd.ticks = 3)
#arrows(20,15,16.5,39.5, col="black", lwd=3)
#arrows(10,85,12,64, col="black", lwd=3)
#text(20, 10, "Effective sample size", cex=2, col = "blue4")
#text(9, 94, "Budget\n($1,000)", cex=2, col = "red3")
mtext("Effective sample size", side = 4, line = 3, cex = 2, col = "blue4")
#text(par("usr")[2]*1.11,mean(par("usr")[3:4]), "Effective sample size", 
#     srt = -90, xpd = TRUE, pos = 4, col = "blue4", cex = 2, adj = 0.5)
mtext("Budget (thousand dollars)", side = 2, line = 3, cex = 2, col = "red3")
text(x = -5, y = 115, labels = "B", xpd = NA, cex = 2)



concat.text<-function(x,y,txt,col) {
    thisx<-x
    for(txtstr in 1:length(txt)) {
        text(thisx,y,txt[txtstr],col=col[txtstr],adj=0)
        thisx<-thisx+strwidth(txt[txtstr])
    }
}



# figure 2 from the paper - this is the main message
important <- aa[(aa$ess == 40 & aa$individuals == 88 & aa$cells == 2250 & aa$reads == 4.5) | 
    (aa$ess == 41 & aa$individuals == 80 & aa$cells == 1250 & aa$reads == 28) | 
    (aa$ess == 42 & aa$individuals == 72 & aa$cells == 1500 & aa$reads == 31) | 
    (aa$ess == 42 & aa$individuals == 64 & aa$cells == 1750 & aa$reads == 35.5) | 
    (aa$ess == 41 & aa$individuals == 56 & aa$cells == 2000 & aa$reads == 47.5),]
important <- important[important$celltype == "cd4",]
important$name <- c("LOW-COVERAGE (~5,000 rpc)", "~28,000 rpc", "~30,000 rpc", "~35,000 rpc", "STANDARD (~50,000 rpc)")


p <- ggbarplot(important, "name", "money", orientation = "vert", fill="green", xlab="Experimental design", ylab="Budget (thousand dollars)", order=c(rev(important$name)))
p + scale_x_discrete(labels=function(x) sub(" ","\n",x,fixed=TRUE)) + 
    ggtitle(expression("Effective sample size"~N[eff]%~~%"40")) +
    theme(axis.title = element_text(size=22, face="bold"),
        axis.text.y = element_text(color="blue", size=25),
        axis.text.x = element_text(color="blue", size=11, face="bold"),
        plot.title = element_text(face="bold", hjust=0.5, size=25)) + 
    geom_text(label=c("N=56\nM=2250\nr=47500",
                      "N=64\nM=1250\nr=35500",
                      "N=72\nM=1500\nr=31000",
                      "N=80\nM=1750\nr=28000",
                      "N=88\nM=2250\nr=4500"), 
                      size=4, position=position_stack(vjust=0.5)) + 
    theme(plot.subtitle=element_text(hjust=0.5, face="italic", color="black"))
#ggsave("figure2.png", width = 8, height = 5, dpi=300)
#p


# THEN ADD THE LABELS AND MERGE TWO PICTURES ONLINE!!!!!


