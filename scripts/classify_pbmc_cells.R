
suppressWarnings(suppressMessages(library(SparseM)))
library(reticulate)
library(Seurat)
library(Matrix)
library(stringr)



args <- commandArgs(trailingOnly=T)

inputdata <- args[1]
outputdata <- args[2]

NR_CHUNK <- 50

pbmc <- readRDS("pbmc.classified.2700.rds")

pbmc_genes <- rownames(pbmc@assays$RNA@counts)


np <- import("numpy")


npz_matrix <- paste0(inputdata, ".MATRIX.npz")
genes <- paste0(inputdata, ".GENES.csv")
observations <- paste0(inputdata, ".OBSERVATIONS.csv")


npz1 <- np$load(npz_matrix)


mat <- new("matrix.coo", ra = as.numeric(npz1$f[["data"]]), ja = as.integer(1 + as.integer(npz1$f[["col"]])), ia = as.integer(1 + as.integer(npz1$f[["row"]])), dimension = as.integer(npz1$f[["shape"]]))
mat <- as(mat, "dgCMatrix")


observations <- read.csv(observations, stringsAsFactors=F)
cells <- observations$row
genes <- read.csv(genes, stringsAsFactors=F)$genes



rownames(mat) <- cells
colnames(mat) <- genes

mat <- t(mat)


igor <- CreateSeuratObject(counts = mat, project = "igor")
igor <- NormalizeData(igor)
igor <- FindVariableFeatures(igor, selection.method = "vst", nfeatures = 2000)


individuals <- unique(as.character(observations$ind_cov))

nr_cells <- length(cells)
igor <- AddMetaData(igor, metadata = observations$ind_cov, col.name = "ind_cov")
igor <- AddMetaData(igor, metadata = sample(1:NR_CHUNK, size=nr_cells, replace=T), col.name="chunk")
igor <- AddMetaData(igor, metadata = "NOTYPE", col.name="celltype")

for (i in 1:NR_CHUNK) {
    igori <- subset(igor, chunk==i)
    igori_cells <- WhichCells(igor, expression=`chunk`==i)
    anchors <- FindTransferAnchors(reference = pbmc, query = igori, dims = 1:30)
    predictions <- TransferData(anchorset = anchors, refdata = pbmc@active.ident, dims = 1:30)
    igori <- AddMetaData(igori, metadata = predictions)
    cts <- igori@meta.data$predicted.id
    igor@meta.data[igori_cells,]$celltype <- cts
}


table(igor@meta.data$celltype)



# TO REMOVE!
saveRDS(igor, "igor.vafli.rds")
write.csv(igor@meta.data$celltype, file="SUKA.CT.csv", quote=F)
quit()


table(observations$ct_cov == igor@meta.data$celltype)


# finding cell type specific gene expression

celltypes <- unique(igor@meta.data$celltype)


ct_names <- list()
ct_names[["B cells"]] <- "b_cells"
ct_names[["CD14+ Monocytes"]] <- "cd14_cells"
ct_names[["CD4 T cells"]] <- "cd4_cells"
ct_names[["CD8 T cells"]] <- "cd8_cells"
ct_names[["Dendritic cells"]] <- "dendritic_cells"
ct_names[["FCGR3A+ Monocytes"]] <- "fcgr3a_cells"
ct_names[["Megakaryocytes"]] <- "megakaryocytes_cells"
ct_names[["NK cells"]] <- "nk_cells"




for (ct in celltypes) {
    ct_name <- ct_names[[ct]]
    ctexpr <- list()
    for (ind in individuals) {
        wcells <- (igor@meta.data$celltype==ct & igor@meta.data$ind_cov==ind)
        wcells_data <- igor@assays$RNA@counts[, wcells]
        if (typeof(wcells_data) == "S4") {
            ctsge <- rowSums(wcells_data)
        } else {
            ctsge <- wcells_data
        }
        ctexpr[[ind]] <- ctsge
    }
    ctexpr <- as.data.frame(ctexpr)
    ctexpr <- ctexpr * 1000000.0 / outer(rep(1, dim(ctexpr)[1]), 0.001 + colSums(ctexpr))
    ctexpr <- log10(1.0 + ctexpr)
    outputdata_ct <- str_replace(outputdata, "CELLTYPE", ct_name)
    write.csv(ctexpr, file=paste0(outputdata_ct), quote=F)
}



